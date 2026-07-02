# /// script
# dependencies = ["psycopg2-binary", "jinja2", "pyyaml"]
# ///
"""
Generate a genome note XML (JATS format) for an OceanGenomes sample.

Queries the OceanOmics PostgreSQL DB and renders the Jinja2 template at
genome_notes_automation/oceangenomes_genome_note.xml.j2

Usage:
  uv run --script generate_genome_note.py OG38
  uv run --script generate_genome_note.py OG38 ~/postgresql_details/oceanomics.cfg

Fields left as PLACEHOLDER_* require manual entry (background text, figure
paths, chromosome INSDC accessions, and publication DOI/date).
"""

from __future__ import annotations

import sys
import os
import re
import json
import subprocess
import configparser
import getpass
from pathlib import Path
from decimal import Decimal
import psycopg2
import yaml
from jinja2 import Environment, FileSystemLoader, StrictUndefined

# ── Tissue code → human-readable ──────────────────────────────────────────────
TISSUE_MAP = {
    "G": "gill", "L": "liver", "M": "muscle", "Sp": "spleen",
    "F": "fin", "K": "kidney", "H": "heart", "B": "blood",
}

# ── Hi-C library code → library type (used in "construction of a X library") ──
HIC_LIB_MAP = {
    "HICL": "Dovetail Hi-C LinkPrep",
    "HICL1": "Dovetail Hi-C LinkPrep",
    "HICL2": "Dovetail Hi-C LinkPrep",
    "HICL3": "Dovetail Hi-C LinkPrep",
    "HCL":  "Dovetail Hi-C",
    "OMNI": "Dovetail Omni-C",
}

# ── Hi-C library type (short name for "a X library") ─────────────────────────
HIC_LIB_TYPE_MAP = {
    "HICL":  "Dovetail Hi-C",
    "HICL1": "Dovetail Hi-C",
    "HICL2": "Dovetail Hi-C",
    "HICL3": "Dovetail Hi-C",
    "HCL":   "Dovetail Hi-C",
    "OMNI":  "Dovetail Omni-C",
}

# ── Hi-C prep module (used in "using the X and Dovetail Library Module") ──────
HIC_PREP_MODULE_MAP = {
    "HICL":  "Dovetail LinkPrep Module",
    "HICL1": "Dovetail LinkPrep Module",
    "HICL2": "Dovetail LinkPrep Module",
    "HICL3": "Dovetail LinkPrep Module",
    "HCL":   "Dovetail Hi-C Module",
    "OMNI":  "Dovetail Omni-C Module",
}

# ── Sequencing instrument → full platform name ────────────────────────────────
INSTRUMENT_MAP = {
    "iSEQ":   "Illumina iSeq 100",
    "ISEQ":   "Illumina iSeq 100",
    "NovaSeq": "Illumina NovaSeq 6000",
    "NOVA":   "Illumina NovaSeq 6000",
    "NextSeq": "Illumina NextSeq 500",
    "NEXT":   "Illumina NextSeq 500",
    "MiSeq":  "Illumina MiSeq",
    "MISEQ":  "Illumina MiSeq",
    "NovaSeqX": "Illumina NovaSeq X",
    "NOVAX":  "Illumina NovaSeq X",
}

# ── RNA tube ID tissue code → human-readable tissue name ─────────────────────
RNA_TISSUE_MAP = {
    "G": "gill",
    "M": "muscle",
    "B": "blood",
    "F": "fin clip",
    "L": "liver",
    "O": "other",
    "E": "eye",
    "H": "heart",
    "S": "skin",
    "W": "whole fish",
    "Sp": "spleen",
    "K": "kidney",
}

# ── Class → BUSCO lineage ──────────────────────────────────────────────────────
BUSCO_LINEAGE_MAP = {
    "Actinopteri":    "actinopterygii_odb10",
    "Actinopterygii": "actinopterygii_odb10",
    "Chondrichthyes": "vertebrata_odb10",
    "Mammalia":       "mammalia_odb10",
    "Reptilia":       "vertebrata_odb10",
    "Aves":           "aves_odb10",
    "Amphibia":       "tetrapoda_odb10",
}


def load_pg_config(path: str) -> dict:
    path = os.path.expanduser(path)
    cfg = configparser.ConfigParser()
    cfg.read(path)
    s = cfg["postgres"]
    return {
        "dbname": s["dbname"], "user": s["user"],
        "password": s["password"], "host": s["host"],
        "port": s.get("port", "5432"),
    }


def fmt(val, decimals: int = 2, suffix: str = "") -> str:
    """Format a numeric DB value to string; return PLACEHOLDER if None."""
    if val is None:
        return "PLACEHOLDER"
    if isinstance(val, Decimal):
        val = float(val)
    if isinstance(val, float):
        return f"{val:.{decimals}f}{suffix}"
    return f"{val}{suffix}"


def tissue_label(code: str | None) -> str:
    if code is None:
        return "PLACEHOLDER"
    return TISSUE_MAP.get(code.strip(), code)


def detect_figures(og_id: str, script_dir: Path) -> dict:
    """
    Scan script_dir/figures/ for known figure files and return a dict of
    {role: relative_path} using paths relative to script_dir (the HTML location).
    """
    figures_dir = script_dir / "figures"
    if not figures_dir.is_dir():
        return {}

    # All globs are prefixed with the OG_ID so mixed-OG figures/ dirs don't cross-contaminate
    p = og_id
    patterns = {
        "genomescope": [f"{p}*genomescope*linear*.png", f"{p}*genomescope*.png"],
        "hic":         [f"{p}*pretext_snapshotFullMap*.png", f"{p}*contact_map*.png", f"{p}*pretext*.png"],
        "merqury":     [f"{p}*spectra-asm*.png", f"{p}*spectra*.png"],
        "snail":       [f"{p}*snail*.png"],
        "blob":        [f"{p}*blob*.png"],
        "cumulative":  [f"{p}*cumulative*.png"],
        "specimen":    [f"{p}*voucher*.png", f"{p}*specimen*.png", f"{p}*photo*.png"],
    }

    found: dict[str, str] = {}
    for role, globs in patterns.items():
        for glob in globs:
            matches = sorted(figures_dir.glob(glob))
            if matches:
                found[role] = str(matches[0].relative_to(script_dir))
                break
    return found


# ── NCBI taxonomy lookup ───────────────────────────────────────────────────────

def query_ncbi_taxonomy(taxid) -> dict:
    """
    Query NCBI datasets for the full taxonomy lineage and authority.

    Two calls:
      1. Fetch the taxon to get its parents list, genus ID, and authority.
      2. Batch-fetch all parent names in one call and build the lineage string
         in order, skipping root (1), cellular organisms (131567), and the genus
         (already present in the species binomial shown separately in the template).
    """
    if not taxid or str(taxid) in ("", "PLACEHOLDER", "None"):
        return {}
    sing = os.environ.get("SING", "/software/projects/pawsey0964/singularity")
    sif = os.path.join(sing, "ncbi-datasets_18.0.2.sif")

    def _run(taxon_arg: str) -> dict:
        r = subprocess.run(
            ["singularity", "run", sif, "datasets", "summary", "taxonomy", "taxon", taxon_arg],
            capture_output=True, text=True, timeout=90,
        )
        return json.loads(r.stdout)

    try:
        # ── Call 1: get the taxon's parents list and authority ──
        data = _run(str(taxid))
        report = data["reports"][0]["taxonomy"]
        parents = report.get("parents", [])
        clf = report.get("classification", {})
        genus_id = clf.get("genus", {}).get("id")
        sci = report.get("current_scientific_name", {})
        authority = sci.get("authority", "").strip().strip("()")

        # Root (1), cellular organisms (131567), and the genus node are excluded
        skip_ids = {1, 131567}
        if genus_id:
            skip_ids.add(int(genus_id))
        lineage_ids = [p for p in parents if p not in skip_ids]

        if not lineage_ids:
            return {"full_taxonomy_string": "", "authority": authority}

        # ── Call 2: batch-resolve all lineage parent names ──
        parent_data = _run(",".join(str(p) for p in lineage_ids))
        id_to_name: dict = {}
        for pr in parent_data.get("reports", []):
            tax = pr["taxonomy"]
            name = tax.get("current_scientific_name", {}).get("name", "")
            if name:
                id_to_name[tax["tax_id"]] = name

        parts = [id_to_name[p] for p in lineage_ids if p in id_to_name]

        return {
            "full_taxonomy_string": "; ".join(parts) + ";",
            "authority": authority,
            "ncbi_order": clf.get("order", {}).get("name", ""),
        }
    except Exception as exc:
        print(f"  WARNING: NCBI taxonomy lookup failed for taxid {taxid}: {exc}", file=sys.stderr)
        return {}


# ── YAML note-input helpers ────────────────────────────────────────────────────

def _ref_label(ref: dict) -> str:
    """'Russell et al., 2010' style label from a ref dict."""
    authors = ref.get("authors", "")
    year = ref.get("year", "")
    parts = [a.strip() for a in authors.split(",") if a.strip()]
    surnames = [p.split()[0] for p in parts if p.split()]
    if len(surnames) == 0:
        label = ""
    elif len(surnames) == 1:
        label = surnames[0]
    elif len(surnames) == 2:
        label = f"{surnames[0]} &amp; {surnames[1]}"
    else:
        label = f"{surnames[0]} et al."
    return f"{label}, {year}" if (label and year) else (label or year)


def convert_background_paragraphs(paragraphs: list, ref_lookup: dict) -> str:
    """Convert YAML background paragraphs to JATS XML <p> elements.

    *text* → <italic>text</italic>
    [ref-id] → <xref ref-type="bibr" rid="ref-id">label</xref>
    """
    result = []
    for para in paragraphs:
        text = " ".join(str(para).split())
        text = re.sub(r'\*([^*]+)\*', r'<italic>\1</italic>', text)
        def _xref(m):
            rid = m.group(1)
            label = ref_lookup.get(rid, rid)
            return f'<xref ref-type="bibr" rid="{rid}">{label}</xref>'
        text = re.sub(r'\[([a-z][a-z0-9-]+)\]', _xref, text)
        result.append(f"      <p>{text}</p>")
    return "\n".join(result)


def build_additional_refs_xml(additional_refs: list) -> str:
    """Render YAML additional_refs as JATS <ref> elements for insertion into <ref-list>."""
    parts = []
    for ref in additional_refs:
        rid     = ref.get("id", "")
        authors = ref.get("authors", "")
        year    = ref.get("year", "")
        title   = ref.get("title", "")
        source  = ref.get("source", "")
        volume  = ref.get("volume", "")
        fpage   = ref.get("fpage", "")
        lpage   = ref.get("lpage", "")

        author_parts = []
        for a in authors.split(","):
            a = a.strip()
            if not a:
                continue
            tokens = a.split()
            if len(tokens) >= 2:
                sn = tokens[0]
                gn = " ".join(tokens[1:])
                author_parts.append(
                    f'<name><surname>{sn}</surname><given-names>{gn}</given-names></name>')
            else:
                author_parts.append(f'<name><surname>{a}</surname></name>')

        pub_type = "journal" if (volume or fpage) else "webpage"
        cit = [f'<element-citation publication-type="{pub_type}">']
        if author_parts:
            cit.append(f'<person-group person-group-type="author">{"".join(author_parts)}</person-group>')
        if year:
            cit.append(f'<year iso-8601-date="{year}">{year}</year>')
        if title:
            cit.append(f'<article-title>{title}</article-title>')
        if source:
            cit.append(f'<source>{source}</source>')
        if volume:
            cit.append(f'<volume>{volume}</volume>')
        if fpage:
            cit.append(f'<fpage>{fpage}</fpage>')
        if lpage:
            cit.append(f'<lpage>{lpage}</lpage>')
        cit.append('</element-citation>')

        label = _ref_label(ref)
        parts.append(
            f'        <ref id="{rid}"><label>{label}</label>{"".join(cit)}</ref>')

    return "\n".join(parts)


def load_note_input(og_id: str, script_dir: Path, metadata_path: str | None = None) -> dict:
    """Load {og_id}_note_input.yaml (or a custom --metadata path) if it exists."""
    if metadata_path:
        p = Path(metadata_path).expanduser()
    else:
        p = script_dir / f"{og_id}_note_input.yaml"
    if not p.exists():
        return {}
    with open(p) as f:
        return yaml.safe_load(f) or {}


def merge_note_input(ctx: dict, note: dict) -> dict:
    """Override ctx fields from note_input YAML, building ref lookup first."""
    if not note:
        ctx.setdefault("additional_refs_xml", "")
        _deduplicate_affiliations(ctx)
        return ctx

    additional_refs = note.get("additional_refs", [])
    ref_lookup = {r["id"]: _ref_label(r) for r in additional_refs if r.get("id")}

    if note.get("background"):
        ctx["background_paragraphs"] = convert_background_paragraphs(
            note["background"], ref_lookup)

    for role in ("collectors", "intro_writers", "curators", "contributors", "oceanomics_division"):
        if note.get(role):
            ctx[role] = note[role]

    # Scalar overrides — any non-empty YAML value wins over DB / MISSING
    _overridable = (
        # specimen / collection
        "sex_chromosomes", "collection_method", "institution", "institution_abbrev", "formal_voucher",
        "voucher_id", "latitude", "longitude",
        # RNA
        "rna_tissues",
        # figures
        "specimen_photo_path",
        # publication (filled at submission time)
        "doi", "pub_year", "pub_month", "pub_day",
    )
    for field in _overridable:
        val = note.get(field)
        # Always apply booleans (False is valid); skip None and empty strings
        if isinstance(val, bool) or (val is not None and val != ""):
            ctx[field] = val

    ctx["additional_refs_xml"] = build_additional_refs_xml(additional_refs)
    _deduplicate_affiliations(ctx)
    return ctx


def _deduplicate_affiliations(ctx: dict) -> None:
    """Assign deduplicated affiliation numbers to all authors in-place.

    Builds ctx['unique_affs'] = [{'id': 'a1', 'label': 1, 'text': '...'}]
    and sets author['aff_num'] on every entry in collectors/intro_writers/curators.
    Identical affiliation strings share the same number.
    """
    seen: dict[str, int] = {}
    unique_affs: list[dict] = []
    for role in ("collectors", "intro_writers", "curators", "contributors"):
        for author in ctx.get(role, []):
            if not isinstance(author, dict):
                continue
            aff = author.get("affiliation", "")
            if aff not in seen:
                num = len(unique_affs) + 1
                seen[aff] = num
                unique_affs.append({"id": f"a{num}", "label": num, "text": aff})
            author["aff_num"] = seen[aff]
    ctx["unique_affs"] = unique_affs


def parse_rna_tissues(og_id: str, tube_ids: list) -> str | None:
    """Derive combined tissue string from RNA library tube IDs.

    Tube ID format: {OG_ID}{TISSUE_CODE}_{rest}, e.g. OG38G_R_KL → gill.
    Returns e.g. "gill and muscle", or None if no tube IDs supplied.
    """
    if not tube_ids:
        return None
    pattern = re.compile(rf'^{re.escape(og_id)}([A-Za-z]+)', re.IGNORECASE)
    seen: list[str] = []
    for tid in sorted(tube_ids):          # sort for deterministic tissue order
        m = pattern.match(tid)
        if not m:
            continue
        code = m.group(1)
        name = RNA_TISSUE_MAP.get(code) or RNA_TISSUE_MAP.get(code.capitalize()) or code.lower()
        if name not in seen:
            seen.append(name)
    if not seen:
        return None
    if len(seen) == 1:
        return seen[0]
    if len(seen) == 2:
        return f"{seen[0]} and {seen[1]}"
    return ", ".join(seen[:-1]) + f" and {seen[-1]}"


def query_genome_note(og_id: str, pg: dict, script_dir: Path | None = None) -> dict:
    figs = detect_figures(og_id, script_dir) if script_dir else {}
    conn = psycopg2.connect(**pg)
    cur = conn.cursor()

    def q(sql: str, params=()) -> list[dict]:
        cur.execute(sql, params)
        cols = [d[0] for d in cur.description]
        return [dict(zip(cols, r)) for r in cur.fetchall()]

    def q1(sql: str, params=()) -> dict | None:
        rows = q(sql, params)
        return rows[0] if rows else None

    # ── Species & sample ──────────────────────────────────────────────────────
    smp = q1("""
        SELECT s.og_id, s.tol_id, s.common_name, s.sex,
               s.date_collected, s.collector,
               s.location, s.latitude_collection, s.longitude_collection,
               s.collection_method, s.preservation_method,
               s.voucher_id, s.tissues, s.ncbi_biosample_id,
               s.nominal_species_id
        FROM sample s WHERE s.og_id = %s
    """, (og_id,))

    sp = q1("""
        SELECT sp.species, sp.afd_common_name, sp.ncbi_taxon_id,
               sp.class, sp.ordr, sp.family, sp.genus, sp.epithet,
               sp.iucn_code
        FROM species sp
        JOIN sample s ON s.nominal_species_id = sp.species
        WHERE s.og_id = %s
    """, (og_id,))

    # ── Sequencing platform ───────────────────────────────────────────────────
    pb_seq = q1("""
        SELECT instrument, run_date
        FROM sequencing
        WHERE og_id = %s AND seq_type = 'PacBio'
        ORDER BY run_date DESC LIMIT 1
    """, (og_id,))

    # ── Raw QC (GenomeScope) ──────────────────────────────────────────────────
    rqc = q1("SELECT * FROM raw_qc WHERE og_id = %s", (og_id,))

    # ── Coverage summary ──────────────────────────────────────────────────────
    cov = q1("SELECT * FROM coverage_summary WHERE og_id = %s", (og_id,))

    # ── HiFi reads (aggregate across runs) ───────────────────────────────────
    hifi_agg = q1("""
        SELECT
            SUM(hifi_reads)           AS total_reads,
            SUM(hifi_yield) / 1e9     AS total_gb,
            MAX(tissue)               AS tissue,
            MAX(lib_code)             AS lib_code,
            STRING_AGG(run_id, ', ' ORDER BY run_id) AS run_ids
        FROM hifi_reads_qc WHERE og_id = %s
    """, (og_id,))

    # ── Hi-C reads (aggregate across lanes) ──────────────────────────────────
    hic_agg = q1("""
        SELECT
            SUM(totalreadspf)         AS total_reads,
            SUM(yield_gb)             AS total_gb,
            MAX(tissue)               AS tissue,
            MAX(lib_code)             AS lib_code,
            STRING_AGG(DISTINCT run_id, ', ' ORDER BY run_id) AS run_ids
        FROM hic_reads_qc WHERE og_id = %s
    """, (og_id,))

    # ── Hi-C sequencing platforms (shallow iSeq + deep NovaSeq/NextSeq) ──────
    hic_seq_runs = q("""
        SELECT instrument, run_date
        FROM sequencing
        WHERE og_id = %s AND seq_type = 'HiC'
        ORDER BY run_date
    """, (og_id,))

    # ── RNA tube IDs → tissue types ───────────────────────────────────────────
    rna_tube_rows = q("""
        SELECT DISTINCT rna_library_tube_id
        FROM sequencing
        WHERE og_id = %s AND rna_library_tube_id IS NOT NULL
    """, (og_id,))
    rna_tube_ids = [r["rna_library_tube_id"] for r in rna_tube_rows]

    # ── SRA accessions ────────────────────────────────────────────────────────
    sra = q("SELECT srr_accession, data_type FROM ref_genomes_sra_uploads WHERE og_id = %s", (og_id,))
    sra_by_type = {}
    for row in sra:
        sra_by_type.setdefault(row["data_type"], []).append(row["srr_accession"])

    # ── Assembly uploads (BioProject / accessions) ────────────────────────────
    uploads = q1("SELECT * FROM ref_genomes_assembly_uploads WHERE og_id = %s", (og_id,))

    # ── Assembly stats — stage 3 (curated) ───────────────────────────────────
    # sum_len = actual assembly bp (correct); num_seqs = scaffold count
    # total_scaffold_length is wrong (do not use)
    hap1 = q1("""
        SELECT * FROM ref_genomes
        WHERE og_id = %s AND stage = 3 AND haplotype = 'hap1'
    """, (og_id,))

    hap2 = q1("""
        SELECT * FROM ref_genomes
        WHERE og_id = %s AND stage = 3 AND haplotype = 'hap2'
    """, (og_id,))

    dual = q1("""
        SELECT * FROM ref_genomes
        WHERE og_id = %s AND stage = 3 AND haplotype = 'dual'
    """, (og_id,))

    busco_lineage = BUSCO_LINEAGE_MAP.get(sp["class"] if sp else "", "vertebrata_odb10")
    if hap1 and hap1.get("dataset"):
        busco_lineage = hap1["dataset"]

    # ── Mitogenome (prefer hifi tech with GenBank accession) ─────────────────
    mito = q1("""
        SELECT length, genbank_accession
        FROM mitogenome_data
        WHERE og_id = %s AND tech = 'hifi'
        ORDER BY (genbank_accession IS NOT NULL) DESC, seq_date DESC
        LIMIT 1
    """, (og_id,))

    cur.close()
    conn.close()

    # ── Build context ─────────────────────────────────────────────────────────
    species_name  = sp["species"]  if sp  else "PLACEHOLDER"
    common_name   = ((sp["afd_common_name"] if sp else None) or (smp["common_name"] if smp else None) or "PLACEHOLDER")
    ncbi_taxid    = sp["ncbi_taxon_id"] if sp else "PLACEHOLDER"
    sp_class      = sp["class"]  if sp else ""
    sp_order      = sp["ordr"]   if sp else ""
    sp_family     = sp["family"] if sp else ""

    # ── NCBI taxonomy (live lookup) ───────────────────────────────────────────
    print(f"  Fetching NCBI taxonomy for taxid {ncbi_taxid}...", file=sys.stderr)
    ncbi_tax = query_ncbi_taxonomy(ncbi_taxid)

    tolid         = smp["tol_id"]           if smp else "PLACEHOLDER"
    _sex_raw = (smp["sex"] or "") if smp else ""
    specimen_sex = {"F": "female", "M": "male"}.get(_sex_raw.strip().upper())  # None if not M/F
    collection_date = (smp["date_collected"].strftime("%Y-%m-%d")
                       if smp and smp["date_collected"] else "PLACEHOLDER")
    coll_location  = smp["location"]    if smp else "PLACEHOLDER"
    latitude       = smp["latitude_collection"]  or "PLACEHOLDER"
    longitude      = smp["longitude_collection"] or "PLACEHOLDER"
    # Strip voyage/trip info after the person's name (e.g. "Glenn Moore, SW Voyage 2023" → "Glenn Moore")
    collector_raw  = smp["collector"] if smp else "PLACEHOLDER"
    collector_name = collector_raw.split(",")[0].strip() if collector_raw else "PLACEHOLDER"
    voucher_id     = smp["voucher_id"]  if smp else "PLACEHOLDER"
    pres_method    = smp["preservation_method"] or "PLACEHOLDER"
    hifi_tissue_code = hifi_agg["tissue"] if hifi_agg else None
    hic_tissue_code  = hic_agg["tissue"]  if hic_agg else None

    # PacBio platform
    pacbio_platform = ("PacBio " + pb_seq["instrument"].strip()
                       if pb_seq and pb_seq["instrument"] else "PLACEHOLDER")

    # Hi-C library name / type / prep module from lib_code
    hic_lib_code = hic_agg["lib_code"] if hic_agg else None
    hic_lib_name = HIC_LIB_MAP.get(hic_lib_code or "", "Hi-C")
    hic_lib_type = HIC_LIB_TYPE_MAP.get(hic_lib_code or "", "Dovetail Hi-C")
    hic_prep_module = HIC_PREP_MODULE_MAP.get(hic_lib_code or "", "Dovetail LinkPrep Module")

    # Separate shallow (iSeq) and deep (NovaSeq/NextSeq) runs from sequencing table
    hic_shallow_platform = None
    hic_deep_platform    = None
    for run in hic_seq_runs:
        instr = run["instrument"] or ""
        mapped = INSTRUMENT_MAP.get(instr, instr)
        if "iSeq" in mapped or "iSEQ" in instr.upper():
            hic_shallow_platform = mapped
        else:
            hic_deep_platform = mapped
    # Fallback: derive deep platform from hic_reads_qc run_id if sequencing table incomplete
    if not hic_deep_platform and hic_agg and hic_agg["run_ids"]:
        prefix = hic_agg["run_ids"].split("_")[0].upper()
        hic_deep_platform = INSTRUMENT_MAP.get(prefix, f"Illumina ({prefix})")

    # GenomeScope
    gs_size  = rqc["genomesize"] if rqc else None
    gs_het   = rqc["heterozygosity"] if rqc else None
    gs_rep   = rqc["repeatsize"]  if rqc else None

    # Coverage
    hifi_gb_val  = cov["total_hifi_yield_gb"] if cov else None
    hic_gb_val   = cov["total_hic_yield_gb"]  if cov else None
    hifi_cov_val = cov["hifi_coverage"]        if cov else None
    hic_cov_val  = cov["hic_coverage"]         if cov else None

    # HiFi reads
    hifi_reads_val = hifi_agg["total_reads"] if hifi_agg else None
    # HiC reads (sum of totalreadspf across lanes)
    hic_reads_val  = hic_agg["total_reads"]  if hic_agg else None

    # SRR
    hifi_srr = ", ".join(sra_by_type.get("hifi", ["PLACEHOLDER"]))
    hic_srr  = ", ".join(sra_by_type.get("hic",  ["PLACEHOLDER"]))

    # BioSamples / BioProjects
    biosample      = uploads["biosample"]               if uploads else "PLACEHOLDER"
    hap1_bioproject = uploads["bioproject_hap1"]        if uploads else "PLACEHOLDER"
    hap2_bioproject = uploads["bioproject_hap2"]        if uploads else "PLACEHOLDER"
    hap1_accession  = uploads["assembly_accession_hap1"] if uploads else "PLACEHOLDER"
    hap2_accession  = uploads["assembly_accession_hap2"] if uploads else "PLACEHOLDER"

    # Assembly names (ToLID-based convention)
    hap1_name = f"{tolid}.hap1" if tolid != "PLACEHOLDER" else "PLACEHOLDER"
    hap2_name = f"{tolid}.hap2" if tolid != "PLACEHOLDER" else "PLACEHOLDER"

    # Mito
    mito_len_kb = "PLACEHOLDER"
    if mito:
        mito_len_kb = f"{mito['length'] / 1e3:.2f}"
    mito_genbank = mito["genbank_accession"] if mito and mito["genbank_accession"] else "PLACEHOLDER"

    # Hap1 stats
    def hap_stats(h: dict | None) -> dict:
        if not h:
            return {k: "PLACEHOLDER" for k in [
                "total_length_mb", "scaffold_n50_mb", "contig_n50_mb",
                "num_scaffolds", "num_gaps", "qv", "completeness", "busco_single",
                "busco_dup", "busco_frag", "busco_missing", "pct_assigned",
                "num_chromosomes", "kmer_completeness",
            ]}
        total_mb = f"{h['sum_len'] / 1e6:.2f}" if h.get("sum_len") else "PLACEHOLDER"
        return {
            "total_length_mb":   total_mb,
            "scaffold_n50_mb":   fmt(h.get("scaffold_n50_size_mb")),
            "contig_n50_mb":     fmt(h.get("contig_n50_size_mb")),
            "num_scaffolds":     str(h["num_seqs"]) if h.get("num_seqs") else "PLACEHOLDER",
            "num_gaps":          str(h["num_gaps"]) if h.get("num_gaps") is not None else "PLACEHOLDER",
            "qv":                fmt(h.get("qv")),
            "completeness":      fmt(h.get("complete") or h.get("completeness")),
            "busco_single":      fmt(h.get("single_copy")),
            "busco_dup":         fmt(h.get("multi_copy")),
            "busco_frag":        fmt(h.get("fragmented")),
            "busco_missing":     fmt(h.get("missing")),
            "pct_assigned":      fmt(h.get("pct_assigned")),
            "num_chromosomes":   str(h["num_chromosomes"]) if h.get("num_chromosomes") else "PLACEHOLDER",
            "kmer_completeness": fmt(h.get("completeness")),
        }

    h1 = hap_stats(hap1)
    h2 = hap_stats(hap2)

    # Combined stats from dual row
    combined_qv            = fmt(dual.get("qv"))           if dual else "PLACEHOLDER"
    combined_completeness  = fmt(dual.get("completeness"))  if dual else "PLACEHOLDER"

    # BUSCO complete % (single + dup)
    def busco_complete(h_stats: dict) -> str:
        try:
            return f"{float(h_stats['busco_single']) + float(h_stats['busco_dup']):.1f}"
        except (ValueError, TypeError):
            return "PLACEHOLDER"

    # Use NCBI order if available (DB ordr field is sometimes an incertae sedis group)
    effective_order = ncbi_tax.get("ncbi_order") or sp_order

    # Higher taxonomy string (brief, for body text and keywords)
    higher_tax_parts = [p for p in [sp_class, effective_order, sp_family] if p]
    higher_taxonomy = "; ".join(higher_tax_parts) if higher_tax_parts else "PLACEHOLDER"

    def missing(description: str, source: str) -> str:
        """Return a clearly labelled placeholder."""
        return f"[MISSING: {description} — {source}]"

    ctx = {
        # ── Publication (filled at submission) ──────────────────────────────
        "doi":           missing("DOI", "assigned by journal at publication"),
        "pub_day":       missing("publication day", "assigned at publication"),
        "pub_month":     missing("publication month", "assigned at publication"),
        "pub_year":      missing("publication year", "assigned at publication"),

        # ── Journal / centre ─────────────────────────────────────────────────
        "centre_email":  "info@oceanomics.au",

        # ── Species ──────────────────────────────────────────────────────────
        "species_name":       species_name,
        "common_name":        common_name,
        "authority":          ncbi_tax.get("authority") or missing(
                                  "taxonomic authority e.g. (Valenciennes, 1839)",
                                  "check FishBase or NCBI taxonomy page for this species"),
        "ncbi_taxid":         ncbi_taxid,
        "order":              effective_order,
        "higher_taxonomy":    higher_taxonomy,
        "full_taxonomy_string": ncbi_tax.get("full_taxonomy_string") or missing(
                                    "full NCBI taxonomy lineage",
                                    f"copy from https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id={ncbi_taxid}"),

        # ── Specimen ─────────────────────────────────────────────────────────
        "tolid":              tolid,
        "specimen_sex":       specimen_sex,  # None when not M/F — template omits it
        "collection_date":    collection_date,
        "collection_location": coll_location or missing("collection location name", "sample.location in DB"),
        "latitude":           latitude if latitude != "PLACEHOLDER" else
                              missing("collection latitude", "sample.latitude_collection is NULL — check field records"),
        "longitude":          longitude if longitude != "PLACEHOLDER" else
                              missing("collection longitude", "sample.longitude_collection is NULL — check field records"),
        "collector_names":    collector_name,
        "institution":        missing("full institution name e.g. Western Australia Museum",
                                      "not stored in DB — check sample.collector or field records"),
        "institution_abbrev": missing("institution abbreviation e.g. WAM",
                                      "not stored in DB — check sample.collector or field records"),
        "formal_voucher":     True,
        "collection_method":  smp["collection_method"] if smp else missing("collection method", "sample.collection_method in DB"),
        "preservation_method": pres_method if pres_method != "PLACEHOLDER" else
                               missing("preservation method", "sample.preservation_method in DB"),
        "voucher_id":         voucher_id if voucher_id else missing("museum voucher ID", "sample.voucher_id in DB"),
"hifi_tissue":        tissue_label(hifi_tissue_code),
        "hic_tissue":         tissue_label(hic_tissue_code),
        "pacbio_platform":    pacbio_platform,
        "hic_lib_name":       hic_lib_name,
        "hic_lib_type":       hic_lib_type,
        "hic_prep_module":    hic_prep_module,
        "hic_shallow_platform": hic_shallow_platform,   # None = no shallow QC run
        "hic_deep_platform":  hic_deep_platform or missing("Hi-C deep sequencing platform", "sequencing table — no HiC runs found"),
        "hic_platform":       f"{hic_lib_name} / {hic_deep_platform}" if hic_deep_platform else "MISSING",

        # ── RNA sequencing (None = no RNA data → section omitted) ───────────
        "rna_tissues":        parse_rna_tissues(og_id, rna_tube_ids),

        # ── Sequencing stats ─────────────────────────────────────────────────
        "hifi_gb":            fmt(hifi_gb_val),
        "hifi_read_count":    f"{int(hifi_reads_val):,}" if hifi_reads_val else "PLACEHOLDER",
        "hic_gb":             fmt(hic_gb_val),
        "hic_read_count":     f"{int(hic_reads_val):,}" if hic_reads_val else "PLACEHOLDER",
        "hifi_biosample":     biosample,
        "hic_biosample":      biosample,
        "hifi_run_accession": hifi_srr,
        "hic_run_accession":  hic_srr,

        # ── GenomeScope ──────────────────────────────────────────────────────
        "genomescope_size_bp":  f"{gs_size:,}" if gs_size else "PLACEHOLDER",
        "genomescope_size_mb":  f"{gs_size / 1e6:.0f}" if gs_size else "PLACEHOLDER",
        "heterozygosity_pct":   fmt(gs_het),
        "repeat_size_bp":       f"{gs_rep:,}" if gs_rep else "PLACEHOLDER",
        "hifi_coverage":        fmt(hifi_cov_val, decimals=0),
        "hic_coverage":         fmt(hic_cov_val,  decimals=0),

        # ── Assembly accessions ──────────────────────────────────────────────
        "hap1_bioproject":      hap1_bioproject or "PLACEHOLDER",
        "hap2_bioproject":      hap2_bioproject or "PLACEHOLDER",
        "hap1_assembly_name":   hap1_name,
        "hap2_assembly_name":   hap2_name,

        # ── Assembly stats (hap1) ────────────────────────────────────────────
        "hap1_total_length_mb":   h1["total_length_mb"],
        "hap1_scaffold_n50_mb":   h1["scaffold_n50_mb"],
        "hap1_contig_n50_mb":     h1["contig_n50_mb"],
        "hap1_num_scaffolds":     h1["num_scaffolds"],
        "hap1_num_gaps":          h1["num_gaps"],
        "hap1_qv":                h1["qv"],
        "hap1_busco_complete_pct": busco_complete(h1),
        "hap1_busco_single":      h1["busco_single"],
        "hap1_busco_dup":         h1["busco_dup"],
        "hap1_busco_frag":        h1["busco_frag"],
        "hap1_busco_missing":     h1["busco_missing"],
        "hap1_pct_assigned":      h1["pct_assigned"],
        "hap1_num_autosomes":     h1["num_chromosomes"],
        "hap1_kmer_completeness": h1["kmer_completeness"],

        # ── Assembly stats (hap2) ────────────────────────────────────────────
        "hap2_total_length_mb":   h2["total_length_mb"],
        "hap2_scaffold_n50_mb":   h2["scaffold_n50_mb"],
        "hap2_contig_n50_mb":     h2["contig_n50_mb"],
        "hap2_num_scaffolds":     h2["num_scaffolds"],
        "hap2_num_gaps":          h2["num_gaps"],
        "hap2_qv":                h2["qv"],
        "hap2_busco_complete_pct": busco_complete(h2),
        "hap2_busco_single":      h2["busco_single"],
        "hap2_busco_dup":         h2["busco_dup"],
        "hap2_pct_assigned":      h2["pct_assigned"],
        "hap2_kmer_completeness": h2["kmer_completeness"],

        # ── Combined ─────────────────────────────────────────────────────────
        "combined_qv":              combined_qv,
        "combined_kmer_completeness": combined_completeness,
        # diploid count = haploid × 2
        "chromosome_count":         str(int(h1["num_chromosomes"]) * 2) if h1["num_chromosomes"] not in ("PLACEHOLDER", None) else "PLACEHOLDER",
        "sex_chromosomes":          missing("sex chromosomes identified e.g. 'Not determined' or 'XY'",
                                             "from manual curation — not in DB"),

        # ── BUSCO lineage ────────────────────────────────────────────────────
        "busco_lineage":            busco_lineage,

        # ── Mitogenome ───────────────────────────────────────────────────────
        "mito_length_kb":           mito_len_kb,
        "mito_genbank_accession":   mito_genbank,

        # ── Authors (ORCIDs not in DB — add manually) ────────────────────────
        "collectors": [
            {
                "surname":    collector_name.split()[-1] if collector_name not in ("PLACEHOLDER", "") else missing("collector surname", "sample.collector in DB"),
                "given_names": " ".join(collector_name.split()[:-1]) if collector_name not in ("PLACEHOLDER", "") else missing("collector given names", "sample.collector in DB"),
                "orcid":      missing(f"ORCID for {collector_name}", "check with collector"),
                "affiliation": missing(f"institution for {collector_name} e.g. Western Australia Museum, Perth, WA, Australia",
                                       "not in DB — check with collector"),
            }
        ],
        "intro_writers": [
            {
                "surname":    missing("background author surname", "person who wrote the species background"),
                "given_names": missing("background author given names", "person who wrote the species background"),
                "orcid":      missing("background author ORCID", "check with author"),
                "affiliation": missing("background author institution", "check with author"),
            }
        ],
        "curators": [
            {
                "surname": missing("curator surname", "one curator — add to YAML"),
                "given_names": missing("curator given names", "one curator — add to YAML"),
                "orcid":   "",
                "affiliation": "OceanOmics Centre, University of Western Australia, Perth, WA 6009, Australia",
            }
        ],
        "contributors": [],
        "oceanomics_division": [],

        # ── Figures (auto-detected from figures/ dir, else MISSING) ─────────────
        "specimen_photo_path":        figs.get("specimen",
                                               missing(f"specimen photo for {tolid}",
                                                       "from WAM or collector — provide image path/URL")),
        "genomescope_plot_path":      figs.get("genomescope",
                                               missing(f"{og_id} GenomeScope linear plot PNG",
                                                       "from genome notes pipeline genomescope output")),
        "hap1_hic_contact_map_path":  figs.get("hic",
                                               missing(f"{tolid}.hap1 Hi-C contact map PNG",
                                                       "from post-curation pipeline PretextView snapshot output")),
        "merqury_kmer_plot_path":     figs.get("merqury",
                                               missing(f"{tolid} MerquryFK spectra-cn plot PNG",
                                                       "from post-curation pipeline MerquryFK output")),
        "hap1_snail_plot_path":       figs.get("snail",
                                               missing(f"{tolid}.hap1 BlobToolKit snail plot PNG",
                                                       "from genome notes pipeline blobtools output")),
        "hap1_blob_plot_path":        figs.get("blob",
                                               missing(f"{tolid}.hap1 BlobToolKit GC-coverage blob plot PNG",
                                                       "from genome notes pipeline blobtools output")),

        # ── Background (overridden by YAML if provided) ───────────────────────
        "background_paragraphs": f"""<!-- MANUAL ENTRY REQUIRED: Write 3-5 paragraphs on {species_name}.
     Include: distribution, ecology, morphology, behaviour, diet, reproduction,
     fisheries/conservation status (IUCN: {sp['iucn_code'] if sp else 'check species table'}).
     Use <p> tags. Use <italic> for species names. Add <xref ref-type="bibr"> for citations. -->
<p>[MISSING: species background text for {species_name}]</p>""",

        # ── Additional refs (populated from YAML additional_refs) ─────────────
        "additional_refs_xml": "",
    }

    return ctx


def _bold(v: str) -> str:
    """Wrap a DB-pulled value in a custom element for review mode."""
    return f"<db-val>{v}</db-val>"


_ATTRIBUTE_FIELDS = {
    # These appear inside XML attribute values in the template — wrapping breaks XML
    "centre_email", "specimen_photo_path", "genomescope_plot_path",
    "hap1_hic_contact_map_path", "merqury_kmer_plot_path",
    "hap1_snail_plot_path", "hap1_blob_plot_path",
    "hap1_bioproject", "hap2_bioproject",
}


def apply_review_highlighting(ctx: dict) -> dict:
    """Return a copy of ctx with all non-MISSING scalar DB values wrapped in <db-val>."""
    out = {}
    for k, v in ctx.items():
        if (k not in _ATTRIBUTE_FIELDS
                and isinstance(v, str)
                and not v.startswith("[MISSING:")
                and v != "PLACEHOLDER"
                and "<" not in v
                and v.strip()):
            out[k] = _bold(v)
        else:
            out[k] = v
    return out


def main():
    flags = sys.argv[1:]
    review_mode = "--review" in flags

    # --metadata /path/to/file.yaml — extract value and remove both from positional args
    metadata_path: str | None = None
    skip_next = False
    for i, f in enumerate(flags):
        if skip_next:
            skip_next = False
            continue
        if f == "--metadata" and i + 1 < len(flags):
            metadata_path = flags[i + 1]
            skip_next = True

    exclude = set()
    for i, f in enumerate(flags):
        if f.startswith("--"):
            exclude.add(i)
            if f == "--metadata" and i + 1 < len(flags):
                exclude.add(i + 1)
    args = [f for i, f in enumerate(flags) if i not in exclude]

    if len(args) < 1:
        print(
            "Usage: uv run --script generate_genome_note.py <OG_ID> [pg_config] "
            "[--review] [--metadata /path/to/OG38_note_input.yaml]",
            file=sys.stderr)
        sys.exit(1)

    og_id   = args[0].strip()
    pg_path = args[1] if len(args) > 1 else "~/postgresql_details/oceanomics.cfg"
    pg = load_pg_config(pg_path)

    script_dir = Path(__file__).parent
    print(f"Querying DB for {og_id}...")
    ctx = query_genome_note(og_id, pg, script_dir=script_dir)

    note = load_note_input(og_id, script_dir, metadata_path)
    if note:
        print(f"  Loaded note input: {metadata_path or script_dir / f'{og_id}_note_input.yaml'}")
    ctx = merge_note_input(ctx, note)

    template_path = script_dir / "oceangenomes_genome_note.xml.j2"
    if not template_path.exists():
        print(f"Template not found: {template_path}", file=sys.stderr)
        sys.exit(1)

    env = Environment(
        loader=FileSystemLoader(str(script_dir)),
        undefined=StrictUndefined,
        trim_blocks=True,
        lstrip_blocks=True,
        autoescape=False,
    )
    template = env.get_template("oceangenomes_genome_note.xml.j2")

    render_ctx = apply_review_highlighting(ctx) if review_mode else ctx
    xml_out = template.render(**render_ctx)

    suffix = "_review" if review_mode else ""
    out_file = script_dir / f"{og_id}_genome_note{suffix}.xml"
    out_file.write_text(xml_out, encoding="utf-8")
    print(f"Written: {out_file}")

    # Report what still needs manual input
    missing_fields = [(k, v) for k, v in ctx.items()
                      if isinstance(v, str) and v.startswith("[MISSING:")]
    if missing_fields:
        print(f"\nFields requiring manual input ({len(missing_fields)}):")
        for k, v in missing_fields:
            print(f"  {k}:\n    {v}")


if __name__ == "__main__":
    main()
