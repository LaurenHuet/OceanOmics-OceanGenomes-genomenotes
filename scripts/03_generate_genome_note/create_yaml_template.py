#!/usr/bin/env python3
# /// script
# requires-python = ">=3.10"
# dependencies = ["psycopg2-binary", "pyyaml"]
# ///
"""
create_yaml_template.py — generate {OG_ID}_note_input.yaml from the database.

All values come from the DB. Fields that are genuinely not in the DB
(background text, sex_chromosomes from curation, ORCIDs) are left as
REPLACE placeholders for manual entry.

Usage:
    uv run --script create_yaml_template.py <OG_ID> [pg_config]

pg_config defaults to ~/postgresql_details/oceanomics.cfg
"""

import sys
import configparser
import textwrap
from pathlib import Path

import psycopg2
import yaml


def load_pg(path: str) -> dict:
    cfg = configparser.ConfigParser()
    cfg.read(Path(path).expanduser())
    return dict(cfg["postgres"])


def query(cur, sql, params=()):
    cur.execute(sql, params)
    cols = [d[0] for d in cur.description]
    return [dict(zip(cols, r)) for r in cur.fetchall()]


def q1(cur, sql, params=()):
    rows = query(cur, sql, params)
    return rows[0] if rows else None


def parse_collector(raw: str):
    """Parse the DB collector field into (given, surname, affiliation).

    The field often contains: "Given Surname, email@example.com, Institution Name"
    We take the first comma-delimited token as the name, strip any email token,
    and use the remaining parts as the affiliation.
    """
    parts = [p.strip() for p in (raw or "").split(",")]
    if not parts or not parts[0]:
        return "REPLACE", "REPLACE", "REPLACE"

    # First token: the person's name
    name_parts = parts[0].split()
    if len(name_parts) >= 2:
        given   = " ".join(name_parts[:-1])
        surname = name_parts[-1]
    else:
        given   = parts[0]
        surname = "REPLACE"

    # Remaining tokens: skip email addresses, join the rest as affiliation
    import re
    aff_parts = [p for p in parts[1:] if p and not re.match(r"[^@]+@[^@]+", p)]
    affiliation = ", ".join(aff_parts) if aff_parts else "REPLACE — check with collector"

    return given, surname, affiliation


def main():
    og_id = sys.argv[1] if len(sys.argv) > 1 else None
    if not og_id:
        sys.exit("Usage: uv run --script create_yaml_template.py <OG_ID> [pg_config]")

    pg_path = sys.argv[2] if len(sys.argv) > 2 else "~/postgresql_details/oceanomics.cfg"
    pg = load_pg(pg_path)
    conn = psycopg2.connect(**pg)
    cur = conn.cursor()

    # ── Sample ────────────────────────────────────────────────────────────────
    smp = q1(cur, """
        SELECT s.og_id, s.tol_id, s.common_name, s.sex,
               s.collector, s.location,
               s.latitude_collection, s.longitude_collection,
               s.collection_method, s.preservation_method, s.voucher_id
        FROM sample s WHERE s.og_id = %s
    """, (og_id,))

    if not smp:
        sys.exit(f"OG_ID {og_id!r} not found in the database.")

    # ── Species ───────────────────────────────────────────────────────────────
    sp = q1(cur, """
        SELECT sp.species, sp.afd_common_name
        FROM species sp
        JOIN sample s ON s.nominal_species_id = sp.species
        WHERE s.og_id = %s
    """, (og_id,))

    species_name   = sp["species"] if sp else "REPLACE"
    common_name    = sp["afd_common_name"] or smp.get("common_name") or "REPLACE"

    # ── Collector ─────────────────────────────────────────────────────────────
    collector_raw  = smp.get("collector") or ""
    coll_given, coll_surname, coll_affiliation = parse_collector(collector_raw)

    conn.close()

    # ── Output path ───────────────────────────────────────────────────────────
    out_path = Path(__file__).parent / f"{og_id}_note_input.yaml"

    # ── Build YAML content as a string (preserves comments and ordering) ──────
    content = textwrap.dedent(f"""\
        # ── Manual input for {og_id} genome note ─────────────────────────────────
        # {species_name}, {common_name}
        # All values below from DB except where marked REPLACE.
        # REPLACE fields need manual entry before running the generator.

        # Write 2–4 paragraphs on the species: distribution, ecology,
        # morphology, behaviour, diet, reproduction, fisheries/conservation.
        # Use *italics* for species names. Use [ref-id] for citations.
        # Add any cited references to additional_refs at the bottom.
        background:
          - "REPLACE"

        # ── Authors ──────────────────────────────────────────────────────────
        # Collector from DB: {collector_raw or 'not in DB'}
        collectors:
          - surname: {coll_surname}
            given_names: {coll_given}
            affiliation: "{coll_affiliation}"
            orcid: ""

        intro_writers:
          - surname: REPLACE
            given_names: REPLACE
            affiliation: "OceanOmics Centre, University of Western Australia, Perth, WA 6009, Australia"
            orcid: ""

        # One curator — the person who performed the genome curation
        curators:
          - surname: REPLACE
            given_names: REPLACE
            affiliation: "OceanOmics Centre, University of Western Australia, Perth, WA 6009, Australia"
            orcid: ""

        # OceanOmics Centre team who contributed to this genome
        contributors:
          - surname: Anstiss
            given_names: Liam
            affiliation: "OceanOmics Centre, University of Western Australia, Perth, WA 6009, Australia"
            orcid: ""
          - surname: de Jong
            given_names: Emma
            affiliation: "OceanOmics Centre, University of Western Australia, Perth, WA 6009, Australia"
            orcid: ""
          - surname: Depiazzi
            given_names: Anna
            affiliation: "OceanOmics Centre, University of Western Australia, Perth, WA 6009, Australia"
            orcid: ""
          - surname: Doran
            given_names: Adrianne
            affiliation: "OceanOmics Centre, University of Western Australia, Perth, WA 6009, Australia"
            orcid: ""
          - surname: Faseeh
            given_names: Ibrahim
            affiliation: "OceanOmics Centre, University of Western Australia, Perth, WA 6009, Australia"
            orcid: ""
          - surname: Missen
            given_names: Laura
            affiliation: "OceanOmics Centre, University of Western Australia, Perth, WA 6009, Australia"
            orcid: ""
          - surname: Peirce
            given_names: Tyler
            affiliation: "OceanOmics Centre, University of Western Australia, Perth, WA 6009, Australia"
            orcid: ""
          - surname: Nester
            given_names: Georgia M.
            affiliation: "OceanOmics Centre, University of Western Australia, Perth, WA 6009, Australia"
            orcid: ""
          - surname: Thorpe
            given_names: Ebony M.
            affiliation: "OceanOmics Centre, University of Western Australia, Perth, WA 6009, Australia"
            orcid: ""
          - surname: Bennett
            given_names: Adam J.
            affiliation: "OceanOmics Centre, University of Western Australia, Perth, WA 6009, Australia"
            orcid: ""
          - surname: Nguyen
            given_names: Olivia
            affiliation: "OceanOmics Centre, University of Western Australia, Perth, WA 6009, Australia"
            orcid: ""
          - surname: Wong
            given_names: Tsz Ching
            affiliation: "OceanOmics Centre, University of Western Australia, Perth, WA 6009, Australia"
            orcid: ""
          - surname: Corrigan
            given_names: Shannon
            affiliation: "OceanOmics Centre, University of Western Australia, Perth, WA 6009, Australia"
            orcid: ""

        # OceanOmics Division — always included
        oceanomics_division:
          - surname: Ayad
            given_names: Marcelle E.
            affiliation: "OceanOmics Centre, University of Western Australia, Perth, WA 6009, Australia"
            orcid: ""
          - surname: Bayer
            given_names: Philipp E.
            affiliation: "OceanOmics Centre, University of Western Australia, Perth, WA 6009, Australia"
            orcid: ""
          - surname: Bunce
            given_names: Michael
            affiliation: "OceanOmics Centre, University of Western Australia, Perth, WA 6009, Australia"
            orcid: ""
          - surname: Fraser
            given_names: Matthew W.
            affiliation: "OceanOmics Centre, University of Western Australia, Perth, WA 6009, Australia"
            orcid: ""
          - surname: Goncalves
            given_names: Priscila
            affiliation: "OceanOmics Centre, University of Western Australia, Perth, WA 6009, Australia"
            orcid: ""
          - surname: Raes
            given_names: Eric J.
            affiliation: "OceanOmics Centre, University of Western Australia, Perth, WA 6009, Australia"
            orcid: ""

        # ── Collection details ────────────────────────────────────────────────
        # sex_chromosomes: result from manual curation — not in DB
        sex_chromosomes: "Not determined"

        # formal_voucher: true = museum accession (e.g. WAM), false = "specimen held at institution"
        formal_voucher: true
        institution: "REPLACE"
        institution_abbrev: "REPLACE"

        # ── Figures ──────────────────────────────────────────────────────────
        # Leave blank — auto-detected from figures/ directory.
        # Drop PNGs named {og_id}*.png into genome_notes_automation/figures/
        specimen_photo_path: ""

        # ── Publication details (fill at submission) ──────────────────────────
        doi: ""
        pub_year: ""
        pub_month: ""
        pub_day: ""

        # ── Species-specific references ───────────────────────────────────────
        # Add any references cited in the background paragraphs above.
        # Format:
        #   - id: ref-smith2010
        #     authors: "Smith A, Jones B"
        #     year: "2010"
        #     title: "Title of paper"
        #     source: "Journal Name"
        #     volume: "10"
        #     fpage: "100"
        #     lpage: "110"
        additional_refs: []
    """)

    out_path.write_text(content)
    print(f"Written: {out_path}")
    print()
    print(f"  Species:    {species_name}")
    print(f"  Common name:{common_name}")
    print(f"  Collector:  {collector_raw or '(not in DB)'}")
    print()
    print("Fields to fill before running the generator:")
    print("  background          — species biology text")
    print("  collectors.affiliation — institution of collector")
    print("  intro_writers       — who wrote the background")
    print("  curators            — genome curation team")
    print("  sex_chromosomes     — from curation output")
    print("  institution / institution_abbrev — voucher holding institution")


if __name__ == "__main__":
    main()
