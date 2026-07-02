# CLAUDE.md — OceanOmics Genome Note Generation

## What this repo is

`OceanOmics-OceanGenomes-genomenotes` is a Nextflow DSL2 pipeline that runs the
**computational steps** for genome notes (blobtools, BUSCO, alignment, window stats, etc.).

The scripts in `scripts/03_generate_genome_note/` handle the **document generation**
side: querying the database, filling a JATS XML template, and rendering to HTML and Word.
These two parts need to be wired together into a single automated pipeline — that is the
main outstanding work described below.

---

## Current manual pipeline (working as of 2026-07)

### Step 1 — Run the Nextflow genomenotes pipeline

Generates blobtools plots, BUSCO, window stats. Output lands in:
```
/scratch/pawsey0964/lhuet/genomenotes/OG{ID}/
```

### Step 2 — Stage figures for the note

Edit `scripts/03_generate_genome_note/stage_figures.conf` (set `OG_ID` and `OUTPUT_DIR`),
then run:
```bash
bash scripts/03_generate_genome_note/stage_figures.sh
```
Downloads GenomeScope plot, Hi-C pretext map, and Merqury spectra plot from Acacia
(`pawsey0964:oceanomics-refassemblies/{OG_ID}/`). Copy blob/snail/cumulative PNGs from
the Nextflow output into the same `figures/` directory.

### Step 3 — Generate YAML metadata template

```bash
uv run --script scripts/03_generate_genome_note/create_yaml_template.py <OG_ID> \
    ~/postgresql_details/oceanomics.cfg
```
Queries the OceanOmics PostgreSQL DB and writes `{OG_ID}_note_input.yaml` pre-filled
with species, common name, and collector from the DB. Manual fields to complete:
- `background` — 2–4 paragraphs on species biology (use `*italics*` for species names)
- `intro_writers` — who wrote the background
- `curators` — who performed genome curation
- `sex_chromosomes` — from curation output
- `formal_voucher`, `institution`, `institution_abbrev` — voucher holding details
- `doi`, `pub_year`, `pub_month`, `pub_day` — fill at submission

### Step 4 — Generate the JATS XML

```bash
uv run --script scripts/03_generate_genome_note/generate_genome_note.py <OG_ID> \
    ~/postgresql_details/oceanomics.cfg \
    --metadata {OG_ID}_note_input.yaml
```
Pulls all assembly stats, taxonomy, sequencing metadata from DB and renders
`oceangenomes_genome_note.xml.j2` → `{OG_ID}_genome_note.xml`.

Figures are auto-detected from a `figures/` directory (prefix-matched on OG_ID to
prevent cross-contamination). Drop PNGs named `{OG_ID}*.png` into `figures/` before
running.

### Step 5 — Render to HTML and Word

```bash
# HTML preview
uv run --script scripts/03_generate_genome_note/render_genome_note.py \
    {OG_ID}_genome_note.xml

# Word document (requires python-docx + pillow; pillow needed for HEIC photos)
uv run --with python-docx --with pillow \
    scripts/03_generate_genome_note/render_genome_note.py \
    {OG_ID}_genome_note.xml --docx
```

**Note:** Voucher photos from iPhone are HEIC format with a .png extension.
Convert before use:
```bash
uv run --with pillow --with pillow-heif python3 - <<'EOF'
import pillow_heif; pillow_heif.register_heif_opener()
from PIL import Image
p = "figures/{OG_ID}_voucher_photo.png"
Image.open(p).convert("RGB").save(p)
EOF
```

### Step 6 — Acacia backup

```bash
DEST="pawsey0964:oceanomics-genomenotes/{OG_ID}"
GN="/scratch/pawsey0964/lhuet/genomenotes"

# Pipeline results (intermediate + final, no raw reads)
rclone copy $GN/blobtoolkit/{OG_ID}_window_stats.tsv  $DEST/blobtoolkit/
rclone copy $GN/{OG_ID}/busco/                        $DEST/busco/
rclone copy $GN/{OG_ID}/blobtools/                    $DEST/blobtools/

# Blast results — zip the per-part TXT files then upload
find $GN -path "*/{OG_ID}_${OG_ID}.part_*/blobtools/blast/results/*.txt" | sort \
    | zip -j /tmp/{OG_ID}_blast_results.zip -@
rclone copy /tmp/{OG_ID}_blast_results.zip $DEST/blobtools/blast/

# Genome note documents
rclone copy {OG_ID}_genome_note.xml   $DEST/genome_note/
rclone copy {OG_ID}_genome_note.docx  $DEST/genome_note/
rclone copy {OG_ID}_genome_note.html  $DEST/genome_note/
rclone copy {OG_ID}_note_input.yaml   $DEST/genome_note/
rclone copy figures/                  $DEST/genome_note/figures/ --include "{OG_ID}*"
```

---

## Author structure (all genome notes)

Authors go in this order in the YAML:
1. `collectors` — specimen collector (from DB)
2. `intro_writers` — who wrote the Background section
3. `curators` — the genome curator (1 person)
4. `contributors` — OceanOmics Centre team (12 people; **Shannon Corrigan always last** — she is centre lead)
5. `oceanomics_division` — OceanOmics Division (6 people, listed in Author Information section only)

Default contributor and division lists are pre-filled by `create_yaml_template.py`.

---

## Future work — Nextflow integration

The goal is to wire `scripts/03_generate_genome_note/` into the Nextflow pipeline so
the genome note document is generated automatically at the end of a pipeline run.

### What needs to happen

1. **Add a `GENERATE_GENOME_NOTE` process** that runs `generate_genome_note.py`
   and `render_genome_note.py` using a container with psycopg2 + jinja2 + python-docx
   + pillow pre-installed.

2. **Figure paths**: The pipeline already produces blob/snail/cumulative PNGs and the
   GenomeScope plot. Pass these as process inputs; the generator currently auto-detects
   them from a `figures/` directory by glob pattern `{OG_ID}*.png`.

3. **YAML metadata**: `create_yaml_template.py` handles all DB-sourced fields
   automatically. The manual fields (background text, sex_chromosomes) still require
   human input — the pipeline can emit a draft with `REPLACE` placeholders, then
   re-render after manual completion.

4. **DB credentials**: Currently uses psycopg2 singularity container:
   ```
   singularity run $SING/psycopg2:0.1.sif python generate_genome_note.py ...
   ```
   Needs to be adapted to the Nextflow `container` directive and secret handling.

5. **Output**: Final outputs:
   - `{OG_ID}_genome_note.xml` — canonical JATS (submit to journal)
   - `{OG_ID}_genome_note.docx` — Word for internal review
   - `{OG_ID}_genome_note.html` — HTML preview

### Suggested Nextflow module location

`modules/local/generate_genome_note/main.nf`

---

## Key files

| File | Purpose |
|------|---------|
| `generate_genome_note.py` | Queries DB + renders JATS XML template |
| `render_genome_note.py` | Converts JATS XML → HTML or DOCX |
| `create_yaml_template.py` | Creates per-OG YAML from DB (run first) |
| `oceangenomes_genome_note.xml.j2` | Jinja2 JATS XML template |
| `GN_Structure_template.docx` | Word template for styling |
| `stage_figures.sh` / `stage_figures.conf` | Downloads figures from Acacia |

## Template notes

- Jinja2 with `StrictUndefined` — every variable must be set or it errors
- `{{ species_name }}` is always wrapped in `<italic>` throughout the template
- `formal_voucher: false` → "specimen held at institution" (no accession number)
- `formal_voucher: true` → full museum accession sentence with `voucher_id`
- Background paragraphs: `*text*` → `<italic>text</italic>`, `[ref-id]` → citation xref

## Completed genome notes

| OG_ID | Species | Status |
|-------|---------|--------|
| OG38 | *Ostorhinchus lineolatus* | Draft complete |
| OG910 | *Choerodon rubescens* (Baldchin Groper) | Background section pending (`REPLACE`) |
