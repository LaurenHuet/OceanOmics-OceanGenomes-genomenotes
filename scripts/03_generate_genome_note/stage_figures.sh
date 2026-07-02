#!/usr/bin/env bash
# stage_figures.sh — stage genome note figures from Acacia
#
# Edit stage_figures.conf to set OG_ID and OUTPUT_DIR, then run:
#   ./stage_figures.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/stage_figures.conf"

RA_BUCKET="pawsey0964:oceanomics-refassemblies"

mkdir -p "${OUTPUT_DIR}"
echo "Staging figures for ${OG_ID} → ${OUTPUT_DIR}"

# Helper: find a file under the OG_ID prefix by regex pattern, copy it flat into OUTPUT_DIR
fetch() {
    local pattern="$1"
    local label="$2"
    local exclude="${3:-____NOMATCH____}"
    local match
    match=$(rclone ls "${RA_BUCKET}/${OG_ID}/" 2>/dev/null \
        | awk '{print $2}' \
        | grep -E "${pattern}" \
        | grep -Ev "${exclude}" \
        | sort | tail -1)
    if [[ -n "${match}" ]]; then
        rclone copyto "${RA_BUCKET}/${OG_ID}/${match}" \
            "${OUTPUT_DIR}/$(basename "${match}")" \
            --local-no-set-modtime
        echo "  [OK]      ${label}: $(basename "${match}")"
    else
        echo "  [MISSING] ${label}"
    fi
}

fetch "genomescope/${OG_ID}_genomescope_linear_plot\.png$" \
      "GenomeScope linear plot"

fetch "pretext_snapshots/hap1/.*\.3\.curated\.hap1\.pretext_snapshotFullMap\.png$" \
      "Hi-C contact map (hap1 FullMap)" \
      "SUPER|unloc"

fetch "merqury/.*\.3\.curated\.spectra-asm\.ln\.png$" \
      "Merqury spectra-asm (ln, curated)"

echo ""
echo "Done. Files in ${OUTPUT_DIR}:"
ls -lh "${OUTPUT_DIR}"/*.png 2>/dev/null | awk '{print "  "$NF, "("$5")"}' \
    || echo "  (none)"
