#!/usr/bin/env bash
# shellcheck shell=bash
set -euo pipefail

#bash 04_get_busco_from_config.sh ../genomenotes_pipeline.conf

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <config.conf>" >&2
  exit 1
fi

CONFIG_FILE="$1"
[[ -f "$CONFIG_FILE" ]] || { echo "Config not found: $CONFIG_FILE" >&2; exit 1; }

# Load config (bash-compatible KEY=VALUE)
set -a
# shellcheck disable=SC1090
source "$CONFIG_FILE"
set +a

# expand {user} placeholders in paths
USER_REAL="${USER:-$(whoami)}"
expand_user() { printf '%s' "${1//\{user\}/$USER_REAL}"; }

SAMPLESHEET="$(expand_user "${SAMPLESHEET}")"
RCLONE_FLAGS="${RCLONE_FLAGS:-}"

[[ -f "$SAMPLESHEET" ]] || { echo "Samplesheet not found: $SAMPLESHEET" >&2; exit 1; }

# Helper: trim CR/LF and leading/trailing whitespace
clean() { printf '%s' "$1" | tr -d '\r' | xargs; }


BUSCO_BUCKET="${BUSCO_BUCKET:?Missing BUSCO_BUCKET in config}"

BUSCO_REQUIRED_SUBSTR_1="${BUSCO_REQUIRED_SUBSTR_1:-busco}"
BUSCO_REQUIRED_SUBSTR_2="${BUSCO_REQUIRED_SUBSTR_2:-hap1.chr_level}"
BUSCO_REQUIRED_SUBSTR_3="${BUSCO_REQUIRED_SUBSTR_3:-full_table.tsv}"

tmp="$(mktemp)"
trap 'rm -f "$tmp" "${tmp}.matched" 2>/dev/null || true' EXIT

declare -A BUSCO_DIR
samples=()

# map sample -> busco_genes dir (col1 -> col5)
while IFS=$'\t' read -r s dir; do
  s="$(clean "$s")"
  dir="$(clean "$dir")"
  [[ -z "$s" ]] && continue
  BUSCO_DIR["$s"]="$(expand_user "$dir")"
  samples+=("$s")
done < <(awk -F, 'NR>1{ gsub(/\r/,"",$1); gsub(/\r/,"",$5); print $1 "\t" $5 }' "$SAMPLESHEET")

[[ ${#samples[@]} -gt 0 ]] || { echo "No samples found in $SAMPLESHEET" >&2; exit 1; }

for s in "${samples[@]}"; do
  target="${BUSCO_DIR[$s]:-}"
  [[ -n "$target" ]] || { echo "Error: no busco_genes entry for ${s} in samplesheet" >&2; exit 1; }

  mkdir -p "$target"

  rclone lsf --recursive ${RCLONE_FLAGS:+$RCLONE_FLAGS} "${BUSCO_BUCKET}/${s}" > "$tmp" 2>/dev/null || true

  grep -i "$BUSCO_REQUIRED_SUBSTR_1" "$tmp" | \
    grep -i "$BUSCO_REQUIRED_SUBSTR_2" | \
    grep -i "$BUSCO_REQUIRED_SUBSTR_3" > "${tmp}.matched" || true

  if [[ ! -s "${tmp}.matched" ]]; then
    echo "Error: no matching BUSCO table found for ${s} on remote ${BUSCO_BUCKET}/${s}" >&2
    exit 1
  fi

  while IFS= read -r relpath; do
    echo "Copying ${BUSCO_BUCKET}/${s}/${relpath} -> ${target}/"
    rclone copy ${RCLONE_FLAGS:+$RCLONE_FLAGS} "${BUSCO_BUCKET}/${s}/${relpath}" "${target}/"
  done < "${tmp}.matched"
done

echo "Done."
