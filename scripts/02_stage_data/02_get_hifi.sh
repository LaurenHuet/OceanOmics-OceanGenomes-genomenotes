#!/usr/bin/env bash
set -euo pipefail

SAMPLESHEET="/scratch/pawsey0964/lhuet/genomenotes/OceanOmics-OceanGenomes-genomenotes/assets/samplesheet.csv"
s3_bucket="pawsey0964:oceanomics-filtered-reads"
tmp_list="$(mktemp)"

[[ -f "$SAMPLESHEET" ]] || { echo "Samplesheet not found: $SAMPLESHEET" >&2; exit 1; }

declare -A HIFI_DIR_MAP
samples=()

while IFS=$'\t' read -r sample hifi_dir; do
  sample="$(printf '%s' "$sample" | tr -d '\r' | tr -cd '[:print:]' | xargs)"
  hifi_dir="$(printf '%s' "$hifi_dir" | tr -d '\r' | xargs)"
  [[ -z "$sample" ]] && continue
  HIFI_DIR_MAP["$sample"]="$hifi_dir"
  samples+=("$sample")
done < <(awk -F, 'NR>1 { gsub(/\r/,"",$1); gsub(/\r/,"",$2); print $1 "\t" $2 }' "$SAMPLESHEET")

[[ ${#samples[@]} -gt 0 ]] || { echo "No samples found in $SAMPLESHEET" >&2; exit 1; }

# only include remote entries that contain "hifi_reads" (case-insensitive)
include=()
for s in "${samples[@]}"; do
  include+=( --include "*${s}*hifi_reads*" )
done

rclone ls "$s3_bucket" "${include[@]}" > "$tmp_list"

while IFS= read -r line || [[ -n "$line" ]]; do
  path="$(printf '%s' "$line" | sed -E 's/^[[:space:]]*[0-9]+[[:space:]]+//; s/\r$//')"
  # ensure it's a HiFi read file
  if ! printf '%s' "$path" | grep -qi 'hifi_reads'; then
    continue
  fi

  og="$(printf '%s' "$path" | grep -o -m1 -E 'OG[0-9]+' || true)"

  target=""
  if [[ -n "$og" && -n "${HIFI_DIR_MAP[$og]:-}" ]]; then
    target="${HIFI_DIR_MAP[$og]}"
  else
    for s in "${samples[@]}"; do
      if [[ "$path" == *"$s"* ]]; then
        target="${HIFI_DIR_MAP[$s]}"
        og="$s"
        break
      fi
    done
  fi

  if [[ -z "$target" ]]; then
    user="${USER:-$(whoami)}"
    target="/scratch/pawsey0964/${user}/genomenotes/${og:-unknown}/hifi"
    echo "Warning: no hifi_dir for remote path; falling back to ${target}" >&2
  fi

  mkdir -p "$target"
  rclone copy "${s3_bucket}/${path}" "${target}/"
done < "$tmp_list"

rm -f "$tmp_list"
echo "Done."
```// filepath: /scratch/pawsey0964/lhuet/genomenotes/OceanOmics-OceanGenomes-genomenotes/scripts/02_stage_data/02_get_hifi.sh
#!/usr/bin/env bash
set -euo pipefail

SAMPLESHEET="/scratch/pawsey0964/lhuet/genomenotes/OceanOmics-OceanGenomes-genomenotes/assets/samplesheet.csv"
s3_bucket="pawsey0964:oceanomics-filtered-reads"
tmp_list="$(mktemp)"

[[ -f "$SAMPLESHEET" ]] || { echo "Samplesheet not found: $SAMPLESHEET" >&2; exit 1; }

declare -A HIFI_DIR_MAP
samples=()

while IFS=$'\t' read -r sample hifi_dir; do
  sample="$(printf '%s' "$sample" | tr -d '\r' | tr -cd '[:print:]' | xargs)"
  hifi_dir="$(printf '%s' "$hifi_dir" | tr -d '\r' | xargs)"
  [[ -z "$sample" ]] && continue
  HIFI_DIR_MAP["$sample"]="$hifi_dir"
  samples+=("$sample")
done < <(awk -F, 'NR>1 { gsub(/\r/,"",$1); gsub(/\r/,"",$2); print $1 "\t" $2 }' "$SAMPLESHEET")

[[ ${#samples[@]} -gt 0 ]] || { echo "No samples found in $SAMPLESHEET" >&2; exit 1; }

# only include remote entries that contain "hifi_reads" (case-insensitive)
include=()
for s in "${samples[@]}"; do
  include+=( --include "*${s}*hifi_reads*" )
done

rclone ls "$s3_bucket" "${include[@]}" > "$tmp_list"

while IFS= read -r line || [[ -n "$line" ]]; do
  path="$(printf '%s' "$line" | sed -E 's/^[[:space:]]*[0-9]+[[:space:]]+//; s/\r$//')"
  # ensure it's a HiFi read file
  if ! printf '%s' "$path" | grep -qi 'hifi_reads'; then
    continue
  fi

  og="$(printf '%s' "$path" | grep -o -m1 -E 'OG[0-9]+' || true)"

  target=""
  if [[ -n "$og" && -n "${HIFI_DIR_MAP[$og]:-}" ]]; then
    target="${HIFI_DIR_MAP[$og]}"
  else
    for s in "${samples[@]}"; do
      if [[ "$path" == *"$s"* ]]; then
        target="${HIFI_DIR_MAP[$s]}"
        og="$s"
        break
      fi
    done
  fi

  if [[ -z "$target" ]]; then
    user="${USER:-$(whoami)}"
    target="/scratch/pawsey0964/${user}/genomenotes/${og:-unknown}/hifi"
    echo "Warning: no hifi_dir for remote path; falling back to ${target}" >&2
  fi

  mkdir -p "$target"
  rclone copy "${s3_bucket}/${path}" "${target}/"
done < "$tmp_list"

rm -f "$tmp_list"
echo "Done."