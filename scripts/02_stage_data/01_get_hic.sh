#!/usr/bin/env bash
set -euo pipefail

SAMPLESHEET="/scratch/pawsey0964/lhuet/genomenotes/OceanOmics-OceanGenomes-genomenotes/assets/samplesheet.csv"
s3_bucket="s3:oceanomics/OceanGenomes/illumina-hic"
tmp_list="$(mktemp)"

[[ -f "$SAMPLESHEET" ]] || { echo "Samplesheet not found: $SAMPLESHEET" >&2; exit 1; }

declare -A HIC_DIR_MAP
samples=()

# build map: sample -> hic_dir (sanitize sample names)
while IFS=$'\t' read -r sample hic_dir; do
  sample="$(printf '%s' "$sample" | tr -d '\r' | tr -cd '[:print:]' | xargs)"
  hic_dir="$(printf '%s' "$hic_dir" | tr -d '\r' | xargs)"
  [[ -z "$sample" ]] && continue
  HIC_DIR_MAP["$sample"]="$hic_dir"
  samples+=("$sample")
done < <(awk -F, 'NR>1 { gsub(/\r/,"",$1); gsub(/\r/,"",$3); print $1 "\t" $3 }' "$SAMPLESHEET")

[[ ${#samples[@]} -gt 0 ]] || { echo "No samples found in $SAMPLESHEET" >&2; exit 1; }

# create include filters and list remote files
include=()
for s in "${samples[@]}"; do include+=( --include "${s}*" ); done
rclone ls "$s3_bucket" "${include[@]}" > "$tmp_list"

# copy each remote path into matching hic_dir
while IFS= read -r line || [[ -n "$line" ]]; do
  path="$(printf '%s' "$line" | sed -E 's/^[[:space:]]*[0-9]+[[:space:]]+//; s/\r$//')"
  og="$(printf '%s' "$path" | grep -o -m1 -E 'OG[0-9]+' || true)"

  target=""
  if [[ -n "$og" && -n "${HIC_DIR_MAP[$og]:-}" ]]; then
    target="${HIC_DIR_MAP[$og]}"
  else
    for s in "${samples[@]}"; do
      if [[ "$path" == *"$s"* ]]; then
        target="${HIC_DIR_MAP[$s]}"
        og="$s"
        break
      fi
    done
  fi

  if [[ -z "$target" ]]; then
    user="${USER:-$(whoami)}"
    target="/scratch/pawsey0964/${user}/genomenotes/${og:-unknown}/hic"
    echo "Warning: no hic_dir for remote path; falling back to ${target}" >&2
  fi

  mkdir -p "$target"
  echo "Copying ${s3_bucket}/${path} -> ${target}/"
  rclone copy "${s3_bucket}/${path}" "${target}/"
done < "$tmp_list"

rm -f "$tmp_list"
echo "Done."
# ...existing code...
```// filepath: /scratch/pawsey0964/lhuet/genomenotes/OceanOmics-OceanGenomes-genomenotes/scripts/02_stage_data/01_get_hic.sh
# ...existing code...
#!/usr/bin/env bash
set -euo pipefail

SAMPLESHEET="/scratch/pawsey0964/lhuet/genomenotes/OceanOmics-OceanGenomes-genomenotes/assets/samplesheet.csv"
s3_bucket="s3:oceanomics/OceanGenomes/illumina-hic"
tmp_list="$(mktemp)"

[[ -f "$SAMPLESHEET" ]] || { echo "Samplesheet not found: $SAMPLESHEET" >&2; exit 1; }

declare -A HIC_DIR_MAP
samples=()

# build map: sample -> hic_dir (sanitize sample names)
while IFS=$'\t' read -r sample hic_dir; do
  sample="$(printf '%s' "$sample" | tr -d '\r' | tr -cd '[:print:]' | xargs)"
  hic_dir="$(printf '%s' "$hic_dir" | tr -d '\r' | xargs)"
  [[ -z "$sample" ]] && continue
  HIC_DIR_MAP["$sample"]="$hic_dir"
  samples+=("$sample")
done < <(awk -F, 'NR>1 { gsub(/\r/,"",$1); gsub(/\r/,"",$3); print $1 "\t" $3 }' "$SAMPLESHEET")

[[ ${#samples[@]} -gt 0 ]] || { echo "No samples found in $SAMPLESHEET" >&2; exit 1; }

# create include filters and list remote files
include=()
for s in "${samples[@]}"; do include+=( --include "${s}*" ); done
rclone ls "$s3_bucket" "${include[@]}" > "$tmp_list"

# copy each remote path into matching hic_dir
while IFS= read -r line || [[ -n "$line" ]]; do
  path="$(printf '%s' "$line" | sed -E 's/^[[:space:]]*[0-9]+[[:space:]]+//; s/\r$//')"
  og="$(printf '%s' "$path" | grep -o -m1 -E 'OG[0-9]+' || true)"

  target=""
  if [[ -n "$og" && -n "${HIC_DIR_MAP[$og]:-}" ]]; then
    target="${HIC_DIR_MAP[$og]}"
  else
    for s in "${samples[@]}"; do
      if [[ "$path" == *"$s"* ]]; then
        target="${HIC_DIR_MAP[$s]}"
        og="$s"
        break
      fi
    done
  fi

  if [[ -z "$target" ]]; then
    user="${USER:-$(whoami)}"
    target="/scratch/pawsey0964/${user}/genomenotes/${og:-unknown}/hic"
    echo "Warning: no hic_dir for remote path; falling back to ${target}" >&2
  fi

  mkdir -p "$target"
  echo "Copying ${s3_bucket}/${path} -> ${target}/"
  rclone copy "${s3_bucket}/${path}" "${target}/"
done < "$tmp_list"

rm -f "$tmp_list"
echo "Done."
