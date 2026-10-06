#!/bin/bash
# RepeatMasker with the common Brassica TE library; keeps the _LTR motif segments.
set -euo pipefail
cd "$(dirname "$0")/.."
ENV="$PWD/envs/rm"
export PATH="$ENV/bin:$PATH"
LIB="$PWD/06_repeatmasker/input/ABC_19k_TElib.fa"
PA="${PA:-16}"          # RepeatMasker parallel jobs; each rmblast job uses several threads
RMOPTS="${RMOPTS:-}"    # extra RepeatMasker flags, e.g. -s for higher sensitivity
STAGE="${1:-validate}"

run_rm () {  # $1 = query fasta, $2 = output dir
  mkdir -p "$2"
  echo "[$(date +%H:%M:%S)] RepeatMasker -pa $PA -lib ABC_19k -gff $RMOPTS  $(basename "$1")"
  RepeatMasker -pa "$PA" -lib "$LIB" -gff -xsmall $RMOPTS -dir "$2" "$1"
  echo "[$(date +%H:%M:%S)] done -> $2"
}

ltr_subset () {  # $1 = .out.gff produced by RepeatMasker -> $2 = _LTR-only gff
  awk -F'\t' '$0 !~ /^#/ && $9 ~ /Motif:[^"]*_LTR"/' "$1" > "$2"
  echo "  $(wc -l < "$2") _LTR features -> $2"
}

case "$STAGE" in
validate)
  OUT="${OUTDIR:-$PWD/06_repeatmasker/validate_rapa}"
  run_rm "$PWD/06_repeatmasker/input/rapa_A01.fa" "$OUT"
  G=$(ls "$OUT"/*.out.gff | head -1)
  ltr_subset "$G" "$OUT/rapa_A01.LTR.gff"
  ;;
nigra)
  run_rm "$PWD/06_repeatmasker/input/nigra_chr.fa" "$PWD/06_repeatmasker/nigra"
  G=$(ls "$PWD"/06_repeatmasker/nigra/*.out.gff | head -1)
  ltr_subset "$G" "$PWD/06_repeatmasker/nigra/nigra.19K.LTR.gff"
  ;;
diploids)
  # All three diploids through one pipeline, so every ratio between them is internally
  # consistent and no cross-pipeline correction is needed.
  for SP in rapa oleracea nigra; do
    D="$PWD/06_repeatmasker/diploids/$SP"
    if [ -s "$D/$SP.19K.LTR.gff" ]; then echo "[skip] $SP already done"; continue; fi
    run_rm "$PWD/06_repeatmasker/input/${SP}_chr.fa" "$D"
    G=$(ls "$D"/*.out.gff | head -1)
    ltr_subset "$G" "$D/$SP.19K.LTR.gff"
  done
  echo "[$(date +%H:%M:%S)] all three diploids complete"
  ;;
*) echo "usage: $0 [validate|nigra|diploids]"; exit 1;;
esac
