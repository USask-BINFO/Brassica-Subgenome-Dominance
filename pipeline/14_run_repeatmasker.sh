#!/bin/bash
# Re-run RepeatMasker on nigra with the SAME TE library as the other five genomes.
#
# Why: nigra's existing LTR calls use a different library (motifs ltr-1_family-*, not in
# ABC_19k_TElib.fa) AND appear to count whole elements rather than terminal-repeat segments:
# 133,665 features averaging 1,148 bp covering 153 Mb (31% of the genome), against rapa 10.2 Mb
# (3.4%) and oleracea 7.6 Mb (1.7%). Neither the library nor the feature definition matches, so
# nigra cannot be compared to the other five as it stands.
#
# Design: VALIDATE FIRST. Run rapa chromosome A01 through this pipeline and check it reproduces
# Sampath's calls for A01 (970 features, 332,379 LTR bp). Only if that matches do the nigra
# numbers mean anything, because it shows the settings are equivalent. Then run nigra B1-B8.
#
# Subset rule, taken from Sampath's files: keep features whose motif name ends in _LTR. That
# selects the long-terminal-repeat segments and excludes _INT internal regions. All 17,074
# features in Sampath's rapa .LTR.gff follow this rule, and all 1,371 distinct motifs it uses
# are present in ABC_19k_TElib.fa (0 missing), which is how the library was confirmed.
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
  ltr_subset "$G" "$OUT/rapa_A01.mine.LTR.gff"
  ;;
nigra)
  run_rm "$PWD/06_repeatmasker/input/nigra_chr.fa" "$PWD/06_repeatmasker/nigra"
  G=$(ls "$PWD"/06_repeatmasker/nigra/*.out.gff | head -1)
  ltr_subset "$G" "$PWD/06_repeatmasker/nigra/nigra.19K.LTR.gff"
  ;;
diploids)
  # All three diploids through ONE pipeline, so every ratio between them is internally
  # consistent and no cross-pipeline correction is needed. my oleracea/rapa can then be
  # checked against Sampath's 1.45 as a ratio-level validation.
  for SP in rapa oleracea nigra; do
    D="$PWD/06_repeatmasker/diploids/$SP"
    if [ -s "$D/$SP.19K.LTR.gff" ]; then echo "[skip] $SP already done"; continue; fi
    run_rm "$PWD/06_repeatmasker/input/${SP}_chr.fa" "$D"
    G=$(ls "$D"/*.out.gff | head -1)
    ltr_subset "$G" "$D/$SP.19K.LTR.gff"
  done
  echo "[$(date +%H:%M:%S)] all three diploids complete"
  ;;
*) echo "usage: $0 [validate|nigra]"; exit 1;;
esac
