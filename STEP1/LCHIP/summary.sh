#!/usr/bin/env bash
set -euo pipefail

M="mutect2"
HP_SDIR="/MitoHPC/scripts"
HP_RDIR="/MitoHPC/RefSeq"
MUTECT2_HEADER="$HP_SDIR/$M.vcf"
HP_DP="${HP_DP:-5}"
HP_RNAME="${HP_RNAME:-hs38DH}"
VCF_DIR="${1:-.}"

cd "$VCF_DIR"
shopt -s nullglob

for VCF in ./*.lchip.filt.vcf.gz; do
  S="$(basename "$VCF" .lchip.filt.vcf.gz)"

  bcftools norm -m-any -f "$HP_RDIR/$HP_RNAME.fa" "$VCF" | \
    "$HP_SDIR/fix${M}Vcf.pl" --file "$HP_RDIR/$HP_RNAME.fa" | \
    bedtools sort -header | \
    "$HP_SDIR/filterVcf.pl" --sample "$S" --source "$M" --header "$MUTECT2_HEADER" --depth "$HP_DP" | \
    "$HP_SDIR/uniqVcf.pl" | \
    bedtools sort -header | \
    bgzip -c > "$S.00.vcf.gz"
done

find . -maxdepth 1 -name "*.00.vcf.gz" -printf "./%f\n" | sort > "./$M.lchip.vcf.list"

if [[ -s "./$M.lchip.vcf.list" ]]; then
  xargs zcat < "./$M.lchip.vcf.list" | \
    grep -v "^##sample=" | \
    "$HP_SDIR/uniq.pl" | \
    bedtools sort -header \
    > "./$M.lchip.filter.1.vcf"
fi
