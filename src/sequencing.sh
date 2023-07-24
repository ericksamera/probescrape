#!/bin/bash

genomic_reference="reference/GCA_016453205.2_ASM1645320v2_genomic.fna"

# ===========================
vcftools \
--vcf non-targets.vcf \
--SNPdensity 1000 \
--out seq-non-targets
# ===========================
vcftools \
--vcf targets.vcf \
--SNPdensity 1000 \
--out seq-targets
# ===========================

python src/combine-snpden.py \
--targets seq-targets.snpden \
--non-targets seq-non-targets.snpden \
--output seq-targets.csv \
--bin_size 1000

head -n 100 seq-targets.csv > short-seq-targets.csv

python src/seqscrape.py \
--targets 'short-seq-targets.csv' \
--reference $genomic_reference