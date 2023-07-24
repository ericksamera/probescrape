#!/bin/bash

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