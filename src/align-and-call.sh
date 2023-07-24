#!/bin/bash

ulimit -Sn 4096
genomic_reference="reference/GCA_016453205.2_ASM1645320v2_genomic.fna"
threads=12

# reference the genome if not already done
minimap_index="REF.mmi"; [ -e $minimap_index ] || minimap2 -d REF.mmi $genomic_reference

for run_type in "targets" "non-targets"; do

        samples_list=$(ls $run_type/*.fna.gz)

        bam_dir="$run_type-bam-files"; [ -d $bam_dir ] || mkdir -p $bam_dir;
        log_dir="logs"; [ -d $log_dir ] || mkdir -p $log_dir;

        # iterate through .fna.gz files
        # generate alignment, output .bam files
        for sample in $samples_list; do
                filename=$(basename "$sample" .fna.gz)
                minimap2 -a -t $threads -x asm5 $minimap_index $sample | samtools sort -o $bam_dir/$filename.bam --write-index -
        done

        # pileup the alignment
        # call variants
        # output to singular bcf file
        bcftools mpileup \
                --threads $threads \
                --output-type u \
                --fasta-ref $genomic_reference \
                $bam_dir/*.bam \
        |\
        bcftools call \
        --threads $threads \
        --ploidy 1 \
        --multiallelic-caller \
        --variants-only \
        --output-type v \
        --output $run_type.vcf
done

# ===========================
vcftools \
--vcf non-targets.vcf \
--SNPdensity 100 \
--out non-targets
# ===========================
vcftools \
--vcf targets.vcf \
--SNPdensity 100 \
--out targets
# ===========================
cat targets/*.fna.gz non-targets/*fna.gz > db_targets-and-non-targets.fna.gz
gunzip db_targets-and-non-targets.fna.gz
#makeblastdb -in db_targets-and-non-targets.fna -parse_seqids -dbtype nucl

python src/combine-snpden.py \
--targets targets.snpden \
--non-targets non-targets.snpden \
--output targets.csv