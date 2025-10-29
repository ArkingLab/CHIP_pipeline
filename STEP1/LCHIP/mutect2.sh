#!/usr/bin/env bash
set -ex 

#1:input CRAM/BAM file (must be indexed)
#2:output prefix

IN=$1
OUT_PREFIX=$2

DIR="."	# must be set to the location of the downloaded cram/vcf file

ref_chr=$DIR/GRCh38_full_analysis_set_plus_decoy_hla.fa
intervals=$DIR/L-CHIP.CDS.merge.bed
ref_pon=$DIR/1000g_pon.hg38_Union_Genes.vcf.gz
ref_germ=$DIR/af-only-gnomad.hg38_Union_Genes.vcf.gz
tmp_dir=~/tmp

test -s ${ref_chr}
test -s ${IN}
mkdir -p ${tmp_dir}

if [ ! -s ${OUT_PREFIX}.filt.vcf.gz ] ; then
  gatk --java-options "-Xmx3g -XX:ParallelGCThreads=1" Mutect2 -R ${ref_chr} -I ${IN} -O ${OUT_PREFIX}.unfiltered.vcf.gz --tmp-dir ${tmp_dir} -L ${intervals} --panel-of-normals ${ref_pon} --germline-resource ${ref_germ} --f1r2-tar-gz ${OUT_PREFIX}.f1r2.tar.gz  --native-pair-hmm-threads 1
  gatk --java-options "-Xmx3g -XX:ParallelGCThreads=1" LearnReadOrientationModel -I ${OUT_PREFIX}.f1r2.tar.gz -O ${OUT_PREFIX}.read-orientation-model.tar.gz
  gatk --java-options "-Xmx3g -XX:ParallelGCThreads=1" FilterMutectCalls -R ${ref_chr} -V ${OUT_PREFIX}.unfiltered.vcf.gz -O ${OUT_PREFIX}.filt.vcf.gz -L ${intervals} --ob-priors ${OUT_PREFIX}.read-orientation-model.tar.gz
fi
