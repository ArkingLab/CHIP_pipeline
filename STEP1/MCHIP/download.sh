wget -N -c ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
wget -N -c ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai
wget -N -c ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.dict

wget -N -c ftp://ftp.ccb.jhu.edu/pub/dpuiu/Homo_sapiens_mito/CHIP/download.sh

wget -N -c ftp://ftp.ccb.jhu.edu/pub/dpuiu/Homo_sapiens_mito/CHIP/1000g_pon.hg38_Union_Genes.vcf.gz
wget -N -c ftp://ftp.ccb.jhu.edu/pub/dpuiu/Homo_sapiens_mito/CHIP/1000g_pon.hg38_Union_Genes.vcf.gz.tbi
wget -N -c ftp://ftp.ccb.jhu.edu/pub/dpuiu/Homo_sapiens_mito/CHIP/af-only-gnomad.hg38_Union_Genes.vcf.gz
wget -N -c ftp://ftp.ccb.jhu.edu/pub/dpuiu/Homo_sapiens_mito/CHIP/af-only-gnomad.hg38_Union_Genes.vcf.gz.tbi
wget -N -c ftp://ftp.ccb.jhu.edu/pub/dpuiu/Homo_sapiens_mito/CHIP/Union_Genes.CDS.bed
wget -N -c ftp://ftp.ccb.jhu.edu/pub/dpuiu/Homo_sapiens_mito/CHIP/README

wget -N -c https://github.com/broadinstitute/gatk/releases/download/4.2.0.0/gatk-4.2.0.0.zip ; unzip gatk-4.2.0.0.zip ; ln -s gatk-4.2.0.0/gatk .

wget -N -c ftp://ftp.ccb.jhu.edu/pub/dpuiu/Homo_sapiens_mito/CHIP/mutect2.sh
wget -N -c ftp://ftp.ccb.jhu.edu/pub/dpuiu/Homo_sapiens_mito/CHIP/NEW_UnfilteredVCF_Files_ManualAlternative.zip ; unzip NEW_UnfilteredVCF_Files_ManualAlternative.zip ; mv NEW_UnfilteredVCF_Files_ManualAlternative/*  .

chmod u+x *.sh *.R
PATH=".":$PATH
