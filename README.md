# CHIP pipeline 

## Setup

 This pipeline is for calling LCHIP and MCHIP from WES/WGS CRAM files. The pipeline has 3 steps in total:

Before STEP1, it is recommended to download MitoHPC (https://github.com/ArkingLab/MitoHPC) and run `init.sh` so the required helper scripts and reference setup are available for this pipeline.

In order to run the full CHIP pipeline, please apply the following changes to scripts within MitoHPC:
1. fixmutect2Vcf.pl: comment out line 46 and 47
2. filterVcf.pl: comment out line 91, add the following immediately after:
```
    if($F[7] eq "." || $F[7] eq "") {
        $F[7] = "SM=$opt{sample}";
    } else {
        $F[7] .= ";SM=$opt{sample}";
    }
```
and comment out line 115 and 116

## Running full MCHIP and LCHIP pipeline

### STEP1

Step1 is running Mutect2 for CHIP genes, to call from CHIP regions

- For this step, first run the download script to download all files, then run `mutect2.sh` script to call CHIP variants.
- After `mutect2.sh` has finished for the cohort, run `summary.sh` to pull the per-sample calls into the summary VCF used as input for STEP2.
- Run both MCHIP and LCHIP

### STEP2
#### STEP2_LCHIP

- Copy the `mutect2.lchip.filter.1.vcf` file from step 1 to this directory.
- Remove header: `tail -n +32 mutect2.lchip.filter.1.vcf > mutect2.lchip.filternh.1.vcf`
- Insert a header row: `(echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE"; cat mutect2.lchip.filternh.1.vcf) > mutect2.lchip.filterh.1.vcf`
- Copy and rename the input file: `mv mutect2.lchip.filterh.1.vcf Input.vcf`
- Make sure the following packages are installed under command line R:
```
library("dplyr")
library("stringr")
library("data.table")
```
##### PathFilt

Download annovar:
```
wget www.openbioinformatics.org/annovar/download/0wgxR2rIVP/annovar.latest.tar.gz
tar -zxvf annovar.latest.tar.gz
cd annovar
chmod +x *.pl
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar refGene humandb/
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar gnomad211_exome humandb/
```
- Change permission for all files: `chmod +x *`
- Copy file over: `cp ../Input.vcf ./`
- Run START bash script first: `. ./Manual_START_UKB_Filter.sh`
- Then the Rscript: `Rscript Manual_UKB_FinalStretch.R`

##### PutativeFilt

```
cp ../Input.vcf ./
chmod +x *
```
For putativeFilt, the script needs sex information for each sample, where column `id` is sample id and column `sex` is `men/women`, and save as `idsex.txt`.

- Run START script first: `. ./Manual_START_UKB_Filter.sh`
- Then the Rscript: `Rscript Manual_UKB_FinalStretch.R`

##### Merger

- Copy the two finalfiltered.csv file generated above and the idsex file over.
- Change the `TNOP` in `Rscript Manual_UKB_FinalStretch.R` to sample size.
- Run: `Rscript MergerOfFinalDF.R`
- `FinalFiltered.tsv` is the final output to be used in STEP3.

#### STEP2_MCHIP

`chmod +x *`
```
cp ../STEP1/mchip_sample/mutect2.mchip.filter.1.vcf ./
tail -n +32 mutect2.mchip.filter.1.vcf > mutect2.mchip.filternh.1.vcf
mv mutect2.mchip.filternh.1.vcf Input.vcf
(echo -e "X.CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE"; cat Input.vcf) > Input.vcf
```
- Run: `. ./Manual_START_UKB_Filter.sh`
- Then: `Rscript Manual_UKB_FinalStretch.R`

Now copy the MCHIP and LCHIP `FinalFiltered.csv` file to STEP3, along with the `idsex.txt` file.

### STEP3

Run `Final_summary.R` and generate files for MCHIP and LCHIP

