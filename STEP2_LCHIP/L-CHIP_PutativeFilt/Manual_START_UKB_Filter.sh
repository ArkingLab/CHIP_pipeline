#!/bin/bash

annovarpath="/annovar"

#echo $1 > PathToFOI.tmp;
#cat $1|awk '/CHROM/{print NR}' > LOC.txt;

Rscript ./Manual_ConversionToAnnovar_UKB.R;

$annovarpath/table_annovar.pl fin_df_Genome.input.txt $annovarpath/humandb/ -buildver hg38 -out PutativeOut -remove -protocol refGene,gnomad211_exome -operation g,f -nastring . -csvout -polish

Rscript ./Manual_UKB_WhiteList.R;

rm PutativeOut.hg38_multianno.csv;
#rm PathToFOI.tmp;
#rm LOC.txt;
rm fin_df_Genome.input.txt;
