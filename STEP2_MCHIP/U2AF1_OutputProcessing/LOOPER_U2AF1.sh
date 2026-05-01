#!/bin/bash

for i in `ls -d */`;
do

echo $i|sed 's/.$//' > ForR.tmp;

ls -l $i|awk '{if($5 < 1) print $5","$9}' > EmptyCheck.csv;

./LOOPING_PuttingU2AF1Together.R;

rm ForR.tmp;

done

for j in `ls|grep "_ZeroFiles.tsv$"`;
do

cat $j >> TheFinalZeroFiles.tsv;

done

./CollapserOfFiles.R;
