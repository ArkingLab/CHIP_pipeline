#!/usr/local/bin/Rscript
library(dplyr, quietly=T)
library(stringr, quietly=T)
library(data.table, quietly=T)
rm(list=ls())

fin_df <- data.frame()

for(i in (list.files()[grepl("_OutU2AF1.tsv$", list.files())])){
  
  int_df <- read.delim(i)
  fin_df <- rbind(fin_df, int_df)
  
}

write.table(fin_df, file = "Final_OutU2AF1.tsv", sep = "\t", row.names = FALSE, quote = FALSE)























