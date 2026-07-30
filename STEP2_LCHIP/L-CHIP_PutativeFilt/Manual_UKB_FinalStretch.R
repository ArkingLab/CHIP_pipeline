#!/usr/local/bin/Rscript
library(dplyr, quietly=T)
library(stringr, quietly=T)
library(data.table, quietly=T)
rm(list=ls())

fin_df <- read.delim("Putative_SemiFilteredDF.tsv", sep = "\t")

#fin_df <- fin_df[!duplicated(fin_df$JoiningVar),]

# Further gene hits -------------------------------------------------------

ATFF <- read.csv("ATFF.csv")

ATFF <- ATFF[ATFF$ATFF %in% c("Frameshift/stop gain/splice site, Nonsynonymous", "Nonsynonymous"),]


table(fin_df$Gene.refGene[!fin_df$Gene.refGene %in% ATFF$Gene])
unique(fin_df$Gene.refGene[!fin_df$Gene.refGene %in% ATFF$Gene])

#Just check that you are not missing something
unique(fin_df$Gene.refGene[!fin_df$Gene.refGene %in% ATFF$Gene])[!unique(fin_df$Gene.refGene[!fin_df$Gene.refGene %in% ATFF$Gene]) %in% read.csv("LCHIPPATH_GeneList_Bick_practical.csv")$Gene]


fin_df <- fin_df[fin_df$Gene.refGene %in% ATFF$Gene,]

fin_df$OnlyNonSyn <- fin_df$Gene.refGene %in% ATFF$Gene[ATFF$ATFF == "Nonsynonymous"]


fin_df[is.na(fin_df$NonsynOI),]
fin_df <- fin_df[!(fin_df$OnlyNonSyn & (!(grepl("^p.[A-Z][0-9]*[A-Z]$", fin_df$NonsynOI) & (!grepl("X", fin_df$NonsynOI))))),]


# Strand Bias only for some -----------------------------------------------

fin_df$SB <- str_extract(fin_df$SAMPLE, "[^:]*$")

nch <- c()
pch <- c()

for(i in 1:length(fin_df$SB)){
  
  nch <- c(nch, i)
  pch <- c(pch, fisher.test(matrix(as.numeric(unlist(str_split(fin_df$SB[i], ","))), 2))[[1]])
  print(i)
  
}

fin_df$p_SB <- pch

fin_df <- fin_df[((fin_df$FILTER != ".")|(fin_df$p_SB >= 0.05)),]

write.table(fin_df, file = "Putative_FinalFiltered.tsv", row.names = FALSE, sep = "\t")

















