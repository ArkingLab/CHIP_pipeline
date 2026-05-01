#!/usr/local/bin/Rscript
library(dplyr, quietly=T)
library(stringr, quietly=T)
library(data.table, quietly=T)
rm(list=ls())

fin_df <- read.delim("Pathogenic_SemiFilteredDF.tsv", sep = "\t")

fin_df <- fin_df[!duplicated(fin_df$JoiningVar),]

fin_df$whitelist <- ifelse(fin_df$transcriptOI %in% (read.csv("NODEL_ManualCheck.csv")%>%
                                                         filter(TorF)%>%
                                                         pull(transcriptOI)), TRUE, fin_df$whitelist)

  
fin_df <- fin_df[fin_df$whitelist,]

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
sum(fin_df$FILTER == ".")
fin_df <- fin_df[((fin_df$FILTER != ".")|(fin_df$p_SB >= 0.05)),]


write.table(fin_df, file = "Pathogenic_FinalFiltered.tsv", row.names = FALSE, sep = "\t")

















