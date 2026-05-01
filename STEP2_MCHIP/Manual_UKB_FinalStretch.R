#!/usr/local/bin/Rscript
library(dplyr, quietly=T)
library(stringr, quietly=T)
library(data.table, quietly=T)
rm(list=ls())

fin_df <- read.delim("SemiFilteredDF.tsv", sep = "\t")

fin_df <- fin_df[!duplicated(fin_df$JoiningVar),]
table(read.csv("NODEL_ManualCheck.csv")$TorF)
fin_df$whitelist <- ifelse(fin_df$transcriptOI %in% (read.csv("NODEL_ManualCheck.csv")%>%
                                                         filter(TorF)%>%
                                                         pull(transcriptOI)), TRUE, fin_df$whitelist)

  
fin_df <- fin_df[fin_df$whitelist,]

# Germline variants -------------------------------------------------------

fin_df$BinomTestGermlineHet <- NA
fin_df$BinomTestGermlineHom <- NA

for(i in 1:length(fin_df$BinomTestGermlineHet)){

    fin_df$BinomTestGermlineHet[i] <- binom.test(fin_df$AltDP[i], fin_df$DP[i], p = 0.5)[[3]]
    fin_df$BinomTestGermlineHom[i] <- binom.test(fin_df$AltDP[i], fin_df$DP[i], p = 0.98)[[3]]

}

fin_df$BinomGermline <- ((fin_df$BinomTestGermlineHet >= 0.01) | (fin_df$BinomTestGermlineHom >= 0.01))



si_df <- table(fin_df$FullVar)%>%
  data.frame()


si_df <- table(fin_df$FullVar[fin_df$BinomGermline])%>%
  data.frame()%>%
  left_join(si_df, ., by = "Var1")

si_df <- si_df%>%
  filter(!is.na(Freq.y))

si_df$P <- si_df$Freq.y/si_df$Freq.x

fin_df$PotentialGermline <- fin_df$FullVar %in% (si_df%>%
                                          filter(Freq.x >= 3)%>%
                                          filter(P >= 0.8)%>%
                                          pull(Var1)%>%
                                          as.character())

sum(fin_df$PotentialGermline)

# Too many indels ---------------------------------------------------------
Indel_last_check <- table(fin_df$SampID)%>%
  data.frame()%>%
  arrange(-Freq)%>%
  filter(Freq >= 4)

Indel_last_check <- left_join(Indel_last_check, (table((fin_df$SampID[(nchar(fin_df$REF)>1|nchar(fin_df$ALT)>1)]))%>%
                                                   data.frame()%>%
                                                   arrange(-Freq)), by = "Var1")

Indel_last_check$Freq.y <- ifelse(is.na(Indel_last_check$Freq.y), 0, Indel_last_check$Freq.y)


fin_df$HighNumberOfIndels <- fin_df$SampID %in% as.character(Indel_last_check$Var1[Indel_last_check$Freq.x == Indel_last_check$Freq.y])

write.table(fin_df, file = "FinalFiltered.tsv", row.names = FALSE, sep = "\t")
















