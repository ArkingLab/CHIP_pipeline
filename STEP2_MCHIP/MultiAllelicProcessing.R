#!/usr/local/bin/Rscript
library(dplyr, quietly=T)
library(stringr, quietly=T)
library(data.table, quietly=T)
rm(list=ls())

fin_df <- read.csv("WITHREF_PotentialVarsMultiallelic.csv")

fin_df$DP <- str_remove(fin_df$SAMPLE, "^[^:]*:[^:]*:[^:]*:")%>%
  str_remove(., ":.*")%>%
  as.character()%>%
  as.numeric()

fin_df <- fin_df[!fin_df$DP < 20,]

fin_df$LeftAltDP <- str_remove(fin_df$AltDP, "^0,")%>%
  str_remove(., ",.*")%>%
  as.character()%>%
  as.numeric()

fin_df$RightAltDP <- str_remove(fin_df$AltDP, "^0,")%>%
  str_remove(., ".*,")%>%
  as.character()%>%
  as.numeric()

fin_df$LeftALT <- str_remove(fin_df$ALT, ",.*")
fin_df$RightALT <- str_remove(fin_df$ALT, ".*,")



fin_df$AltAltDP <- ifelse(fin_df$LeftAltDP > fin_df$RightAltDP, fin_df$RightAltDP, fin_df$LeftAltDP)
fin_df$AltALT <- ifelse(fin_df$LeftAltDP > fin_df$RightAltDP, fin_df$RightALT, fin_df$LeftALT)


fin_df$RefAltDP <- ifelse(fin_df$LeftAltDP > fin_df$RightAltDP, fin_df$LeftAltDP, fin_df$RightAltDP)
fin_df$RefALT <- ifelse(fin_df$LeftAltDP > fin_df$RightAltDP, fin_df$LeftALT, fin_df$RightALT)

fin_df <- fin_df[fin_df$AltAltDP >= 3,]

(nchar(fin_df$RefALT)>1)|(nchar(fin_df$AltALT)>1)


fin_df <- fin_df[!((nchar(fin_df$RefALT)>1|nchar(fin_df$AltALT)>1) & (fin_df$AltAltDP < 5)),]



table(paste0(fin_df$X.CHROM, "_", fin_df$POS))%>%
  sort()


nchar(fin_df$RefALT) - nchar(fin_df$REF)
nchar(fin_df$AltALT) - nchar(fin_df$REF)

sort(nchar(fin_df$RefALT))
sort(nchar(fin_df$AltALT))


write.csv(fin_df, file = "ProcessedMultiallelic.csv", row.names = FALSE)
















