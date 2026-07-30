#!/usr/local/bin/Rscript
library(dplyr, quietly=T)
library(stringr, quietly=T)
library(data.table, quietly=T)
rm(list=ls())

Putative_DF <- read.delim("Putative_FinalFiltered.tsv", sep = "\t")

Putative_DF <- Putative_DF%>%
  select(-sex)
Pathogenic_DF <- read.delim("Pathogenic_FinalFiltered.tsv", sep = "\t")

VectToMod <- read.csv("PathogenicNODEL_ManualCheck.csv")%>%
  filter(TorF == TRUE)%>%
  pull(transcriptOI)%>%
  as.character()

Pathogenic_DF$TruncComb <- ((Pathogenic_DF$transcriptOI %in% VectToMod)|Pathogenic_DF$wl.splice|Pathogenic_DF$wl.lof)

Pathogenic_DF <- Pathogenic_DF[!((Pathogenic_DF$Gene.refGene == "NOTCH1")&(Pathogenic_DF$TruncComb)&(!grepl(":exon34:", Pathogenic_DF$transcriptOI))),]
Pathogenic_DF <- Pathogenic_DF[!((Pathogenic_DF$Gene.refGene == "NOTCH2")&(Pathogenic_DF$TruncComb)&(!grepl(":exon34:", Pathogenic_DF$transcriptOI))),]
Pathogenic_DF <- Pathogenic_DF[!((Pathogenic_DF$Gene.refGene == "CXCR4")&(Pathogenic_DF$TruncComb)&(!grepl(":exon2:", Pathogenic_DF$transcriptOI))),]
#Pathogenic_DF <- Pathogenic_DF[!((Pathogenic_DF$Gene.refGene == "IRF8")&(Pathogenic_DF$TruncComb)&(!grepl(":any:", Pathogenic_DF$transcriptOI))),]#any
#Pathogenic_DF <- Pathogenic_DF[!((Pathogenic_DF$Gene.refGene == "CCND3")&(Pathogenic_DF$TruncComb)&(!grepl(":any:", Pathogenic_DF$transcriptOI))),]#any



table(nchar(Putative_DF$SampID))
table(nchar(Pathogenic_DF$SampID))
fin_df <- bind_rows(Pathogenic_DF, Putative_DF)

fin_df <- fin_df[!grepl("weak_evidence", fin_df$FILTER),]

GenderDF <- fread("idsex.txt")%>%#running the gender again to modify everything all at once (it does not apply the gender adjustment to the adjusted value, no worries)
  select(id, sex)

GenderDF$id <- as.character(GenderDF$id)

fin_df$SampID <- as.character(fin_df$SampID)

fin_df <- left_join(fin_df, GenderDF, by = c("SampID" = "id"))
table(fin_df$X.CHROM)
fin_df$VAF_adj <- ifelse(fin_df$X.CHROM == "X" & fin_df$sex == "men", fin_df$newVAF/2, fin_df$newVAF)
table(fin_df$X.CHROM[fin_df$VAF_adj != fin_df$newVAF])

#fin_df <- fin_df[!duplicated(fin_df$JoiningVar),]

fin_df <- fin_df[!(fin_df$VAF_adj < 2),]

fin_df <- fin_df[!((nchar(fin_df$REF)>1|nchar(fin_df$ALT)>1) & (fin_df$VAF_adj < 10)),]

# Too frequent ------------------------------------------------------------

TFDFTR <- table(fin_df$FullVar)%>%
  data.frame()

MisPath <- read.csv("LCHIPPATH_CHIP_missense_vars_Bick_practical.csv")

fin_df$RepPathVar <- paste0(fin_df$Gene.refGene, " ", fin_df$NonsynOI) %in% paste0(MisPath$Gene, " p.", MisPath$AAChange)

TNOP <- 12345 #Set this to the total number of individuals assessed
fin_df$PTF <- fin_df$FullVar %in% as.character(TFDFTR$Var1[(TFDFTR$Freq/TNOP*100)>1])

fin_df <- fin_df[!((fin_df$RepPathVar == FALSE)&fin_df$PTF),]

# Too large VAF -----------------------------------------------------------

fin_df <- fin_df[!((fin_df$VAF_adj>=35)&(fin_df$RepPathVar==FALSE)),]

# Germline variants -------------------------------------------------------

#This section is not applied because the VAF_adj >= 35% was already applied; it is more of a check
fin_df$BinomTestGermlineHet <- NA
fin_df$BinomTestGermlineHom <- NA

fin_df = fin_df %>% filter(!is.na(AltDP))

for(i in 1:length(fin_df$BinomTestGermlineHet)){

  fin_df$BinomTestGermlineHet[i] <- binom.test(fin_df$AltDP[i], fin_df$DP[i], p = 0.5)[[3]]
  fin_df$BinomTestGermlineHom[i] <- binom.test(fin_df$AltDP[i], fin_df$DP[i], p = 0.98)[[3]]

}

fin_df$BinomGermline <- ifelse(((fin_df$sex == "men") & (fin_df$X.CHROM == "X")), (fin_df$BinomTestGermlineHom >= 0.01), ((fin_df$BinomTestGermlineHet >= 0.01) | (fin_df$BinomTestGermlineHom >= 0.01)))


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
table(fin_df$transcriptOI[fin_df$PotentialGermline])

# Too many gene level hits ---------------------------------------------------------

DFTRM <- table(fin_df$Gene.refGene)%>%
  data.frame()%>%
  mutate(P = Freq/TNOP*100)


FLOGTR_DF <- data.frame()

for(k in as.character(DFTRM$Var1[DFTRM$P>0.05])){

  GeneIntDf <- table(fin_df$FullVar[fin_df$Gene.refGene == k])%>%
    data.frame()%>%
    filter(Freq>1)
  
  if(length(as.character(GeneIntDf$Var1))==0){
    next
  }
  
  GeneIntDf$OutVar <- (GeneIntDf$Freq>quantile(GeneIntDf$Freq)[3][[1]])
  GeneIntDf$Gene <- k
  
  FLOGTR_DF <- rbind(FLOGTR_DF, GeneIntDf)
  
  print(k)
  
}
sum(fin_df$RepPathVar)
fin_df <- fin_df[(!(fin_df$FullVar %in% as.character(FLOGTR_DF$Var1[(FLOGTR_DF$Freq>10)&FLOGTR_DF$OutVar])))|fin_df$RepPathVar,]
sum(fin_df$RepPathVar)

length(unique(fin_df$SampID))/TNOP*100
sort(table(fin_df$Gene.refGene))

write.table(fin_df, file = "FinalFiltered.tsv", row.names = FALSE, sep = "\t")































