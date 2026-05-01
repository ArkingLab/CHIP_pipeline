library(dplyr, quietly=T)
library(stringr, quietly=T)
library(data.table, quietly=T)
rm(list=ls())


# Extra script for germline -----------------------------------------------
#!!!recheck the script to match the names you are using
fin_df$BinomGermline <- ifelse(((fin_df$GENDER == "M") & (fin_df$X.CHROM == "X")), (fin_df$BinomTestGermlineHom >= 0.01), ((fin_df$BinomTestGermlineHet >= 0.01) | (fin_df$BinomTestGermlineHom >= 0.01)))


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

fin_df <- fin_df[(!fin_df$PotentialGermline),]


# Remove samples with too many indels or hits in general ------------------
#!!!recheck the script to match the names you are using

Indel_last_check <- table(fin_df$SampID)%>%
  data.frame()%>%
  arrange(-Freq)%>%
  filter(Freq >= 4)

Indel_last_check <- left_join(Indel_last_check, (table((fin_df$SampID[(nchar(fin_df$REF)>1|nchar(fin_df$ALT)>1)]))%>%
                                                   data.frame()%>%
                                                   arrange(-Freq)), by = "Var1")

Indel_last_check$Freq.y <- ifelse(is.na(Indel_last_check$Freq.y), 0, Indel_last_check$Freq.y)


fin_df$HighNumberOfIndels <- fin_df$SampID %in% as.character(Indel_last_check$Var1[Indel_last_check$Freq.x == Indel_last_check$Freq.y])


# Remove potential MPNs ---------------------------------------------------
#Remember to add JAK2/CALR/MPL
#!!!recheck the script to match the names you are using
#!!!recheck the script to match the names you are using
PotentialJAKMPNVect <- fin_df$SampID[(fin_df$Gene.refGene == "JAK2")]
PotentialCALRMPNVect <- fin_df$SampID[(fin_df$Gene.refGene == "CALR")]
PotentialMPLMPNVect <- fin_df$SampID[(fin_df$Gene.refGene == "MPL")]

MPN_OUT <- (PhenoDF[PhenoDF$SAMPLE %in% c(PotentialJAKMPNVect, PotentialCALRMPNVect, PotentialMPLMPNVect),])%>%
  select("SAMPLE", "HCT", "Platelet")


MPN_OUT$OutCol <- (MPN_OUT$Platelet > 450)|(MPN_OUT$HCT > 48)|(is.na(MPN_OUT$Platelet))|(is.na(MPN_OUT$HCT))

MPN_OUT$OutCol <- ifelse(is.na(MPN_OUT$OutCol), TRUE, MPN_OUT$OutCol)

OutVect <- MPN_OUT$SAMPLE[MPN_OUT$OutCol]

PhenoDF <- PhenoDF[!(PhenoDF$SAMPLE %in% OutVect),]


#Remember to keep only Full Consent or Not for Profit in the case of ARIC
#Remember to check that the sample has been sequenced and analyzed












