#!/usr/local/bin/Rscript
library(dplyr, quietly=T)
library(stringr, quietly=T)
library(data.table, quietly=T)
rm(list=ls())

BadVect2 <- read.csv("EmptyCheck.csv", header = FALSE)
BadVect2 <- BadVect2[BadVect2$V2 != "",]
BadVect2 <- BadVect2[!is.na(BadVect2$V2),]

print("The following is the distribution of bits from the removed files")
print(table(BadVect2$V1))

ForR.tmp <- as.character(read.delim("ForR.tmp", header = FALSE)$V1)

PRUA <- data.frame()

FullLength <- length(list.files(ForR.tmp))
countU2AF1 <- 0
for(i in list.files(ForR.tmp)){
  
  countU2AF1 <- countU2AF1 + 1
  if(i %in% BadVect2$V2){
    next
  }
  
  int_df <- read.delim(paste0(ForR.tmp, "/", i), sep = "\t", header = FALSE)
  colnames(int_df) <- c("X.CHROM", "POS", "REF", "ALT", "aa_change", "FULL_DP", "ALT_DP")
  int_df$SAMPLE <- i
  PRUA <- rbind(PRUA, int_df)
  
  #if((round(countU2AF1/FullLength*100, 0) == 25)|(round(countU2AF1/FullLength*100, 0) == 50)|(round(countU2AF1/FullLength*100, 0) == 75)){
    print(countU2AF1/FullLength*100)  
  #}
  
}



PRUA <- PRUA[PRUA$FULL_DP>0,]

PRUA$CombID <- paste(PRUA$SAMPLE, PRUA$aa_change,  sep = ":")

sPRUA <- PRUA%>%
  group_by(CombID)%>%
  summarise(sFULL_DP = sum(FULL_DP), sALT_DP = sum(ALT_DP),
            mFULL_DP = mean(FULL_DP), mALT_DP = mean(ALT_DP))%>%
  ungroup()

sPRUA$FullVar <- str_remove(sPRUA$CombID, ".*:")
sPRUA$SAMPLE <- str_remove(sPRUA$CombID, ":.*")
sPRUA$LargeQ <- ((sPRUA$FullVar %in% "R156Q")&(sPRUA$sALT_DP>0))
sPRUA$LargeH <- ((sPRUA$FullVar %in% "R156H")&(sPRUA$sALT_DP>0))
SWLQvect <- (sPRUA$SAMPLE[sPRUA$LargeQ])[(sPRUA$SAMPLE[sPRUA$LargeQ]) %in% (sPRUA$SAMPLE[sPRUA$LargeH])]
sPRUA$LargeInterQ <- (sPRUA$CombID %in% paste0(SWLQvect, ":", "R156Q"))&(sPRUA$sALT_DP>0)
sPRUA$LargeInterH <- (sPRUA$CombID %in% paste0(SWLQvect, ":", "R156H"))&(sPRUA$sALT_DP>0)
sPRUA <- sPRUA[(sPRUA$LargeInterH == FALSE),]
sPRUA <- sPRUA[!((sPRUA$LargeQ == TRUE)&(sPRUA$LargeInterQ == FALSE)),]

sPRUA <- sPRUA[(sPRUA$mALT_DP >= 3),]
sPRUA <- sPRUA[(sPRUA$mFULL_DP >= 20),]

print("The following have both R156Q and R156H")
print(sPRUA$CombID[duplicated(str_replace(sPRUA$CombID, "R156Q", "R156H"))])


sPRUA$mVAF <- sPRUA$mALT_DP/sPRUA$mFULL_DP*100


write.table(sPRUA, file = paste0(ForR.tmp, "_OutU2AF1.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

write.table(BadVect2, file = paste0(ForR.tmp, "_ZeroFiles.tsv"), sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)


