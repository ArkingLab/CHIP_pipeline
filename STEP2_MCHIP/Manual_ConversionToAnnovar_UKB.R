#!/usr/local/bin/Rscript
library(dplyr, quietly = TRUE)
library(stringr, quietly = TRUE)
library(data.table, quietly=T)
rm(list=ls())

fin_df <- read.delim("Input.vcf", skip = 0)

print(paste0("Unfiltered number of lines is ", length(fin_df$X.CHROM)))


# Look at multiallelic ----------------------------------------------------

MultDF <- fin_df[grepl("multiallelic", fin_df$FILTER),]

MultDF$X.CHROM <- str_remove(MultDF$X.CHROM, "chr")

MultDF$CombID <- paste(MultDF$X.CHROM, MultDF$POS, MultDF$REF, MultDF$ALT, sep = ":")

Mutect2All <- ((unique(MultDF$FILTER)%>%
                  paste(., collapse = ";")%>%
                  str_split(., ";"))[[1]])%>%
  unique()

Mutect2_OUT <- Mutect2All[!Mutect2All %in% c("PASS", "weak_evidence", "germline", "multiallelic")]

Mutect2_OUT <- paste(Mutect2_OUT, collapse = "|")

if(Mutect2_OUT == ""){
}else{
  MultDF <- MultDF[!grepl(Mutect2_OUT, MultDF$FILTER),]
}

MultDF$AltDP <- str_remove(MultDF$SAMPLE, "^[^:]*:")%>%
  str_remove(., ":.*")

MultDF <- MultDF[grepl("^0,|,0,|,0$", MultDF$AltDP),]

write.csv(MultDF, file = "WITHREF_PotentialVarsMultiallelic.csv", row.names = FALSE)

MultDF <- MultDF[grepl(",0,|,0$", MultDF$AltDP),]

write.csv(MultDF, file = "WITHoutREF_PotentialVarsMultiallelic.csv", row.names = FALSE)

# Back to generating the ANNOVAR input ------------------------------------

fin_df <- fin_df[!grepl("multiallelic", fin_df$FILTER),]

print(paste0("Without multiallelic number of lines is ", length(fin_df$X.CHROM)))

fin_df$REF <- ifelse(fin_df$REF == "TRUE", "T", fin_df$REF)
fin_df$ALT <- ifelse(fin_df$ALT == "TRUE", "T", fin_df$ALT)

# Mutect2 FILTER ----------------------------------------------------------

Mutect2All <- ((unique(fin_df$FILTER)%>%
                  paste(., collapse = ";")%>%
                  str_split(., ";"))[[1]])%>%
  unique()

Mutect2_OUT <- Mutect2All[!Mutect2All %in% c("PASS", "weak_evidence", "germline")]

Mutect2_OUT <- paste(Mutect2_OUT, collapse = "|")

if(Mutect2_OUT == ""){
}else{
  fin_df <- fin_df[!grepl(Mutect2_OUT, fin_df$FILTER),]
}

print(paste0("After mutect2 filter number of lines is ", length(fin_df$X.CHROM)))

# Getting the SAMPLE field ------------------------------------------------
unique(fin_df$FORMAT)

fin_df$AltDP <- str_remove(fin_df$SAMPLE, "^[^:]*:")%>%
  str_remove(., ":.*")%>%
  str_remove(., ".*,")%>%
  as.character()%>%
  as.numeric()

fin_df$AltF1R2 <- str_remove(fin_df$SAMPLE, "^[^:]*:[^:]*:[^:]*:[^:]*:")%>%
  str_remove(., ":.*")%>%
  str_remove(., ".*,")%>%
  as.character()%>%
  as.numeric()

fin_df$AltF2R1 <- str_remove(fin_df$SAMPLE, "^[^:]*:[^:]*:[^:]*:[^:]*:[^:]*:")%>%
  str_remove(., ":.*")%>%
  str_remove(., ".*,")%>%
  as.character()%>%
  as.numeric()

fin_df$DP <- str_remove(fin_df$SAMPLE, "^[^:]*:[^:]*:[^:]*:")%>%
  str_remove(., ":.*")%>%
  as.character()%>%
  as.numeric()


#Getting VAF
fin_df$newVAF <- str_remove(fin_df$SAMPLE, "^[^:]*:[^:]*:")%>%
  str_remove(., ":.*")%>%
  as.character()%>%
  as.numeric()


fin_df$newVAF <- fin_df$newVAF*100


# Sequencing depth-based filtering ----------------------------------------

fin_df <- fin_df[!fin_df$DP < 20,]

print(paste0("After DP under 20 number of lines is ", length(fin_df$X.CHROM)))

fin_df <- fin_df[!(fin_df$AltDP < 3),]
print(paste0("After AltDP under 3 number of lines is ", length(fin_df$X.CHROM)))

fin_df <- fin_df[!((nchar(fin_df$REF)>1|nchar(fin_df$ALT)>1) & (fin_df$AltDP < 5)),]
print(paste0("After MNV Alt DP under 5 number of lines is ", length(fin_df$X.CHROM)))

fin_df <- fin_df[!(fin_df$AltF1R2 < 1),]
fin_df <- fin_df[!(fin_df$AltF2R1 < 1),]

print(paste0("After F1R2 and F2R1 under 1 number of lines is ", length(fin_df$X.CHROM)))

fin_df <- fin_df[!(fin_df$newVAF < 2),]

print(paste0("After VAF under 2 number of lines is ", length(fin_df$X.CHROM)))

fin_df <- fin_df[!((nchar(fin_df$REF)>1|nchar(fin_df$ALT)>1) & (fin_df$newVAF < 10)),]

print(paste0("After MNV VAF under 10 number of lines is ", length(fin_df$X.CHROM)))

fin_df$FullVar <- paste(fin_df$X.CHROM, fin_df$POS, fin_df$REF, fin_df$ALT, sep = ":")

si_df <- table(fin_df$FullVar)%>%
  data.frame()%>%
  arrange(-Freq)

write.table(fin_df, file = "FirstFilt.tsv", row.names = FALSE, sep = "\t")



fin_df$End <- fin_df$POS+(nchar(fin_df$REF)-1)

fin_df$X.CHROM <- str_remove(fin_df$X.CHROM, "chr")

fin_df <- fin_df%>%
  select(X.CHROM, POS, End, REF, ALT)%>%
  distinct(X.CHROM, POS, End, REF, ALT)

write.table(fin_df, file = "fin_df_Genome.input.txt", quote = FALSE, row.names = FALSE,
            col.names = FALSE, sep = "\t")
















