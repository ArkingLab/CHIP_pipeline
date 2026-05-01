#!/usr/local/bin/Rscript
library(dplyr, quietly = TRUE)
library(stringr, quietly = TRUE)
library(data.table, quietly=T)
rm(list=ls())

fin_df <- read.delim("Input.vcf", skip = 1)#skip the VCF lines until the X.CHROM line

colnames(fin_df) <- c("X.CHROM",  "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE")

print(paste0("Unfiltered number of lines is ", length(fin_df$X.CHROM)))

# Back to generating the ANNOVAR input ------------------------------------

fin_df$REF <- ifelse(fin_df$REF == "TRUE", "T", fin_df$REF)
fin_df$ALT <- ifelse(fin_df$ALT == "TRUE", "T", fin_df$ALT)

# Mutect2 FILTER ----------------------------------------------------------

Mutect2All <- ((unique(fin_df$FILTER)%>%
                  paste(., collapse = ";")%>%
                  str_split(., ";"))[[1]])%>%
  unique()

Mutect2_OUT <- Mutect2All[!Mutect2All %in% c("PASS", "weak_evidence", "germline", ".")]

Mutect2_OUT <- paste(Mutect2_OUT, collapse = "|")

if(Mutect2_OUT == ""){
}else{
  fin_df <- fin_df[!grepl(Mutect2_OUT, fin_df$FILTER),]
}

print(paste0("After mutect2 filter number of lines is ", length(fin_df$X.CHROM)))

# Getting the SAMPLE field ------------------------------------------------
unique(fin_df$FORMAT)

fin_df <- fin_df[!grepl(",", str_remove(fin_df$SAMPLE, "^[^:]*:[^:]*:")%>%
                          str_remove(., ":.*")%>%
                          as.character()),]


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

write.table(fin_df, file = "FirstFilt.tsv", row.names = FALSE, sep = "\t")

fin_df$End <- fin_df$POS+(nchar(fin_df$REF)-1)

fin_df$X.CHROM <- str_remove(fin_df$X.CHROM, "chr")

fin_df <- fin_df%>%
  select(X.CHROM, POS, End, REF, ALT)%>%
  distinct(X.CHROM, POS, End, REF, ALT)

write.table(fin_df, file = "fin_df_Genome.input.txt", quote = FALSE, row.names = FALSE,
            col.names = FALSE, sep = "\t")
















