#!/usr/local/bin/Rscript
library(dplyr, quietly=T)
library(stringr, quietly=T)
library(data.table, quietly=T)
rm(list=ls())


gList<-read.csv("LCHIPPATH_GeneList_Bick_practical.csv")


vars1<-read.csv("PutativeOut.hg38_multianno.csv")

vars1$non_cancer_AF_popmax <- ifelse(vars1$non_cancer_AF_popmax == ".", 0, as.numeric(vars1$non_cancer_AF_popmax))

vars1 <- vars1[vars1$non_cancer_AF_popmax < 0.001,]
vars1 <- vars1[(vars1$ExonicFunc.refGene != "synonymous SNV"),]

vars1$BACKUP_Gene.refGene <- vars1$Gene.refGene

vars1$Gene.refGene <- ifelse(grepl("KIR2DL3", vars1$Gene.refGene), "KIR2DL3", vars1$Gene.refGene)
vars1$Gene.refGene <- ifelse(grepl("MEF2B", vars1$Gene.refGene), "MEF2B", vars1$Gene.refGene)
vars1$Gene.refGene <- ifelse(grepl("CXCL16|MED11", vars1$Gene.refGene), "MED11", vars1$Gene.refGene)
vars1$Gene.refGene <- ifelse(grepl("PIK3CA", vars1$Gene.refGene), "PIK3CA", vars1$Gene.refGene)
vars1$Gene.refGene <- ifelse(grepl("NTRK1", vars1$Gene.refGene), "NTRK1", vars1$Gene.refGene)
vars1$Gene.refGene <- ifelse(grepl("TSC2", vars1$Gene.refGene), "TSC2", vars1$Gene.refGene)
vars1$Gene.refGene <- ifelse(grepl("H1-4", vars1$Gene.refGene), "H1-4", vars1$Gene.refGene)

vars1$Gene.refGene <- str_remove_all(vars1$Gene.refGene, ";.*")

table(vars1$Gene.refGene, vars1$BACKUP_Gene.refGene)%>%
  data.frame()%>%
  filter(Freq > 0)%>%
  mutate(Var1 = as.character(Var1), Var2 = as.character(Var2))%>%
  filter(Var1 != Var2)

vars<-transform(vars1)



colnames=c("Chr","Start","End","Ref","Alt","Func.refGene","Gene.refGene","GeneDetail.refGene","ExonicFunc.refGene","AAChange.refGene")
all(names(vars)[1:10] == colnames)
rm(colnames)


vars$Gene.refGene[!vars$Gene.refGene%in%gList$Gene]
varsOI<-vars[vars$Gene.refGene%in%gList$Gene,]


varsOI.func<-varsOI[varsOI$Gene.refGene%in%gList$Gene&
                      (grepl("exonic",varsOI$Func.refGene,fixed=T)|
                         grepl("splicing",varsOI$Func.refGene,fixed=T)),]

varsOI.func<-merge(varsOI.func,gList[,c("Gene", "Accession")],by.x="Gene.refGene", by.y="Gene")


#Func.refGene # (exonic, splicing or some wierd merged combination)
#GeneDetail.refGene # NM_017940:exon16:c.1380-2A>G
#ExonicFunc.refGene #"nonsynonymous SNV"
#AAChange.refGene #"ASXL1:NM_015338:exon11:c.G1718A:p.R573Q"
#Accession          "NM_015338"
extractTranscript<-function(GeneDetail,AAChange,Accession){
  splice<-grep(Accession,strsplit(GeneDetail,",")[[1]], value=T, fixed=T)
  nonsyn<-grep(Accession,strsplit(AAChange,",")[[1]], value=T, fixed=T)
  if(length(splice)>0){return(splice)} else {
    if(length(nonsyn)>0){return(nonsyn)} else{
      return("nan")}
  }
}

#Apply function
aa<-apply(varsOI.func[,c('GeneDetail.refGene', 'AAChange.refGene','Accession')],
          1, function(x) {extractTranscript(x[1],x[2],x[3]) })
aa2<-lapply(aa, `[[`, 1)
varsOI.func$transcriptOI<-unlist(aa2)
varsOI.func<-varsOI.func[varsOI.func$transcriptOI!="nan",]
rm(aa,aa2)

extractNonsyn<-function(transcriptOI){
  protChange<-grep("p.",strsplit(transcriptOI,":")[[1]], value=T,fixed=T)
  if(length(protChange)>0){return(gsub("p\\.","",protChange))} else {return("nan")}
}


varsOI.func$NonsynOI <- NA


for(i in 1:length(varsOI.func$NonsynOI)){
  
  varsOI.func$NonsynOI[i] <- extractNonsyn(((varsOI.func[,c('transcriptOI')])[[i]]))  
  
}




# Columns to be removed if they exist
if ("aalen" %in% colnames(varsOI.func)) {
  varsOI.func <- varsOI.func[,-c("aalen", "aapos", "aafirst10pctPeptide", "aalast10pctPeptide", "LOFfirst10pct", "LOFlast10pct")]
}

varsOI.func$NonsynOI <- ifelse(varsOI.func$GeneDetail.refGene != "." & varsOI.func$NonsynOI == "nan", str_extract(varsOI.func$transcriptOI, paste0(varsOI.func$Accession, "[^;]*;"))%>%str_remove(";")%>%str_remove_all(".*:"),
       paste0("p.", varsOI.func$NonsynOI))



varsOI.func$JoiningVar <- paste(varsOI.func$Chr, varsOI.func$Start, varsOI.func$Ref, varsOI.func$Alt, sep = ":")
vars1$JoiningVar <- paste(vars1$Chr, vars1$Start, vars1$Ref, vars1$Alt, sep = ":")


vars1 <- vars1%>%
  mutate(Accession = NA,
         transcriptOI = NA,
         NonsynOI = NA,
         whitelist = NA,
         wl.mis = NA,
         wl.lof = NA,
         wl.splice = NA,
         wl.exception = NA,
         manualreview = NA)

vars1 <- vars1[!vars1$JoiningVar %in% varsOI.func$JoiningVar,]

ColVect <- colnames(varsOI.func)

vars1 <- vars1%>%
  select(ColVect)

varsOI.func <- rbind(varsOI.func, vars1)

fin_df <- read.delim("FirstFilt.tsv", sep = "\t")

fin_df <- fin_df[!grepl("multiallelic", fin_df$FILTER),]


fin_df$REF <- ifelse(fin_df$REF == "TRUE", "T", fin_df$REF)
fin_df$ALT <- ifelse(fin_df$ALT == "TRUE", "T", fin_df$ALT)


fin_df$X.CHROM <- str_remove(fin_df$X.CHROM, "chr")
fin_df$JoiningVar <- paste(fin_df$X.CHROM, fin_df$POS, fin_df$REF, fin_df$ALT, sep = ":")

fin_df <- left_join(fin_df, varsOI.func, by = c("JoiningVar"))
fin_df <- fin_df[!is.na(fin_df$Accession),]



# Removing duplicated records ---------------------------------------------

fin_df$JoiningVar <- paste(fin_df$SampID, fin_df$Chr, fin_df$Start, fin_df$End, fin_df$Ref, fin_df$Alt, sep = ":")

#fin_df <- fin_df[!duplicated(fin_df$JoiningVar),]


print(paste0("There are ", sum(duplicated(fin_df$JoiningVar)), " duplicates, remember duplicates are not automatically removed here."))

# Making the semi filtered files ------------------------------------------

write.table(fin_df, file = "Putative_SemiFilteredDF.tsv", quote = FALSE, row.names = FALSE, sep = "\t")


















