#!/usr/local/bin/Rscript
library(dplyr, quietly=T)
library(stringr, quietly=T)
library(data.table, quietly=T)
rm(list=ls())


gList<-read.csv("GeneList_Bick_practical.csv")

whitelist.mis<-read.csv("CHIP_missense_vars_Bick_practical.csv")

whitelist.splice<-read.csv("CHIP_splice_vars_Bick_practical.csv")

whitelist.LoF<-read.csv("CHIP_nonsense_FS_vars_Bick_practical.csv")



vars1<-read.csv("Out.hg38_multianno.csv")

vars1$non_cancer_AF_popmax <- ifelse(vars1$non_cancer_AF_popmax == ".", 0, as.numeric(vars1$non_cancer_AF_popmax))

vars1 <- vars1[vars1$non_cancer_AF_popmax < 0.001,]
vars1 <- vars1[(vars1$ExonicFunc.refGene != "synonymous SNV"),]


vars1$Gene.refGene <- str_remove_all(vars1$Gene.refGene, ";.*")
vars<-transform(vars1)

colnames=c("Chr","Start","End","Ref","Alt","Func.refGene","Gene.refGene","GeneDetail.refGene","ExonicFunc.refGene","AAChange.refGene")
all(names(vars)[1:10] == colnames)
rm(colnames)



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



#Annotate and filter
varsOI.func<-transform(varsOI.func,whitelist=F,
                       wl.mis=F,wl.lof=F,wl.splice=F,wl.exception=F,
                       manualreview=F)
#1) Handle missense vars
vmis<-varsOI.func$ExonicFunc.refGene=="nonsynonymous SNV"
vmis_wl<-paste(varsOI.func$Gene.refGene,varsOI.func$NonsynOI,sep="_")%in% paste(whitelist.mis$Gene,whitelist.mis$AAChange,sep="_")
varsOI.func[vmis&vmis_wl,"whitelist"]=T
varsOI.func[vmis&vmis_wl,"wl.mis"]=T

#2) Handle LoF and frame shift vars
vlof<-grepl("X",varsOI.func$NonsynOI, fixed=T) #stop gain or stop loss
vFS<-grepl("fs",varsOI.func$NonsynOI, fixed=T) #frameshift
vLOFgene<-varsOI.func$Gene.refGene%in%whitelist.LoF$Gene #genes are Lof Genes
varsOI.func[(vlof|vFS)&vLOFgene,"whitelist"]=T
varsOI.func[(vlof|vFS)&vLOFgene,"wl.lof"]=T

#3) Handle Splice vars
#possible splicing variant
vSplice<-grepl("splicing",varsOI.func$Func.refGene, fixed=T)
#genes are splice Genes
vSplicegene<-varsOI.func$Gene.refGene%in%whitelist.splice$Gene
#confirm that the splicing refers to the correct transcript
vSpliceCorrectTranscript<-apply(varsOI.func[,c('GeneDetail.refGene','Accession')],
                                1, function(x) {grepl(x[2],x[1],fixed=T)})
varsOI.func[vSplice&vSplicegene&vSpliceCorrectTranscript,"whitelist"]=T
varsOI.func[vSplicegene&vSplicegene&vSpliceCorrectTranscript,"wl.splice"]=T
#If not the correct transcript, flag for manual review
varsOI.func[(vSplice&vSplicegene) & (!vSpliceCorrectTranscript),"manualreview"]=T

#4) Handle the following exceptions
#ASXL1	Frameshift/nonsense/splice-site in exon 11-12
#LoF
vlof<-grepl("X",varsOI.func$NonsynOI, fixed=T) #stop gain or stop loss
vFS<-grepl("fs",varsOI.func$NonsynOI, fixed=T) #frameshift
vexon11<-grepl("exon11",varsOI.func$transcriptOI, fixed=T)
vexon12<-grepl("exon12",varsOI.func$transcriptOI, fixed=T)
asxl1Exception<-(varsOI.func$Gene.refGene=="ASXL1")&(vlof|vFS)&(vexon11|vexon12)
varsOI.func[asxl1Exception,"whitelist"]=T
varsOI.func[asxl1Exception,"wl.lof"]=T
varsOI.func[asxl1Exception,"wl.exception"]=T
#Splicing
vSplice<-grepl("splicing",varsOI.func$Func.refGene, fixed=T)
vSpliceCorrectTranscript<-apply(varsOI.func[,c('GeneDetail.refGene','Accession')],
                                1, function(x) {grepl(x[2],x[1],fixed=T)})
vexon11splice<-grepl("exon11",varsOI.func$GeneDetail.refGene, fixed=T)
vexon12splice<-grepl("exon12",varsOI.func$GeneDetail.refGene, fixed=T)
asxl1ExceptionSplice<-(varsOI.func$Gene.refGene=="ASXL1")&
  vSplice&vSpliceCorrectTranscript&
  (vexon11splice|vexon12splice)
varsOI.func[asxl1ExceptionSplice,"whitelist"]=T
varsOI.func[asxl1ExceptionSplice,"wl.splice"]=T
varsOI.func[asxl1ExceptionSplice,"wl.exception"]=T



#CALR	Frameshift/nonsense/splice-site in exon 9
#LoF
vlof<-grepl("X",varsOI.func$NonsynOI, fixed=T) #stop gain or stop loss
vFS<-grepl("fs",varsOI.func$NonsynOI, fixed=T) #frameshift
vexon9<-grepl("exon9",varsOI.func$transcriptOI, fixed=T)
CALRException<-(varsOI.func$Gene.refGene=="CALR")&(vlof|vFS)&(vexon9)
varsOI.func[CALRException,"whitelist"]=T
varsOI.func[CALRException,"wl.lof"]=T
varsOI.func[CALRException,"wl.exception"]=T
#Splicing
vSplice<-grepl("splicing",varsOI.func$Func.refGene, fixed=T)
vSpliceCorrectTranscript<-apply(varsOI.func[,c('GeneDetail.refGene','Accession')],
                                1, function(x) {grepl(x[2],x[1],fixed=T)})
vexon9splice<-grepl("exon9",varsOI.func$GeneDetail.refGene, fixed=T)
CALRExceptionSplice<-(varsOI.func$Gene.refGene=="CALR")&
  vSplice&vSpliceCorrectTranscript&
  (vexon9splice)
varsOI.func[CALRExceptionSplice,"whitelist"]=T
varsOI.func[CALRExceptionSplice,"wl.splice"]=T
varsOI.func[CALRExceptionSplice,"wl.exception"]=T


#ASXL2	Frameshift/nonsense/splice-site in exon 11-12
#Lof
asxl2Exception<-(varsOI.func$Gene.refGene=="ASXL2")&(vlof|vFS)&(vexon11|vexon12)
varsOI.func[asxl2Exception,"whitelist"]=T
varsOI.func[asxl2Exception,"wl.lof"]=T
varsOI.func[asxl2Exception,"wl.exception"]=T
#Splice
asxl2ExceptionSplice<-(varsOI.func$Gene.refGene=="ASXL2")&
  vSplice&vSpliceCorrectTranscript&
  (vexon11splice|vexon12splice)
varsOI.func[asxl2ExceptionSplice,"whitelist"]=T
varsOI.func[asxl2ExceptionSplice,"wl.splice"]=T
varsOI.func[asxl2ExceptionSplice,"wl.exception"]=T


#PPM1D	Frameshift/nonsense in exon 5 or 6
vlof<-grepl("X",varsOI.func$NonsynOI, fixed=T) #stop gain or stop loss
vFS<-grepl("fs",varsOI.func$NonsynOI, fixed=T) #frameshift
vexon5<-grepl("exon5",varsOI.func$transcriptOI, fixed=T)
vexon6<-grepl("exon6",varsOI.func$transcriptOI, fixed=T)
ppm1dException<-(varsOI.func$Gene.refGene=="PPM1D")&(vlof|vFS)&(vexon5|vexon6)
varsOI.func[ppm1dException,"whitelist"]=T
varsOI.func[ppm1dException,"wl.lof"]=T
varsOI.func[ppm1dException,"wl.exception"]=T

#TET2	missense mutations in catalytic domains (p.1104-1481 and 1843-2002)
TETidx<-which(varsOI.func$Gene.refGene=="TET2"&
                varsOI.func$ExonicFunc.refGene=="nonsynonymous SNV"&
                nchar(varsOI.func$NonsynOI)==6)
for(i in TETidx){
  AApos<-as.numeric(substr(varsOI.func$NonsynOI[i],2,5))
  if((AApos>=1104&AApos<=1481)|(AApos>=1843&AApos<=2002))
  {
    varsOI.func[i,"whitelist"]=T
    varsOI.func[i,"wl.mis"]=T
    varsOI.func[i,"wl.exception"]=T
  }
}


  

#TET2	missense mutations in catalytic domains (p.1104-1481 and 1843-2002)
# TETidx<-which(varsOI.func$Gene.refGene=="TET2"&
#                 varsOI.func$ExonicFunc.refGene=="nonsynonymous SNV")
# for(i in TETidx){
#   AApos<-str_remove((varsOI.func$NonsynOI[i]), "[A-Z]$")%>%
#     str_remove(., ".*[A-Z]")%>%
#     as.character()%>%
#     as.numeric()
#   if((AApos>=1104&AApos<=1481)|(AApos>=1843&AApos<=2002))
#   {
#     varsOI.func[i,"whitelist"]=T
#     varsOI.func[i,"wl.mis"]=T
#     varsOI.func[i,"wl.exception"]=T
#   }
# }


#ZBTB33	first missense mutations in functional domains
ZBTBfirstidx<-which(varsOI.func$Gene.refGene=="ZBTB33"&
                 varsOI.func$ExonicFunc.refGene=="nonsynonymous SNV"&
               nchar(varsOI.func$NonsynOI)==3)
for(i in ZBTBfirstidx){
  AApos<-as.numeric(substr(varsOI.func$NonsynOI[i],2,2))
  if((AApos>=1&AApos<=120)|(AApos>=330&AApos<=600))
  {
    varsOI.func[i,"whitelist"]=T
    varsOI.func[i,"wl.mis"]=T
    varsOI.func[i,"wl.exception"]=T
  }
}

#ZBTB33	second missense mutations in functional domains
ZBTBsecondidx<-which(varsOI.func$Gene.refGene=="ZBTB33"&
                      varsOI.func$ExonicFunc.refGene=="nonsynonymous SNV"&
                    nchar(varsOI.func$NonsynOI)==4)
for(i in ZBTBsecondidx){
  AApos<-as.numeric(substr(varsOI.func$NonsynOI[i],2,3))
  if((AApos>=1&AApos<=120)|(AApos>=330&AApos<=600))
  {
    varsOI.func[i,"whitelist"]=T
    varsOI.func[i,"wl.mis"]=T
    varsOI.func[i,"wl.exception"]=T
  }
}

#ZBTB33	third missense mutations in functional domains
ZBTBthirdidx<-which(varsOI.func$Gene.refGene=="ZBTB33"&
                      varsOI.func$ExonicFunc.refGene=="nonsynonymous SNV"&
                    nchar(varsOI.func$NonsynOI)==5)
for(i in ZBTBthirdidx){
  AApos<-as.numeric(substr(varsOI.func$NonsynOI[i],2,4))
  if((AApos>=1&AApos<=120)|(AApos>=330&AApos<=600))
  {
    varsOI.func[i,"whitelist"]=T
    varsOI.func[i,"wl.mis"]=T
    varsOI.func[i,"wl.exception"]=T
  }
}


#ZBTB33	missense mutations in functional domains
# ZBTBidx<-which(varsOI.func$Gene.refGene=="ZBTB33"&
#                 varsOI.func$ExonicFunc.refGene=="nonsynonymous SNV")
# for(i in ZBTBidx){
#   AApos<-str_remove((varsOI.func$NonsynOI[i]), "[A-Z]$")%>%
#     str_remove(., ".*[A-Z]")%>%
#     as.character()%>%
#     as.numeric()
#   if((AApos>=1&AApos<=150)|(AApos>=330))
#   {
#     varsOI.func[i,"whitelist"]=T
#     varsOI.func[i,"wl.mis"]=T
#     varsOI.func[i,"wl.exception"]=T
#   }
# }

#CBL	RING finger missense p.381-421
CBLidx<-which(varsOI.func$Gene.refGene=="CBL"&
                varsOI.func$ExonicFunc.refGene=="nonsynonymous SNV"&
                nchar(varsOI.func$NonsynOI)==5)
for(i in CBLidx){
  AApos<-as.numeric(substr(varsOI.func$NonsynOI[i],2,4))
  if(AApos>=381&AApos<=421)
  {
    varsOI.func[i,"whitelist"]=T
    varsOI.func[i,"wl.mis"]=T
    varsOI.func[i,"wl.exception"]=T
  }
}

#CBL	RING finger missense p.381-421
# CBLidx<-which(varsOI.func$Gene.refGene=="CBL"&
#                 varsOI.func$ExonicFunc.refGene=="nonsynonymous SNV")
# for(i in CBLidx){
#   AApos<-str_remove((varsOI.func$NonsynOI[i]), "[A-Z]$")%>%
#     str_remove(., ".*[A-Z]")%>%
#     as.character()%>%
#     as.numeric()
#   if(AApos>=381&AApos<=421)
#   {
#     varsOI.func[i,"whitelist"]=T
#     varsOI.func[i,"wl.mis"]=T
#     varsOI.func[i,"wl.exception"]=T
#   }
# }

#CBLB	RING finger missense p.372-412
CBLBidx<-which(varsOI.func$Gene.refGene=="CBLB"&
                 varsOI.func$ExonicFunc.refGene=="nonsynonymous SNV"&
                 nchar(varsOI.func$NonsynOI)==5)
for(i in CBLBidx){
  AApos<-as.numeric(substr(varsOI.func$NonsynOI[i],2,4))
  if(AApos>=372&AApos<=412)
  {
    varsOI.func[i,"whitelist"]=T
    varsOI.func[i,"wl.mis"]=T
    varsOI.func[i,"wl.exception"]=T
  }
}


#CBLB	RING finger missense p.372-412
# CBLBidx<-which(varsOI.func$Gene.refGene=="CBLB"&
#                  varsOI.func$ExonicFunc.refGene=="nonsynonymous SNV")
# for(i in CBLBidx){
#   AApos<-str_remove((varsOI.func$NonsynOI[i]), "[A-Z]$")%>%
#     str_remove(., ".*[A-Z]")%>%
#     as.character()%>%
#     as.numeric()
#   if(AApos>=372&AApos<=412)
#   {
#     varsOI.func[i,"whitelist"]=T
#     varsOI.func[i,"wl.mis"]=T
#     varsOI.func[i,"wl.exception"]=T
#   }
# }




#5) flag remaining exceptions for manual review

vlof<-grepl("X",varsOI.func$NonsynOI, fixed=T) #stop gain or stop loss
vFS<-grepl("fs",varsOI.func$NonsynOI, fixed=T) #frameshift
vSplice<-grepl("splicing",varsOI.func$Func.refGene)

#GATA3	Frameshift/nonsense/splice-site ZNF domain
varsOI.func[(vlof|vFS|vSplice)&(varsOI.func$Gene.refGene=="GATA3"),"manualreview"]=T

#ZBTB33	INDEL 519_520del but could also be others
varsOI.func[varsOI.func$ExonicFunc.refGene=="nonframeshift substitution" &
              (varsOI.func$Gene.refGene=="ZBTB33"),"manualreview"]=T

#CREBBP	S1680del
varsOI.func[varsOI.func$ExonicFunc.refGene=="nonframeshift substitution" &
              (varsOI.func$Gene.refGene=="CREBBP"),"manualreview"]=T
#CSF3R	truncating c.741-791
varsOI.func[(vlof|vFS|vSplice)&(varsOI.func$Gene.refGene=="CSF3R"),"manualreview"]=T
#DNMT3A	F732del,	F752del
varsOI.func[varsOI.func$ExonicFunc.refGene=="nonframeshift substitution" &
              (varsOI.func$Gene.refGene=="DNMT3A"),"manualreview"]=T
#EP300	VF1148_1149del
varsOI.func[varsOI.func$ExonicFunc.refGene=="nonframeshift substitution" &
              (varsOI.func$Gene.refGene=="EP300"),"manualreview"]=T
#FLT3	FY590-591GD,	del835
varsOI.func[varsOI.func$ExonicFunc.refGene=="nonframeshift substitution"
            &(varsOI.func$Gene.refGene=="FLT3"),"manualreview"]=T
#JAK2	del/ins537-539L, del/ins538-539L, del/ins540-543MK, del/ins540-544MK, del/ins541-543K, del542-543, del543-544,	ins11546-547
varsOI.func[varsOI.func$ExonicFunc.refGene=="nonframeshift substitution"
            &(varsOI.func$Gene.refGene=="JAK2"),"manualreview"]=T
#KDM6A	del419
varsOI.func[(varsOI.func$ExonicFunc.refGene=="nonframeshift substitution")
            &(varsOI.func$Gene.refGene=="KDM6A"),"manualreview"]=T
#KIT	ins503,	del560,	del579,	del551-559
varsOI.func[varsOI.func$ExonicFunc.refGene=="nonframeshift substitution"
            &(varsOI.func$Gene.refGene=="KIT"),"manualreview"]=T
#MPL	del513 W515-518KT
varsOI.func[varsOI.func$ExonicFunc.refGene=="nonframeshift substitution"
            &(varsOI.func$Gene.refGene=="MPL"),"manualreview"]=T
#NPM1	Frameshift p.W288fs (insertion at c.859_860, 860_861, 862_863, 863_864)
varsOI.func[varsOI.func$ExonicFunc.refGene=="frameshift substitution"
            &(varsOI.func$Gene.refGene=="NPM1"),"manualreview"]=T

# Removed length of protein code 

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
print(paste0("VCF main has ", length(rownames(fin_df)), " lines"))
fin_df <- fin_df[!grepl("multiallelic", fin_df$FILTER),]
print(paste0("VCF main without multiallelic has ", length(rownames(fin_df)), " lines"))

fin_df$REF <- ifelse(fin_df$REF == "TRUE", "T", fin_df$REF)
fin_df$ALT <- ifelse(fin_df$ALT == "TRUE", "T", fin_df$ALT)


fin_df$X.CHROM <- str_remove(fin_df$X.CHROM, "chr")
fin_df$JoiningVar <- paste(fin_df$X.CHROM, fin_df$POS, fin_df$REF, fin_df$ALT, sep = ":")


fin_df <- left_join(fin_df, varsOI.func, by = c("JoiningVar"))
print(paste0("VCF main after left joining ", length(rownames(fin_df)), " lines"))
fin_df <- fin_df[!is.na(fin_df$whitelist),]
# Removing duplicated records ---------------------------------------------

fin_df$SampID <- str_extract(fin_df$INFO, "SM=.*")%>%
  str_remove(., "SM=")%>%
  str_remove_all(., "_.*")%>%
  as.character()

fin_df$JoiningVar <- paste(fin_df$SampID, fin_df$Chr, fin_df$Start, fin_df$End, fin_df$Ref, fin_df$Alt, sep = ":")
si_df <- fin_df[fin_df$JoiningVar %in% fin_df$JoiningVar[duplicated(fin_df$JoiningVar)],]

#fin_df <- fin_df[!duplicated(fin_df$JoiningVar),]


print(paste0("There are ", sum(duplicated(fin_df$JoiningVar)), " duplicates, remember duplicates are not automatically removed here."))

# Making the semi filtered files ------------------------------------------

ManCheck <- fin_df[fin_df$manualreview,]%>%
  distinct(Gene.refGene, transcriptOI)


write.csv(ManCheck, file = "ManualCheck.csv", row.names = FALSE)


fin_df$DiscordantLof <- ((fin_df$wl.lof) == (fin_df$ExonicFunc.refGene %in% c("stopgain", "stoploss", "frameshift substitution", "frameshift insertion", "frameshift deletion")))

fin_df$DiscordantLof <- !(fin_df$DiscordantLof)
fin_df$DiscordantLof <- ((fin_df$Gene.refGene %in% c(whitelist.LoF$Gene, "ASXL1", "ASXL2", "NPM1", "PPM1D", "CALR", "GATA3", "CSF3R", "SRSF2", "BRAF")) & fin_df$DiscordantLof)

DiscordantLofCheck <- fin_df[(fin_df$DiscordantLof),]%>%
  distinct(Gene.refGene, transcriptOI)

DiscordantLofCheck <- DiscordantLofCheck[!DiscordantLofCheck$transcriptOI %in% ManCheck$transcriptOI,]

write.csv(DiscordantLofCheck, file = "DiscordantLofCheck.csv", row.names = FALSE)


  
NewGenesCheck_DF <- (fin_df[fin_df$Gene.refGene %in% c("ATM","CALR","MYD88","PIGA","STAT3","NXF1","SRCAP","YLPM1","ZNF318","ZBTB33"),])
  
  
  
#write.csv(NewGenesCheck_DF, file = "NewGenesCheck_DF.csv", row.names = FALSE)
write.table(NewGenesCheck_DF, file = "NewGenesCheck_DF.tsv", quote = FALSE, row.names = FALSE, sep = "\t")
  
#write.csv(fin_df, file = "SemiFilteredDF.csv", row.names = FALSE)
write.table(fin_df, file = "SemiFilteredDF.tsv", quote = FALSE, row.names = FALSE, sep = "\t")


















