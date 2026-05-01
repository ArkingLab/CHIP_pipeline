## Convert MCHIP to long format and generate final YES/NO file for MCHIP and LCHIP

## packages
library(data.table)
library(tidyverse)
library(magrittr)
library(lubridate)
library(stringr)

## Read min MCHIP and LCHIP final filtered tsv files
mchip = fread('MCHIP_FinalFiltered.tsv')
lchip = fread('LCHIP_FinalFiltered.tsv')
## both need to fix SampID, as mchip some has INDEL, l chip all NA
mchip = mchip %>%
  mutate(
    SampID = str_extract(INFO, "(?<=SM=)[^;]+")
  )
lchip = lchip %>%
  mutate(
    SampID = str_extract(INFO, "(?<=SM=)[^;]+")
  )
# Read in sex information for filtering downstream, can use idsex.txt used in earlier CHIP step
sex = fread("idsex.txt")
# Read in ARIC and UKB germline list
aric_germline = fread('ARIC_germline.txt') %>% as.data.frame() %>% pull(x)
## For UKB, use the dp5 germline list with 261 variants
#write.table(germline_list_ukb_dp5, file='UKB_germline.txt', sep="\t", row.names=FALSE, quote = FALSE)
ukb_germline = fread('UKB_germline.txt') %>% as.data.frame() %>% pull(x)
# union genes
genes = fread("FullGeneList.csv") %>% as.data.frame()

######################## MCHIP process for diff dp ##############################

######## adjustment for MCHIP ################

#### adjust AF for men&chrX
mchip1 = mchip %>% 
  mutate(id = SampID, eid = SampID) %>%
  left_join(
    sex,
    by = "id") %>%
  mutate(VAF_adj = case_when(
    X.CHROM == "X" & sex == "men" ~ newVAF/2, 
    TRUE ~ newVAF
  ))

#### flag VAF_adj < 2 and Multinucleotide variants with VAF < 10 
mchip2 = 
  mchip1 |>
  mutate(
    # `VAF_adj` < 2
    filter_vaf = ifelse(VAF_adj < 2, 1, 0), 
    # multiallelic variants with VAF < 10 
    filter_mnv = ifelse(
      (nchar(REF) != 1 | nchar(ALT) != 1) & VAF_adj < 10, 1, 0
    ), 
    # temporary variable to apply for germline filtering  
    SubBinomGermline = 
      ifelse(
        (sex == "men" & X.CHROM == "X"),
        (BinomTestGermlineHom >= 0.01),
        (BinomTestGermlineHet >= 0.01) | (BinomTestGermlineHom >= 0.01)
      )
  )

#### filter out those above
mchip3 =
  mchip2 |>
  filter(
    filter_vaf == 0 & filter_mnv == 0
  )

#### summarize
mchip3 |> 
  summarise(
    n = n_distinct(id),
    obs = n()             
  )

############ Now apply different alt DP cutoffs, smallest is 3 already, need to apply 456 ##################
mchip_dp3 = mchip3  

mchip_dp4 = 
  mchip3  |>
  filter(AltDP >= 4)

mchip_dp5 = 
  mchip3  |>
  filter(AltDP >= 5)

mchip_dp6 = 
  mchip3  |>
  filter(AltDP >= 6)

########## Start with germline definition for each dp ##########################
for (dp in 3:6) {
  
  # get original dataset (chip_dp*)
  mchip <- get(paste0("mchip_dp", dp))
  
  # ---- main filtering (optional print) ----
  mchip %>%
    filter(SubBinomGermline != BinomGermline) %>%
    select(FullVar, transcriptOI, eid)
  
  # ---- frequency calculations ----
  freq_all <- mchip %>%
    group_by(FullVar) %>%
    summarise(n = n(), .groups = "drop")
  
  freq_sub <- mchip %>%
    filter(SubBinomGermline == TRUE) %>%
    group_by(FullVar) %>%
    summarise(n = n(), .groups = "drop")
  
  freq <- freq_all %>%
    left_join(freq_sub, by = "FullVar")
  
  print(freq %>% count(is.na(n.y)))
  
  freq_nonmiss <- freq %>%
    filter(!is.na(n.y)) %>%
    mutate(prop = n.y / n.x)
  
  print(summary(freq_nonmiss$prop))
  
  # ---- germline list ----
  germline_list <- freq_nonmiss %>%
    filter(n.x >= 3 & prop >= 0.8) %>%
    pull(FullVar)
  
  print(germline_list)
  
  # ---- annotate dataset ----
  mchip <- mchip %>%
    mutate(
      SubPotentialGermline = FullVar %in% germline_list
    )
  
  print(table(mchip$SubPotentialGermline))
  print(table(mchip$PotentialGermline,
              mchip$SubPotentialGermline,
              useNA = "ifany"))
  
  mchip <- mchip %>%
    mutate(
      filter_germline = ifelse(
        SubPotentialGermline == TRUE |
          transcriptOI %in% aric_germline |
          FullVar %in% ukb_germline,
        1, 0
      )
    )
  
  print(table(mchip$filter_germline, useNA = "ifany"))
  
  # ---- save outputs with dp-specific names ----
  assign(paste0("mchip_dp", dp), mchip)
  assign(paste0("freq_dp", dp), freq)
  assign(paste0("freq_nonmiss_dp", dp), freq_nonmiss)
  assign(paste0("germline_list_dp", dp), germline_list)
}
# filter out
for (dp in 3:6) {
  
  # get processed dataset (mchip_dp*)
  mchip <- get(paste0("mchip_dp", dp))
  
  # filter
  mchip_filtered <- mchip %>%
    filter(filter_germline == 0)
  
  # summary
  summary_stats <- mchip_filtered %>%
    summarise(
      n = n_distinct(eid),
      obs = n()
    )
  
  print(paste0("DP = ", dp))
  print(summary_stats)
  
  # assign outputs
  assign(paste0("mchip_filtered_dp", dp), mchip_filtered)
  assign(paste0("summary_dp", dp), summary_stats)
}
# union genes
gene_list_all = genes |> pull(Union_Genes) 
for (dp in 3:6) {
  
  # get filtered dataset
  mchip_filtered <- get(paste0("mchip_filtered_dp", dp))
  
  # subset to union gene list
  mchip_union <- mchip_filtered %>%
    filter(Gene.refGene %in% gene_list_all)
  
  # summary
  summary_stats <- mchip_union %>%
    summarise(
      n = n_distinct(eid),
      obs = n()
    )
  
  print(paste0("DP = ", dp))
  print(summary_stats)
  
  # assign outputs
  assign(paste0("mchip_union_dp", dp), mchip_union)
  assign(paste0("summary_union_dp", dp), summary_stats)
}
# indels, use cutoff of 4
for (dp in 3:6) {
  
  # get union dataset
  mchip_union <- get(paste0("mchip_union_dp", dp))
  
  # all variants per sample
  check_indel_all <- mchip_union %>%
    group_by(SampID) %>%
    summarise(n = n(), .groups = "drop") %>%
    arrange(desc(n)) %>%
    filter(n >= 4)
  
  # MNV/indel-like variants per sample
  check_indel_mnv <- mchip_union %>%
    filter(nchar(REF) > 1 | nchar(ALT) > 1) %>%
    group_by(SampID) %>%
    summarise(n = n(), .groups = "drop") %>%
    arrange(desc(n)) %>%
    filter(n >= 4)
  
  # compare all vs MNV/indel-like hits
  check_indel <- check_indel_all %>%
    left_join(check_indel_mnv, by = "SampID") %>%
    mutate(
      n.y = ifelse(is.na(n.y), 0, n.y)
    )
  
  # add high_indel flag
  mchip_union <- mchip_union %>%
    mutate(
      high_indel = ifelse(
        SampID %in% check_indel$SampID[check_indel$n.x == check_indel$n.y],
        1, 0
      )
    )
  
  print(paste0("DP = ", dp))
  print(sum(mchip_union$HighNumberOfIndels, na.rm = TRUE))
  print(sum(mchip_union$high_indel, na.rm = TRUE))
  
  # assign outputs
  assign(paste0("check_indel_all_dp", dp), check_indel_all)
  assign(paste0("check_indel_mnv_dp", dp), check_indel_mnv)
  assign(paste0("check_indel_dp", dp), check_indel)
  assign(paste0("mchip_union_dp", dp), mchip_union)
}
# new variables, final section
for (dp in 3:6) {
  
  chip_combined <- get(paste0("mchip_union_dp", dp))
  
  chip_all <- chip_combined %>%
    mutate(
      VAF_raw = newVAF,
      VAF = VAF_adj
    ) %>%
    select(-newVAF) %>%
    group_by(eid) %>%
    mutate(
      maxVAF_ind = max(VAF_adj, na.rm = TRUE),
      chip_count = n()
    ) %>%
    ungroup() %>%
    group_by(eid, Gene.refGene) %>%
    mutate(
      maxVAF_chip = max(VAF_adj, na.rm = TRUE)
    ) %>%
    ungroup() %>%
    arrange(eid)
  
  chip_count_dist <- chip_all %>%
    group_by(chip_count) %>%
    summarise(n = n(), .groups = "drop")
  
  gene_count_dist <- chip_all %>%
    group_by(Gene.refGene) %>%
    summarise(
      n_id = n_distinct(eid),
      n_obs = n(),
      .groups = "drop"
    ) %>%
    arrange(desc(n_id))
  
  print(paste0("DP = ", dp))
  print(chip_count_dist)
  print(gene_count_dist)
  
  assign(paste0("mchip_all_dp", dp), chip_all)
  assign(paste0("chip_count_dist_dp", dp), chip_count_dist)
  assign(paste0("gene_count_dist_dp", dp), gene_count_dist)
}

####### After all run, save into diff dp MCHIP files
write.table(mchip_all_dp3, file='MCHIP_dp3.txt', sep="\t", row.names=FALSE, quote = FALSE)
write.table(mchip_all_dp4, file='MCHIP_dp4.txt', sep="\t", row.names=FALSE, quote = FALSE)
write.table(mchip_all_dp5, file='MCHIP_dp5.txt', sep="\t", row.names=FALSE, quote = FALSE)
write.table(mchip_all_dp6, file='MCHIP_dp6.txt', sep="\t", row.names=FALSE, quote = FALSE)



######################## M and L CHIP yes no ###################################

## Read in list of full samples
sample = fread('SPARK.iWGS_v1.1.mastertable.2023_03.tsv')
## select variables to keep
sample = sample %>% select(spid,sfid)

## Now mutate variables into whether each sample has MCHIP or LCHIP
sample1 = sample %>%
  mutate(MCHIP_dp3 = ifelse(spid %in% mchip_all_dp3$eid,1,0),
         MCHIP_dp4 = ifelse(spid %in% mchip_all_dp4$eid,1,0),
         MCHIP_dp5 = ifelse(spid %in% mchip_all_dp5$eid,1,0),
         MCHIP_dp6 = ifelse(spid %in% mchip_all_dp6$eid,1,0),
         LCHIP = ifelse(spid %in% lchip$SampID,1,0)
         )

# qc
table(sample1$MCHIP_dp3,useNA = 'always')
table(sample1$MCHIP_dp4,useNA = 'always')
table(sample1$MCHIP_dp5,useNA = 'always')
table(sample1$MCHIP_dp6,useNA = 'always')
table(sample1$LCHIP,useNA = 'always')

## Save the Yes No file
write.table(sample1, file='M_L_CHIP_yesno_050126.txt', sep="\t", row.names=FALSE, quote = FALSE)


