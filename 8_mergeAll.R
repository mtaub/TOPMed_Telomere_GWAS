
library(data.table)
library(dplyr)

### will want to get full set of all positions that cross 5e-8 in any population
## eventually will filter to 5e-9 for final table

baseLocation<-"/dcl01/mathias1/data/telomere_mtaub/gwas/results/"

allFileBase<-"ALL_GWASResults/allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop"
whiteFileBase<-"White_GWASResults/White_allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop"
blackFileBase<-"Black_GWASResults/Black_allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop"
asianFileBase<-"Asian_GWASResults/Asian_allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop"
HLFileBase<-"HispanicLatino_GWASResults/HL_allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop"
brazilianFileBase<-"Brazilian_GWASResults/Brazilian_allChrs_telomere_adjagesexbatchPCs_minDP0_BRAVODepthDrop"
samoanFileBase<-"Samoan_GWASResults/Samoan_allChrs_telomere_adjagesexseqctrbatchPCs_minDP0_BRAVODepthDrop"

allFileBases<-c(allFileBase, whiteFileBase, blackFileBase, asianFileBase, HLFileBase, brazilianFileBase, samoanFileBase)

########################################
## all results that have p < 5e-8 in any primary analysis
########################################

allSigIDs <- NULL
for (currBase in allFileBases){
  currRes<-read.csv(gzfile(paste0(baseLocation, currBase, "_p_lt_5e-8.csv.gz")), header = TRUE, stringsAsFactors = FALSE)
  allSigIDs <- union(allSigIDs, currRes$snpID)
}

allRes<-fread(paste0("gunzip -c ", baseLocation, allFileBase, ".csv.gz"), header=TRUE, sep=",")
allResSub <- allRes %>% filter(snpID %in% allSigIDs) %>% select(snpID, chr, pos, ref, alt, n.obs.ALL=n.obs, freq.ALL=freq, MAC.ALL=MAC, Score.ALL = Score, Score.SE.ALL=Score.SE,Score.Stat.ALL = Score.Stat,Score.pval.ALL = Score.pval,Est.ALL=Est,Est.SE.ALL=Est.SE,PVE.ALL=PVE, filt.ALL=filt, low_depth_rate.ALL=low_depth_rate)

for (currBase in allFileBases[-1]){
  currRes <- fread(paste0("gunzip -c ", baseLocation, currBase, ".csv.gz"), header=TRUE, sep=",")
  currResSub <- currRes %>% filter(snpID %in% allSigIDs) %>% select(-c("variant.id", "allele.index")) %>% select(snpID, chr, pos, ref, alt, n.obs:low_depth_rate)
  colnames(currResSub)[-c(1:5)] <- paste(colnames(currResSub)[-c(1:5)], sub("\\_.*", "", currBase), sep=".")  
  allResSub <- allResSub %>% left_join(currResSub, by=c("snpID", "chr", "pos", "ref", "alt"))
}

write.csv(allResSub, file=paste0(baseLocation, "MERGED_allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop_p_lt_5e-8.csv"), row.names=FALSE, quote=FALSE)
system(paste0("gzip -f ", baseLocation, "MERGED_allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop_p_lt_5e-8.csv"))

# add OASIS annotation
allAnno <- read.csv(gzfile("/dcl01/mathias1/data/telomere_mtaub/gwas/results/OASISAnno/allAnno.csv.gz"), header= TRUE, stringsAsFactors = FALSE)
allAnno <- allAnno %>% rename(snpID = SNPname)
allResSub <- allResSub %>% left_join(allAnno, by = "snpID")

write.csv(allResSub, file=paste0(baseLocation, "MERGED_allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop_p_lt_5e-8_OASISAnno.csv"), row.names=FALSE, quote=FALSE)
system(paste0("gzip -f ", baseLocation, "MERGED_allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop_p_lt_5e-8_OASISAnno.csv"))


########################################
## want to pull all ancestry results for all 59 variants in Table 1
########################################

allSNPResSub<-read.csv(file = "/dcl01/mathias1/data/telomere_mtaub/gwas/results/ALL_GWASResults/conditional/allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop_p_lt_5e-9_conditionalAllSNPResSub.csv", header=TRUE, stringsAsFactors = FALSE)


for (currBase in allFileBases[-1]){
  currRes <- fread(paste0("gunzip -c ", baseLocation, currBase, ".csv.gz"), header=TRUE, sep=",")
  currResSub <- currRes %>% filter(snpID %in% allSNPResSub$snpID) %>% select(-c("variant.id", "allele.index")) %>% select(snpID, chr, pos, ref, alt, n.obs:low_depth_rate)
  colnames(currResSub)[-c(1:5)] <- paste(colnames(currResSub)[-c(1:5)], sub("\\_.*", "", currBase), sep=".")  
  allSNPResSub <- allSNPResSub %>% left_join(currResSub, by=c("snpID", "chr", "pos", "ref", "alt"))
}

write.csv(allSNPResSub, file=paste0(baseLocation, "MERGED_allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop_p_lt_5e-9_conditionalAllSNPResSub.csv"), row.names=FALSE, quote=FALSE)
system(paste0("gzip -f ", baseLocation, "MERGED_allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop_p_lt_5e-9_conditionalAllSNPResSub.csv"))

# merge in OASIS annotation here
allAnno <- read.csv(gzfile("/dcl01/mathias1/data/telomere_mtaub/gwas/results/OASISAnno/allAnno.csv.gz"), header= TRUE, stringsAsFactors = FALSE)
allAnno <- allAnno %>% rename(snpID = SNPname)
allSNPResSub <- allSNPResSub %>% left_join(allAnno, by = "snpID")

write.csv(allSNPResSub, file=paste0(baseLocation, "MERGED_allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop_p_lt_5e-9_conditionalAllSNPResSub_OASISAnno.csv"), row.names=FALSE, quote=FALSE)
system(paste0("gzip -f ", baseLocation, "MERGED_allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop_p_lt_5e-9_conditionalAllSNPResSub_OASISAnno.csv"))

########################################
## want to pull all results with p < 1e-5 from pooled analysis (not all analyses)
########################################

allRes<-fread(paste0("gunzip -c ", baseLocation, allFileBase, ".csv.gz"), header=TRUE, sep=",")
allResSub <- allRes %>% filter(Score.pval < 1e-5) %>% select(snpID, chr, pos, ref, alt, n.obs.ALL=n.obs, freq.ALL=freq, MAC.ALL=MAC, Score.ALL = Score, Score.SE.ALL=Score.SE,Score.Stat.ALL = Score.Stat,Score.pval.ALL = Score.pval,Est.ALL=Est,Est.SE.ALL=Est.SE,PVE.ALL=PVE, filt.ALL=filt, low_depth_rate.ALL=low_depth_rate)

for (currBase in allFileBases[-1]){
  currRes <- fread(paste0("gunzip -c ", baseLocation, currBase, ".csv.gz"), header=TRUE, sep=",")
  currResSub <- currRes %>% filter(snpID %in% allResSub$snpID) %>% select(-c("variant.id", "allele.index")) %>% select(snpID, chr, pos, ref, alt, n.obs:low_depth_rate)
  colnames(currResSub)[-c(1:5)] <- paste(colnames(currResSub)[-c(1:5)], sub("\\_.*", "", currBase), sep=".")  
  allResSub <- allResSub %>% left_join(currResSub, by=c("snpID", "chr", "pos", "ref", "alt"))
}

write.csv(allResSub, file=paste0(baseLocation, "MERGED_allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop_p_lt_1e-5.csv"), row.names=FALSE, quote=FALSE)
system(paste0("gzip -f ", baseLocation, "MERGED_allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop_p_lt_1e-5.csv"))

# merge in OASIS annotation here
allAnno <- read.csv(gzfile("/dcl01/mathias1/data/telomere_mtaub/gwas/results/OASISAnno/allAnno.csv.gz"), header= TRUE, stringsAsFactors = FALSE)
allAnno <- allAnno %>% rename(snpID = SNPname)
allResSub <- allResSub %>% left_join(allAnno, by = "snpID")

write.csv(allResSub, file=paste0(baseLocation, "MERGED_allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop_p_lt_1e-5_OASISAnno.csv"), row.names=FALSE, quote=FALSE)
system(paste0("gzip -f ", baseLocation, "MERGED_allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop_p_lt_1e-5_OASISAnno.csv"))

## 


########################################
## want to pull all ancestry results for caviar identified credible set SNPs
########################################

caviarSNPs <- read.csv("/dcl01/mathias1/data/telomere_mtaub/gwas/results/CAVIAR/obfc1_credible_sets-2020_05_28.csv", stringsAsFactors = FALSE, header = FALSE)
colnames(caviarSNPs) <- c("rsNumCAVIAR", "snpID")
caviarSNPs <- caviarSNPs %>% mutate(round = ifelse(rsNumCAVIAR == "rs9420907", "Primary", ifelse(rsNumCAVIAR == "rs111447985", "Round1", ifelse(rsNumCAVIAR == "rs112163720", "Round2", ifelse(rsNumCAVIAR == "rs3752946", "Round3", NA)))))

caviarRes<-fread(paste0("gunzip -c ", baseLocation, allFileBase, ".csv.gz"), header=TRUE, sep=",")
caviarResSub <- caviarRes %>% filter(snpID %in% caviarSNPs$snpID) %>% select(-c("variant.id", "allele.index")) %>% select(snpID, chr, pos, ref, alt, n.obs:low_depth_rate)
colnames(caviarResSub)[-c(1:5)] <- paste(colnames(caviarResSub)[-c(1:5)], "ALL", sep=".") 

for (currBase in allFileBases[-1]){
  currRes <- fread(paste0("gunzip -c ", baseLocation, currBase, ".csv.gz"), header=TRUE, sep=",")
  currResSub <- currRes %>% filter(snpID %in% caviarSNPs$snpID) %>% select(-c("variant.id", "allele.index")) %>% select(snpID, chr, pos, ref, alt, n.obs:low_depth_rate)
  colnames(currResSub)[-c(1:5)] <- paste(colnames(currResSub)[-c(1:5)], sub("\\_.*", "", currBase), sep=".")  
  caviarResSub <- caviarResSub %>% left_join(currResSub, by=c("snpID", "chr", "pos", "ref", "alt"))
}


write.csv(caviarResSub, file=paste0(baseLocation, "CAVIAR/MERGED_allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop_CAVIAR_OBFC1.csv"), row.names=FALSE, quote=FALSE)

caviarResSub <- read.csv(paste0(baseLocation, "CAVIAR/MERGED_allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop_CAVIAR_OBFC1.csv"), header = TRUE, stringsAsFactors = FALSE)
caviarResSub<- caviarResSub %>% left_join(caviarSNPs)

## merge together Oasis annotation for all regions around all loci
## read in and merge in conditional results
filesToRead <- system(paste0("ls ", baseLocation, "CAVIAR/*.csv.gz"), intern = TRUE)
filesToRead <- rev(sort(filesToRead))[1:4]
names(filesToRead) <- c("Primary", "Round1", "Round2", "Round3")

for(currRound in names(filesToRead)){
  currRes <- fread(paste0("gunzip -c ", filesToRead[currRound]), header=TRUE, sep=",")
  currResSub <- currRes %>% filter(snpID %in% caviarSNPs$snpID) %>% select(-c("variant.id", "allele.index")) %>% select(snpID, chr, pos, ref, alt, n.obs:PVE)
  colnames(currResSub)[-c(1:5)] <- paste(colnames(currResSub)[-c(1:5)], sub("\\_.*", "", currRound), sep=".")  
  caviarResSub <- caviarResSub %>% left_join(currResSub, by=c("snpID", "chr", "pos", "ref", "alt"))
}


# merge in OASIS annotation here
allAnno <- read.csv(gzfile("/dcl01/mathias1/data/telomere_mtaub/gwas/results/OASISAnno/allAnno.csv.gz"), header= TRUE, stringsAsFactors = FALSE)
allAnno <- allAnno %>% rename(snpID = SNPname)
caviarResSub <- caviarResSub %>% left_join(allAnno, by = "snpID")


write.csv(caviarResSub, file=paste0(baseLocation, "CAVIAR/MERGED_allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop_CAVIAR_OBFC1_OASISAnno.csv"), row.names=FALSE, quote=FALSE)

##########

filesToRead <- system("ls /dcl01/mathias1/data/telomere_mtaub/gwas/results/OASISAnno", intern=TRUE)
allRes <- NULL
for(currFile in filesToRead){
  currRes <- read.csv(gzfile(paste0("/dcl01/mathias1/data/telomere_mtaub/gwas/results/OASISAnno/", currFile)), header=TRUE, stringsAsFactors = FALSE)
  allRes <- rbind(allRes, currRes)
}

write.csv(allRes, file = "/dcl01/mathias1/data/telomere_mtaub/gwas/results/OASISAnno/allAnno.csv", row.names = FALSE, quote = FALSE)
system("gzip -f /dcl01/mathias1/data/telomere_mtaub/gwas/results/OASISAnno/allAnno.csv")


##############################
###### OLD CODE
##############################

## has not been re-done since getting final set of variants to filter out (should be a super-set in any case)
## merge together Oasis annotation

allResSub <- read.csv(gzfile(paste0(baseLocation, "MERGED_allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop_p_lt_5e-8.csv.gz")))

allFileBaseOasis<-"ALL_GWASResults/OASIS_allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop_p_lt_5e-8"
whiteFileBaseOasis<-"White_GWASResults/OASIS_White_allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop_p_lt_5e-8"
blackFileBaseOasis<-"Black_GWASResults/OASIS_Black_allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop_p_lt_5e-8"
asianFileBaseOasis<-"Asian_GWASResults/OASIS_Asian_allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop_p_lt_5e-8"
HLFileBaseOasis<-"HispanicLatino_GWASResults/OASIS_HL_allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop_p_lt_5e-8"
brazilianFileBaseOasis<-"Brazilian_GWASResults/OASIS_Brazilian_allChrs_telomere_adjagesexbatchPCs_minDP0_BRAVODepthDrop_p_lt_5e-8"
samoanFileBaseOasis<-"Samoan_GWASResults/OASIS_Samoan_allChrs_telomere_adjagesexseqctrbatchPCs_minDP0_BRAVODepthDrop_p_lt_5e-8"
allFileBasesOasis<-c(allFileBaseOasis, whiteFileBaseOasis, blackFileBaseOasis, asianFileBaseOasis, HLFileBaseOasis, brazilianFileBaseOasis, samoanFileBaseOasis)

allOasis<-read.csv(paste0(baseLocation, allFileBasesOasis[1], ".csv"), header=TRUE, stringsAsFactors = FALSE)
for(currBase in allFileBasesOasis[-1]){
  currRes <- read.csv(paste0(baseLocation, currBase, ".csv"), header=TRUE, stringsAsFactors = FALSE)
  currRes <- currRes %>% filter(!SNPname %in% allOasis$SNPname)
  allOasis <-rbind(allOasis, currRes)
}
allOasis <- allOasis %>% rename(snpID = SNPname)

allResSub <- allResSub %>% left_join(allOasis, by="snpID")

write.csv(allResSub, file=paste0(baseLocation, "MERGED_allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop_p_lt_5e-8_OASISAnno.csv"), row.names=FALSE, quote=FALSE)
system(paste0("gzip -f ", baseLocation, "MERGED_allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop_p_lt_5e-8_OASISAnno.csv"))

