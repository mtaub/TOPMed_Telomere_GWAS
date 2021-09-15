
library(tidyverse)
library(readxl)

## starting from single chromosome files, want to merge them all together
## also want to filter out positions identified as pathogenic by Rasika in BRAVO

##############################
## SNPs to drop from BRAVO
##############################

setwd("/Users/mtaub/Research/OneDrive - Johns Hopkins/telomere/gwas")

bravoRes<-NULL
for(currChr in c(paste("Chr", c(1:22, "X")))){
  bravoRes<-rbind(bravoRes, read_excel("./results/ALL_GWASResults/OASIS QC/Chromosome_BRAVO_confirm.xlsx", sheet=currChr))
}
bravoResExtra<-read_excel("./results/ALL_GWASResults/OASIS QC/Chromosome_BRAVO_confirm.xlsx", sheet="Chr N")
bravoResExtra2<-read_excel("./results/ALL_GWASResults/OASIS QC/Chromosome_BRAVO_confirm.xlsx", sheet="Extras looked up")
bravoResConditional<-read_excel("./results/ALL_GWASResults/OASIS QC/Chromosome_BRAVO_confirm.xlsx", sheet="Conditional Analysis Lookups")
bravoResStratified<-read_excel("./results/ALL_GWASResults/OASIS QC/Chromosome_BRAVO_confirm.xlsx", sheet="Stratified Anlaysis Lookups")

snpIDsToDrop<-bravoRes %>% filter(grepl("fail", BRAVO) | grepl("not found", BRAVO)) %>% pull(SNPname)
snpIDsToDrop <- c(snpIDsToDrop, bravoResExtra$snpID[grepl("fail", bravoResExtra$BRAVO)], bravoResExtra2$SNPname[grepl("fail", bravoResExtra2$BRAVO)], bravoResConditional$VarID[grepl("fail", bravoResConditional$BRAVO)], bravoResStratified$`x...1`[grepl("fail", bravoResStratified$BRAVO)])
write(snpIDsToDrop, file="./results/ALL_GWASResults/snpIDsToDrop_BRAVO.txt")

system("scp './results/ALL_GWASResults/snpIDsToDrop_BRAVO.txt' jhp:/dcl01/mathias1/data/telomere_mtaub/gwas/results/ALL_GWASResults/")


##############################
## Merge all chromosomes and drop SNPs
##############################

## First on local machine, copy files to cluster:
# rsync -hv --progress ./ALL_GWASResults/chr*_telomere_adjagesexstudyseqctrbatchPCs_minDP0.csv.gz jhp:/dcl01/mathias1/data/telomere_mtaub/gwas/results/ALL_GWASResults
# rsync -hv --progress ./Black_GWASResults/Black_chr*_telomere_adjagesexstudyseqctrbatchPCs_minDP0.csv.gz jhp:/dcl01/mathias1/data/telomere_mtaub/gwas/results/Black_GWASResults
# rsync -hv --progress ./White_GWASResults/White_chr*_telomere_adjagesexstudyseqctrbatchPCs_minDP0.csv.gz jhp:/dcl01/mathias1/data/telomere_mtaub/gwas/results/White_GWASResults
#rsync -hv --progress ./Asian_GWASResults/Asian_chr*_telomere_adjagesexstudyseqctrbatchPCs_minDP0.csv.gz jhp:/dcl01/mathias1/data/telomere_mtaub/gwas/results/Asian_GWASResults
#rsync -hv --progress ./Brazilian_GWASResults/Brazilian_chr*_telomere_adjagesexbatchPCs_minDP0.csv.gz jhp:/dcl01/mathias1/data/telomere_mtaub/gwas/results/Brazilian_GWASResults
#rsync -hv --progress ./Samoan_GWASResults/Samoan_chr*_telomere_adjagesexseqctrbatchPCs_minDP0.csv.gz jhp:/dcl01/mathias1/data/telomere_mtaub/gwas/results/Samoan_GWASResults
#rsync -hv --progress ./HispanicLatino_GWASResults/HL_chr*_telomere_adjagesexstudyseqctrbatchPCs_minDP0.csv.gz jhp:/dcl01/mathias1/data/telomere_mtaub/gwas/results/HispanicLatino_GWASResults


## now on JHPCE cluster
library(dplyr)
library(data.table)

mergeChrs<-function(pathToResults, popTag, baseFileName){
  system(paste0("gunzip -c ", pathToResults, "/", popTag, "chr1_telomere_", baseFileName, ".csv.gz > ", pathToResults, "/", popTag, "allChrs_telomere_", baseFileName, ".csv"))
  for (currChr in c(2:22, "X")){
    system(paste0("gunzip -c ", pathToResults, "/", popTag, "chr", currChr, "_telomere_", baseFileName, ".csv.gz | tail -n +2 >> ", pathToResults, "/", popTag, "allChrs_telomere_", baseFileName, ".csv"))
  }
  system(paste0("gzip ", pathToResults, "/", popTag, "allChrs_telomere_", baseFileName, ".csv"))
}

mergeChrs("/dcl01/mathias1/data/telomere_mtaub/gwas/results/Samoan_GWASResults", "Samoan_", "adjagesexseqctrbatchPCs_minDP0")
mergeChrs("/dcl01/mathias1/data/telomere_mtaub/gwas/results/White_GWASResults", "White_", "adjagesexstudyseqctrbatchPCs_minDP0")
mergeChrs("/dcl01/mathias1/data/telomere_mtaub/gwas/results/Asian_GWASResults", "Asian_", "adjagesexstudyseqctrbatchPCs_minDP0")
mergeChrs("/dcl01/mathias1/data/telomere_mtaub/gwas/results/HispanicLatino_GWASResults", "HL_", "adjagesexstudyseqctrbatchPCs_minDP0")
mergeChrs("/dcl01/mathias1/data/telomere_mtaub/gwas/results/Black_GWASResults", "Black_", "adjagesexstudyseqctrbatchPCs_minDP0")
mergeChrs("/dcl01/mathias1/data/telomere_mtaub/gwas/results/Brazilian_GWASResults", "Brazilian_", "adjagesexbatchPCs_minDP0")
#mergeChrs("/dcl01/mathias1/data/telomere_mtaub/gwas/results/ALL_GWASResults", "", "adjagesexstudyseqctrbatchPCs_minDP0")


## Identify SNPs with low depth to filter out
#allLoci <- fread(paste0("gunzip -c /dcl01/mathias1/data/telomere_mtaub/gwas/results/ALL_GWASResults/allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0.csv.gz"), header=TRUE, sep=",")
#snpIDsToDrop_lowDepth<- allLoci %>% filter(low_depth_rate > 0.1) %>% pull(snpID)
#write(snpIDsToDrop_lowDepth, file="/dcl01/mathias1/data/telomere_mtaub/gwas/results/ALL_GWASResults/snpIDsToDrop_lowDepth.txt")


createCleanResults<-function(pathToResults, popTag, baseFileName){
  ## Full data set
  # 162526185 variants
  allLoci<-fread(paste0("gunzip -c ", pathToResults, "/", popTag, "allChrs_telomere_", baseFileName, ".csv.gz"), header=TRUE, sep=",")
  
  ## remove SNPs identified in BRAVO
  snpIDsToDrop <- scan("/dcl01/mathias1/data/telomere_mtaub/gwas/results/ALL_GWASResults/snpIDsToDrop_BRAVO.txt", what="character")
  snpIDsToDrop_lowDepth <- scan("/dcl01/mathias1/data/telomere_mtaub/gwas/results/ALL_GWASResults/snpIDsToDrop_lowDepth.txt", what="character")
  allLociSub<-allLoci %>% filter(!snpID %in% c(snpIDsToDrop, snpIDsToDrop_lowDepth))
  
  
  ## write clean data set for input into OASIS, with just columns needed for OASIS (to make smaller file)
  #write.csv(allLociSub %>% select(variant.id, chr, pos, allele.index, n.obs, freq, Score, Score.SE, Score.Stat, Score.pval, ref, alt), file=paste0(pathToResults, "/", popTag, "allChrs_telomere_", baseFileName, "_BRAVODepthDrop_forOASIS.csv"), row.names=FALSE, quote=FALSE)
  #system(paste0("gzip -f ", pathToResults, "/", popTag, "allChrs_telomere_", baseFileName, "_BRAVODepthDrop_forOASIS.csv"))

  ## write clean data set for input into OASIS, with just columns needed for OASIS (to make smaller file)
  write.csv(allLociSub %>% select(CHR=chr, POS=pos, REF=ref, ALT=alt, AF=freq, N=n.obs, BETA=Est, SE=Est.SE, P=Score.pval), file=paste0(pathToResults, "/", popTag, "allChrs_telomere_", baseFileName, "_BRAVODepthDrop_forOASISNewFormat.csv"), row.names=FALSE, quote=FALSE)
  system(paste0("gzip -f ", pathToResults, "/", popTag, "allChrs_telomere_", baseFileName, "_BRAVODepthDrop_forOASISNewFormat.csv"))
  
  ## write all sig positions, clean data set
  write.csv(allLociSub %>% filter(Score.pval < 5e-8), file=paste0(pathToResults, "/", popTag, "allChrs_telomere_", baseFileName, "_BRAVODepthDrop_p_lt_5e-8.csv"), row.names=FALSE, quote=FALSE)
  system(paste0("gzip -f ", pathToResults, "/", popTag, "allChrs_telomere_", baseFileName, "_BRAVODepthDrop_p_lt_5e-8.csv"))
  
  ## now write full file with all results and all columns, with dropped variants
  write.csv(allLociSub, file=paste0(pathToResults, "/", popTag, "allChrs_telomere_", baseFileName, "_BRAVODepthDrop.csv"), row.names=FALSE, quote=FALSE)
  system(paste0("gzip -f ", pathToResults, "/", popTag, "allChrs_telomere_", baseFileName, "_BRAVODepthDrop.csv"))
  
  
}

createCleanResults("/dcl01/mathias1/data/telomere_mtaub/gwas/results/ALL_GWASResults", "", "adjagesexstudyseqctrbatchPCs_minDP0")
createCleanResults("/dcl01/mathias1/data/telomere_mtaub/gwas/results/White_GWASResults", "White_", "adjagesexstudyseqctrbatchPCs_minDP0")
createCleanResults("/dcl01/mathias1/data/telomere_mtaub/gwas/results/Asian_GWASResults", "Asian_", "adjagesexstudyseqctrbatchPCs_minDP0")
createCleanResults("/dcl01/mathias1/data/telomere_mtaub/gwas/results/HispanicLatino_GWASResults", "HL_", "adjagesexstudyseqctrbatchPCs_minDP0")
createCleanResults("/dcl01/mathias1/data/telomere_mtaub/gwas/results/Black_GWASResults", "Black_", "adjagesexstudyseqctrbatchPCs_minDP0")
createCleanResults("/dcl01/mathias1/data/telomere_mtaub/gwas/results/Brazilian_GWASResults", "Brazilian_", "adjagesexbatchPCs_minDP0")
createCleanResults("/dcl01/mathias1/data/telomere_mtaub/gwas/results/Samoan_GWASResults", "Samoan_", "adjagesexseqctrbatchPCs_minDP0")


################################
######## OLD CODE
################################



## select some random variants for Rasika:
#forRasBravoCheck<-rbind(allLoci %>% filter(MAC < 30, Score.pval < 5e-8) %>% sample_n(20), allLoci %>% filter(MAC < 30, Score.pval > 5e-8, Score.pval < 1e-5) %>% sample_n(20), allLoci %>% filter(MAC < 30, Score.pval > 1e-5) %>% sample_n(20))
#write.csv(forRasBravoCheck, file="~/Research/telomere/newManuscript/results/forRasBravoCheck.csv", row.names=FALSE, quote=FALSE)

## after reading in full new minDP0 data set, need to pull positions with p<5e-8 to compare to previous list (so we know what still needs to be looked up in BRAVO)
write.csv(allLoci %>% filter(Score.pval < 5e-8), file=paste0(pathToResults, "/", popTag, "allChrs_telomere_", baseFileName, "_p_lt_5e-8.csv"), row.names=FALSE, quote=FALSE)
system(paste0("gzip ", pathToResults, "/", popTag, "allChrs_telomere_", baseFileName, "_p_lt_5e-8.csv"))

## read in old Bravo data
newRes<-read.csv(gzfile("~/Research/OneDrive/telomere/gwas/results/ALL_GWASResults/allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_p_lt_5e-8.csv.gz"))
write.csv(newRes %>% filter(!snpID %in% c(bravoRes$SNPname, bravoResExtra$snpID)), file="~/Research/OneDrive/telomere/gwas/results/ALL_GWASResults/OASIS QC/moreToCheck_minDP0.csv", row.names=FALSE, quote=FALSE)
