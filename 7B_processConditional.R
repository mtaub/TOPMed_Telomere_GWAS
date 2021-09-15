

library(Gviz)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(Homo.sapiens)
library(data.table)
library(dplyr)
library(RColorBrewer)
library(tidyr)

highlightColors<-brewer.pal(n = 8, name = "Dark2")

snpIDsToDrop <- scan("/dcl01/mathias1/data/telomere_mtaub/gwas/results/ALL_GWASResults/snpIDsToDrop_BRAVO.txt", what="character")
all( c("4:163128904:T:G", "8:73046129:G:A", "4:171697861:A:C", "7:125154536:A:G", "22:40475409:T:C", "20:63667021:A:T", "10:103911091:A:T") %in% snpIDsToDrop)
snpIDsToDrop_lowDepth <- scan("/dcl01/mathias1/data/telomere_mtaub/gwas/results/ALL_GWASResults/snpIDsToDrop_lowDepth.txt", what="character")

pathToResults<-"/dcl01/mathias1/data/telomere_mtaub/gwas/results/ALL_GWASResults/"

condRes<-read.csv("/dcl01/mathias1/data/telomere_mtaub/gwas/results/ALL_GWASResults/conditional/allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop_p_lt_5e-9_peakSNPsByChromosome_round7PostBRAVO.csv", header=TRUE, stringsAsFactors = FALSE)

## need to fix the row for chr22
newPeaks<- read.csv("/dcl01/mathias1/data/telomere_mtaub/gwas/results/ALL_GWASResults/allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop_p_lt_5e-8_peakSNPsByLocus.csv", header=TRUE, stringsAsFactors = FALSE)
cond22<-fread("gunzip -c /dcl01/mathias1/data/telomere_mtaub/gwas/results/ALL_GWASResults/conditional/chr22_telomere_adjagesexstudyseqctrbatchPCs_minDP0_cond_22-50532618-T-C.csv.gz", header= TRUE, sep=",")
cond22Sub<-cond22 %>% filter(!snpID %in% c(snpIDsToDrop, snpIDsToDrop_lowDepth), Score.pval < 0.1)
cond22SubTop<-cond22Sub %>% arrange(Score.pval) %>% head(1)

condRes[condRes$chr == 22,c("pos_round1", "snpID_round1", "Score.pval_round1")] <- newPeaks[newPeaks$chr == 22,c("pos", "snpID", "Score.pval")]
condRes[condRes$chr == 22,c("pos_round2", "snpID_round2", "Score.pval_round2")] <- cond22SubTop[,c("pos", "snpID", "Score.pval")]
condRes[condRes$chr == 22,c("toCondition_round2")] <- "No more rounds"
condRes[condRes$chr == 22,c("pos_round3", "snpID_round3", "Score.pval_round3","toCondition_round3","pos_round4", "snpID_round4", "Score.pval_round4","toCondition_round4"  )] <- NA

write.csv(condRes, file="/dcl01/mathias1/data/telomere_mtaub/gwas/results/ALL_GWASResults/conditional/allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop_p_lt_5e-9_peakSNPsByChromosome_round7PostBRAVO_chr22Fix.csv", row.names=FALSE, quote=TRUE)


setsByRound<-condRes %>% rename(toCondition_round1 = snpID_round1) %>% select(chr, starts_with("toCondition_")) %>% pivot_longer(-chr, names_to="round", values_to ="toCondition") %>% mutate(round = sub("toCondition_", "", round))

toProcess<-setsByRound %>% filter(!toCondition == "No more rounds", !is.na(toCondition)) %>% arrange(chr, desc(round)) %>% filter(!duplicated(chr)) %>% arrange(as.numeric(chr))

allSNPRes<-do.call(rbind, lapply(toProcess$toCondition, makePlotAndPullNext))

write.csv(allSNPRes, file="/dcl01/mathias1/data/telomere_mtaub/gwas/results/ALL_GWASResults/conditional/allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop_p_lt_5e-9_conditionalAllSNPRes.csv", row.names=FALSE, quote=TRUE)

allSNPResSub <- allSNPRes %>% filter(Score.pval < 5e-9 | round == "Primary")
allSNPResSub <- allSNPResSub %>% mutate(roundName = ifelse(round == "Primary", "Primary", "Cond"), roundNum = ifelse(round == "Primary", 0, as.numeric(sub("cond", "", round))))
allSNPResSub <- allSNPResSub %>% pivot_wider(names_from = roundName, values_from = c(Score, Score.SE, Score.Stat, Score.pval, Est, Est.SE, PVE, round, roundNum))
allSNPResSub <- allSNPResSub %>% mutate(roundNum_Cond = ifelse(is.na(roundNum_Cond), 0, roundNum_Cond)) %>% arrange(as.numeric(chr), roundNum_Cond)

write.csv(allSNPResSub, file = "/dcl01/mathias1/data/telomere_mtaub/gwas/results/ALL_GWASResults/conditional/allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop_p_lt_5e-9_conditionalAllSNPResSub.csv", row.names=FALSE, quote=TRUE)




## want this function to:
## pull results for all SNPs from that chr that were ever conditioned on from the main results
## pull complete results from preceding round of conditional for next highest SNP
## create final version of manhattan-style plot
makePlotAndPullNext <- function(condSNP){
  gc()
  
  chr<-paste0("chr", gsub("\\:.*", "", condSNP))
  condSNPPos<-as.numeric(strsplit(condSNP, split=":")[[1]][2])
  gen <- "hg38"
  
  #Genome Axis Track
  gtrack <- GenomeAxisTrack(littleTicks=TRUE)
  
  nonCond<-fread(paste0("gunzip -c ", pathToResults, "/", chr, "_telomere_adjagesexstudyseqctrbatchPCs_minDP0.csv.gz"), header=TRUE, sep=",")
  nonCondSub<-nonCond %>% filter(!snpID %in% c(snpIDsToDrop, snpIDsToDrop_lowDepth), Score.pval < 0.1)

  condSNPVec <- unlist(strsplit(condSNP, split=";"))
  nSNPs<-length(condSNPVec)
  SNPTags<-gsub("\\:", "\\-", condSNPVec)
  SNPPos<-as.numeric(sapply(strsplit(SNPTags, split="-"), function(x) x[2]))

  ## pull results for all conditioned SNPs from main results
  allSNPRes<- nonCond %>% filter(snpID %in% condSNPVec)
  allSNPRes$round <- "Primary"
  
  ## for now this will work, but will need to update once there is more than one file per chr
  ## read in appropriate conditional data files
  ## need to modify so that it sorts snps by position since that is the order in the file name
  condResList <- NULL
  for (i in 1:nSNPs){
    gc()
    currCondRes<-fread(paste0("gunzip -c ", pathToResults, "/conditional/", chr, "_telomere_adjagesexstudyseqctrbatchPCs_minDP0_cond_", paste(SNPTags[1:i], collapse="_"), ".csv.gz"), header=TRUE, sep=",")
    currCondResSub<-currCondRes %>% filter(!snpID %in% c(snpIDsToDrop, snpIDsToDrop_lowDepth), Score.pval < 0.1, !snpID %in% condSNPVec[1:i])
    condResList <- c(condResList, list(currCondResSub))
    allSNPRes <- rbind(allSNPRes, data.frame(currCondResSub %>% filter(Score.pval == min(Score.pval)), round = paste0("cond", i)))
  }
  yLimMaxCond <- max(sapply(condResList, function(x) max(-log10(x$Score.pval))))
  yLimTop<-ceiling(max(max(-log10(nonCondSub$Score.pval)), yLimMaxCond))
    
  dtrack_nonCond <- DataTrack(data=-log10(nonCondSub$Score.pval), start=nonCondSub$pos-1, end= nonCondSub$pos, genome="hg38", chromosome=chr, name= paste0(chr, " Non conditional"), type=c("p"), col="black", baseline = -log10(5e-9), col.baseline="chartreuse3", ylim=c(0, yLimTop), grid = TRUE, grid.col = "gray84")
  
  dtrack_condList <- NULL
  for (i in 1:nSNPs){
    currCondResSub <- condResList[[i]]
    dtrack_currCondResSub <- DataTrack(data=-log10(currCondResSub$Score.pval), start=currCondResSub$pos-1, end= currCondResSub$pos, genome="hg38", chromosome=chr, name= paste0("Conditional on: ", paste(condSNPVec[1:i], collapse=";\n")), type=c("p"), col="black", baseline = -log10(5e-9), col.baseline="chartreuse3", ylim=c(0, yLimTop), grid = TRUE, grid.col = "gray84")
    dtrack_condList <- c(dtrack_condList, list(dtrack_currCondResSub))
  }
  

  ht <- HighlightTrack(trackList=c(list(dtrack_nonCond), dtrack_condList, list(gtrack)), start=c(SNPPos-10000), width=20000, chromosome=chr, col=highlightColors[1:nSNPs])
  

  #Create PDF
  png(paste0("/dcl01/mathias1/data/telomere_mtaub/gwas/results/ALL_GWASResults/conditional/", chr, "_cond_", paste(SNPTags[1:i], collapse="_"), "_Final.png"),height=1200,width=1600)
  #Plot Data
  plotTracks(list(ht),background.title="darkgray")
  #End
  dev.off()
  
  return(allSNPRes)
}


