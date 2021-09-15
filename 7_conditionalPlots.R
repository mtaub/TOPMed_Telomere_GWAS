

library(Gviz)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(Homo.sapiens)
library(data.table)
library(dplyr)
library(RColorBrewer)

highlightColors<-brewer.pal(n = 8, name = "Dark2")

snpIDsToDrop <- scan("/dcl01/mathias1/data/telomere_mtaub/gwas/results/ALL_GWASResults/snpIDsToDrop_BRAVO.txt", what="character")
all( c("4:163128904:T:G", "8:73046129:G:A", "4:171697861:A:C", "7:125154536:A:G", "22:40475409:T:C", "20:63667021:A:T", "10:103911091:A:T") %in% snpIDsToDrop)
snpIDsToDrop_lowDepth <- scan("/dcl01/mathias1/data/telomere_mtaub/gwas/results/ALL_GWASResults/snpIDsToDrop_lowDepth.txt", what="character")

pathToResults<-"/dcl01/mathias1/data/telomere_mtaub/gwas/results/ALL_GWASResults/"


peakPerChr<- read.csv("/dcl01/mathias1/data/telomere_mtaub/gwas/results/ALL_GWASResults/allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop_p_lt_5e-9_peakSNPsByChromosome.csv", header=TRUE, stringsAsFactors = FALSE)


nextToPull<-do.call(rbind, lapply(peakPerChr$snpID, makePlotAndPullNext))
makePlotAndPullNext("22:40128558:G:A")

peakPerChr <- peakPerChr %>% select(chr, pos_round1=pos, snpID_round1 = snpID, Score.pval_round1 = Score.pval)
nextToPull <- nextToPull %>% select(chr, pos_round2 = pos, snpID_round2 = snpID, Score.pval_round2 = Score.pval)
peakPerChr <- left_join(peakPerChr, nextToPull)
peakPerChr <- peakPerChr %>% mutate(toCondition_round2 = ifelse(Score.pval_round2 < 5e-9, paste(snpID_round1, snpID_round2, sep=";"), "No more rounds"))
write.csv(peakPerChr, file="/dcl01/mathias1/data/telomere_mtaub/gwas/results/ALL_GWASResults/conditional/allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop_p_lt_5e-9_peakSNPsByChromosome_round2.csv", row.names=FALSE, quote=TRUE)

peakPerChr2<-read.csv("/dcl01/mathias1/data/telomere_mtaub/gwas/results/ALL_GWASResults/conditional/allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop_p_lt_5e-9_peakSNPsByChromosome_round2PostBRAVO.csv", header=TRUE, stringsAsFactors = FALSE)
peakPerChr2Sub <- peakPerChr2 %>% filter(!toCondition_round2 == "No more rounds")
nextToPull2<-do.call(rbind, lapply(peakPerChr2Sub$toCondition_round2, makePlotAndPullNext))

nextToPull2 <- nextToPull2 %>% select(chr, pos_round3 = pos, snpID_round3 = snpID, Score.pval_round3 = Score.pval)
peakPerChr2 <- left_join(peakPerChr2, nextToPull2)
peakPerChr2 <- peakPerChr2 %>% mutate(toCondition_round3 = ifelse(Score.pval_round3 < 5e-9, paste(snpID_round1, snpID_round2, snpID_round3, sep=";"), "No more rounds"))
write.csv(peakPerChr2, file="/dcl01/mathias1/data/telomere_mtaub/gwas/results/ALL_GWASResults/conditional/allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop_p_lt_5e-9_peakSNPsByChromosome_round3.csv", row.names=FALSE, quote=TRUE)

peakPerChr3<-read.csv("/dcl01/mathias1/data/telomere_mtaub/gwas/results/ALL_GWASResults/conditional/allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop_p_lt_5e-9_peakSNPsByChromosome_round3PostBRAVO.csv", header=TRUE, stringsAsFactors = FALSE)
peakPerChr3Sub <- peakPerChr3 %>% filter(!toCondition_round3 == "No more rounds")
nextToPull3<-do.call(rbind, lapply(peakPerChr3Sub$toCondition_round3, makePlotAndPullNext))

nextToPull3 <- nextToPull3 %>% select(chr, pos_round4 = pos, snpID_round4 = snpID, Score.pval_round4 = Score.pval) %>% mutate(chr=as.character(chr))
peakPerChr3 <- left_join(peakPerChr3, nextToPull3)
peakPerChr3 <- peakPerChr3 %>% mutate(toCondition_round4 = ifelse(Score.pval_round4 < 5e-9, paste(snpID_round1, snpID_round2, snpID_round3, snpID_round4, sep=";"), "No more rounds"))
write.csv(peakPerChr3, file="/dcl01/mathias1/data/telomere_mtaub/gwas/results/ALL_GWASResults/conditional/allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop_p_lt_5e-9_peakSNPsByChromosome_round4.csv", row.names=FALSE, quote=TRUE)

peakPerChr4<-read.csv("/dcl01/mathias1/data/telomere_mtaub/gwas/results/ALL_GWASResults/conditional/allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop_p_lt_5e-9_peakSNPsByChromosome_round4PostBRAVO.csv", header=TRUE, stringsAsFactors = FALSE)
peakPerChr4Sub <- peakPerChr4 %>% filter(!toCondition_round4 == "No more rounds", chr %in% c(3,5,14,16,18,20))
nextToPull41<-do.call(rbind, lapply(peakPerChr4Sub$toCondition_round4, makePlotAndPullNext))
peakPerChr4Sub2 <- peakPerChr4 %>% filter(!toCondition_round4 == "No more rounds", chr %in% c(4,10))
nextToPull42<-do.call(rbind, lapply(peakPerChr4Sub2$toCondition_round4, makePlotAndPullNext))

nextToPull4<-rbind(nextToPull41, nextToPull42)

nextToPull4 <- nextToPull4 %>% select(chr, pos_round5 = pos, snpID_round5 = snpID, Score.pval_round5 = Score.pval) %>% mutate(chr=as.character(chr))
peakPerChr4 <- left_join(peakPerChr4, nextToPull4)
peakPerChr4 <- peakPerChr4 %>% mutate(toCondition_round5 = ifelse(Score.pval_round5 < 5e-9, paste(snpID_round1, snpID_round2, snpID_round3, snpID_round4, snpID_round5, sep=";"), "No more rounds"))
write.csv(peakPerChr4, file="/dcl01/mathias1/data/telomere_mtaub/gwas/results/ALL_GWASResults/conditional/allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop_p_lt_5e-9_peakSNPsByChromosome_round5.csv", row.names=FALSE, quote=TRUE)


peakPerChr5<-read.csv("/dcl01/mathias1/data/telomere_mtaub/gwas/results/ALL_GWASResults/conditional/allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop_p_lt_5e-9_peakSNPsByChromosome_round5PostBRAVO.csv", header=TRUE, stringsAsFactors = FALSE)
peakPerChr5Sub <- peakPerChr5 %>% filter(!toCondition_round5 == "No more rounds")
nextToPull5<-do.call(rbind, lapply(peakPerChr5Sub$toCondition_round5, makePlotAndPullNext))

nextToPull5 <- nextToPull5 %>% select(chr, pos_round6 = pos, snpID_round6 = snpID, Score.pval_round6 = Score.pval) %>% mutate(chr=as.character(chr))
peakPerChr5 <- left_join(peakPerChr5, nextToPull5)
peakPerChr5 <- peakPerChr5 %>% mutate(toCondition_round6 = ifelse(Score.pval_round6 < 5e-9, paste(snpID_round1, snpID_round2, snpID_round3, snpID_round4, snpID_round5, snpID_round6, sep=";"), "No more rounds"))
write.csv(peakPerChr5, file="/dcl01/mathias1/data/telomere_mtaub/gwas/results/ALL_GWASResults/conditional/allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop_p_lt_5e-9_peakSNPsByChromosome_round6.csv", row.names=FALSE, quote=TRUE)

peakPerChr6<-read.csv("/dcl01/mathias1/data/telomere_mtaub/gwas/results/ALL_GWASResults/conditional/allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop_p_lt_5e-9_peakSNPsByChromosome_round6PostBRAVO.csv", header=TRUE, stringsAsFactors = FALSE)
peakPerChr6Sub <- peakPerChr6 %>% filter(!toCondition_round6 == "No more rounds")
nextToPull6<-do.call(rbind, lapply(peakPerChr6Sub$toCondition_round6, makePlotAndPullNext))

nextToPull6 <- nextToPull6 %>% select(chr, pos_round7 = pos, snpID_round7 = snpID, Score.pval_round7 = Score.pval) %>% mutate(chr=as.character(chr))
peakPerChr6 <- left_join(peakPerChr6, nextToPull6)
peakPerChr6 <- peakPerChr6 %>% mutate(toCondition_round7 = ifelse(Score.pval_round7 < 5e-9, paste(snpID_round1, snpID_round2, snpID_round3, snpID_round4, snpID_round5, snpID_round6, snpID_round7, sep=";"), "No more rounds"))
write.csv(peakPerChr6, file="/dcl01/mathias1/data/telomere_mtaub/gwas/results/ALL_GWASResults/conditional/allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop_p_lt_5e-9_peakSNPsByChromosome_round7.csv", row.names=FALSE, quote=TRUE)


peakPerChr7<-read.csv("/dcl01/mathias1/data/telomere_mtaub/gwas/results/ALL_GWASResults/conditional/allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop_p_lt_5e-9_peakSNPsByChromosome_round7PostBRAVO.csv", header=TRUE, stringsAsFactors = FALSE)
peakPerChr7Sub <- peakPerChr7 %>% filter(!toCondition_round7 == "No more rounds")
nextToPull7<-do.call(rbind, lapply(peakPerChr7Sub$toCondition_round7, makePlotAndPullNext))

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

  ## for now this will work, but will need to update once there is more than one file per chr
  ## read in appropriate conditional data files
  ## need to modify so that it sorts snps by position since that is the order in the file name
  condResList <- NULL
  for (i in 1:nSNPs){
    gc()
    currCondRes<-fread(paste0("gunzip -c ", pathToResults, "/conditional/", chr, "_telomere_adjagesexstudyseqctrbatchPCs_minDP0_cond_", paste(SNPTags[1:i], collapse="_"), ".csv.gz"), header=TRUE, sep=",")
    currCondResSub<-currCondRes %>% filter(!snpID %in% c(snpIDsToDrop, snpIDsToDrop_lowDepth), Score.pval < 0.1, !snpID %in% condSNPVec[1:i])
    condResList <- c(condResList, list(currCondResSub))
  }
  yLimMaxCond <- max(sapply(condResList, function(x) max(-log10(x$Score.pval))))
  yLimTop<-ceiling(max(max(-log10(nonCondSub$Score.pval)), yLimMaxCond))
    
  dtrack_nonCond <- DataTrack(data=-log10(nonCondSub$Score.pval), start=nonCondSub$pos-1, end= nonCondSub$pos, genome="hg38", chromosome=chr, name= paste0(chr, " Non conditional"), type=c("p","g"), col="black", baseline = -log10(5e-9), col.baseline="blue", ylim=c(0, yLimTop))
  
  dtrack_condList <- NULL
  for (i in 1:nSNPs){
    currCondResSub <- condResList[[i]]
    dtrack_currCondResSub <- DataTrack(data=-log10(currCondResSub$Score.pval), start=currCondResSub$pos-1, end= currCondResSub$pos, genome="hg38", chromosome=chr, name= paste0(chr, " Conditional on ", paste(condSNPVec[1:i], collapse=";")), type=c("p","g"), col="black", baseline = -log10(5e-9), col.baseline="blue", ylim=c(0, yLimTop))
    dtrack_condList <- c(dtrack_condList, list(dtrack_currCondResSub))
  }
  

  ht <- HighlightTrack(trackList=c(list(dtrack_nonCond), dtrack_condList, list(gtrack)), start=c(SNPPos-10000), width=20000, chromosome=chr, col=highlightColors[1:nSNPs])
  

  #Create PDF
  png(paste0("/dcl01/mathias1/data/telomere_mtaub/gwas/results/ALL_GWASResults/conditional/", chr, "_cond_", paste(SNPTags[1:i], collapse="_"), ".png"),height=1200,width=1600)
  #Plot Data
  plotTracks(list(ht),background.title="darkgray")
  #End
  dev.off()
  
  newPeak <- condResList[[nSNPs]] %>% filter(Score.pval == min(Score.pval)) %>% select(snpID, chr, pos, Score.pval)
  return(newPeak)
}


oldPeakPerChr<-read.csv("/dcl01/mathias1/data/telomere_mtaub/gwas/results/ALL_GWASResults/conditional/allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop_p_lt_5e-9_peakSNPsByChromosome_round2.csv", header=TRUE, stringsAsFactors = FALSE)

## need to pull new second round SNPs for chr4 and chr8
## for both, want to take the second best one based on Rasika's BRAVO assessment
cond1<-fread(paste0("gunzip -c ", pathToResults, "/conditional/chr4_*.csv.gz"), header=TRUE, sep=",")
cond1Sub<-cond1 %>% filter(!snpID %in% c(snpIDsToDrop, snpIDsToDrop_lowDepth), Score.pval < 0.1)
toInclude<-(cond1Sub %>% arrange(Score.pval) %>% select(snpID, chr, pos, Score.pval, freq, MAC) %>% head())[2,]

cond1_chr8<-fread(paste0("gunzip -c ", pathToResults, "/conditional/chr8_*.csv.gz"), header=TRUE, sep=",")
cond1_chr8Sub<-cond1_chr8 %>% filter(!snpID %in% c(snpIDsToDrop, snpIDsToDrop_lowDepth), Score.pval < 0.1)
cond1_chr8Sub %>% arrange(Score.pval) %>% select(snpID, chr, pos, Score.pval, freq, MAC) %>% head()
toInclude_chr8<-(cond1_chr8Sub %>% arrange(Score.pval) %>% select(snpID, chr, pos, Score.pval, freq, MAC) %>% head())[2,]

oldPeakPerChr[oldPeakPerChr$chr == 4, "pos_round2"]<-toInclude[,"pos"]
oldPeakPerChr[oldPeakPerChr$chr == 4, "snpID_round2"]<-toInclude[,"snpID"]
oldPeakPerChr[oldPeakPerChr$chr == 4, "Score.pval_round2"]<-toInclude[,"Score.pval"]
oldPeakPerChr[oldPeakPerChr$chr == 4, "toCondition_round2"]<-paste(oldPeakPerChr[oldPeakPerChr$chr == 4, "snpID_round1"], toInclude[,"snpID"], sep=";")

oldPeakPerChr[oldPeakPerChr$chr == 8, "pos_round2"]<-toInclude_chr8[,"pos"]
oldPeakPerChr[oldPeakPerChr$chr == 8, "snpID_round2"]<-toInclude_chr8[,"snpID"]
oldPeakPerChr[oldPeakPerChr$chr == 8, "Score.pval_round2"]<-toInclude_chr8[,"Score.pval"]
oldPeakPerChr[oldPeakPerChr$chr == 8, "toCondition_round2"]<-paste(oldPeakPerChr[oldPeakPerChr$chr == 8, "snpID_round1"], toInclude_chr8[,"snpID"], sep=";")

write.csv(oldPeakPerChr, file="/dcl01/mathias1/data/telomere_mtaub/gwas/results/ALL_GWASResults/conditional/allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop_p_lt_5e-9_peakSNPsByChromosome_round2PostBRAVO.csv", row.names=FALSE, quote=TRUE)


## need to pull new third round SNPs for chr4 and chr7
## for both, want to take the second best one based on Rasika's BRAVO assessment

oldPeakPerChrRound3<-read.csv("/dcl01/mathias1/data/telomere_mtaub/gwas/results/ALL_GWASResults/conditional/allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop_p_lt_5e-9_peakSNPsByChromosome_round3.csv", header=TRUE, stringsAsFactors = FALSE)

cond1<-fread(paste0("gunzip -c ", pathToResults, "/conditional/chr4_telomere_adjagesexstudyseqctrbatchPCs_minDP0_cond_4-163144568-A-G_4-163155406-G-A.csv.gz"), header=TRUE, sep=",")
cond1Sub<-cond1 %>% filter(!snpID %in% c(snpIDsToDrop, snpIDsToDrop_lowDepth), Score.pval < 0.1)
cond1Sub %>% arrange(Score.pval) %>% select(snpID, chr, pos, Score.pval, freq, MAC) %>% head()
toInclude<-(cond1Sub %>% arrange(Score.pval) %>% select(snpID, chr, pos, Score.pval, freq, MAC) %>% head())[1,]

cond1_chr7<-fread(paste0("gunzip -c ", pathToResults, "/conditional/chr7_telomere_adjagesexstudyseqctrbatchPCs_minDP0_cond_7-124812616-C-T_7-129041243-T-C.csv.gz"), header=TRUE, sep=",")
cond1_chr7Sub<-cond1_chr7 %>% filter(!snpID %in% c(snpIDsToDrop, snpIDsToDrop_lowDepth), Score.pval < 0.1)
cond1_chr7Sub %>% arrange(Score.pval) %>% select(snpID, chr, pos, Score.pval, freq, MAC) %>% head()
toInclude_chr7<-(cond1_chr7Sub %>% arrange(Score.pval) %>% select(snpID, chr, pos, Score.pval, freq, MAC) %>% head())[1,]

oldPeakPerChrRound3[oldPeakPerChrRound3$chr == 4, "pos_round3"]<-toInclude[,"pos"]
oldPeakPerChrRound3[oldPeakPerChrRound3$chr == 4, "snpID_round3"]<-toInclude[,"snpID"]
oldPeakPerChrRound3[oldPeakPerChrRound3$chr == 4, "Score.pval_round3"]<-toInclude[,"Score.pval"]
oldPeakPerChrRound3[oldPeakPerChrRound3$chr == 4, "toCondition_round3"]<-paste(oldPeakPerChrRound3[oldPeakPerChrRound3$chr == 4, "snpID_round1"], oldPeakPerChrRound3[oldPeakPerChrRound3$chr == 4, "snpID_round2"], toInclude[,"snpID"], sep=";")

oldPeakPerChrRound3[oldPeakPerChrRound3$chr == 7, "pos_round3"]<-toInclude_chr7[,"pos"]
oldPeakPerChrRound3[oldPeakPerChrRound3$chr == 7, "snpID_round3"]<-toInclude_chr7[,"snpID"]
oldPeakPerChrRound3[oldPeakPerChrRound3$chr == 7, "Score.pval_round3"]<-toInclude_chr7[,"Score.pval"]
oldPeakPerChrRound3[oldPeakPerChrRound3$chr == 7, "toCondition_round3"]<-paste(oldPeakPerChrRound3[oldPeakPerChrRound3$chr == 7, "snpID_round1"], oldPeakPerChrRound3[oldPeakPerChrRound3$chr == 7, "snpID_round2"],toInclude_chr7[,"snpID"], sep=";")

write.csv(oldPeakPerChrRound3, file="/dcl01/mathias1/data/telomere_mtaub/gwas/results/ALL_GWASResults/conditional/allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop_p_lt_5e-9_peakSNPsByChromosome_round3PostBRAVO.csv", row.names=FALSE, quote=TRUE)

## check chr18
cond1_chr18<-fread(paste0("gunzip -c ", pathToResults, "/conditional/chr18_telomere_adjagesexstudyseqctrbatchPCs_minDP0_cond_18-666625-C-T_18-676473-C-T.csv.gz"), header=TRUE, sep=",")
cond1_chr18Sub<-cond1_chr18 %>% filter(!snpID %in% c(snpIDsToDrop, snpIDsToDrop_lowDepth), Score.pval < 0.1)
cond1_chr18Sub %>% arrange(Score.pval) %>% select(snpID, chr, pos, Score.pval, freq, MAC) %>% head()

cond1_chr18<-fread(paste0("gunzip -c ", pathToResults, "/conditional/chr18_telomere_adjagesexstudyseqctrbatchPCs_minDP0_cond_18-666625-C-T.csv.gz"), header=TRUE, sep=",")
cond1_chr18Sub<-cond1_chr18 %>% filter(!snpID %in% c(snpIDsToDrop, snpIDsToDrop_lowDepth), Score.pval < 0.1)
cond1_chr18Sub %>% arrange(Score.pval) %>% select(snpID, chr, pos, Score.pval, freq, MAC) %>% head()


oldPeakPerChrRound4<-read.csv("/dcl01/mathias1/data/telomere_mtaub/gwas/results/ALL_GWASResults/conditional/allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop_p_lt_5e-9_peakSNPsByChromosome_round4.csv", header=TRUE, stringsAsFactors = FALSE)

cond4<-fread(paste0("gunzip -c ", pathToResults, "/conditional/chr22_telomere_adjagesexstudyseqctrbatchPCs_minDP0_cond_22-40128558-G-A_22-40023952-C-T_22-50532618-T-C.csv.gz"), header=TRUE, sep=",")
cond4Sub<-cond4 %>% filter(!snpID %in% c(snpIDsToDrop, snpIDsToDrop_lowDepth), Score.pval < 0.1)
cond4Sub %>% arrange(Score.pval) %>% select(snpID, chr, pos, Score.pval, freq, MAC) %>% head()
toInclude<-(cond4Sub %>% arrange(Score.pval) %>% select(snpID, chr, pos, Score.pval, freq, MAC) %>% head())[1,]

oldPeakPerChrRound4[oldPeakPerChrRound4$chr == 22, "pos_round4"]<-toInclude[,"pos"]
oldPeakPerChrRound4[oldPeakPerChrRound4$chr == 22, "snpID_round4"]<-toInclude[,"snpID"]
oldPeakPerChrRound4[oldPeakPerChrRound4$chr == 22, "Score.pval_round4"]<-toInclude[,"Score.pval"]
oldPeakPerChrRound4[oldPeakPerChrRound4$chr == 22, "toCondition_round4"]<-"No more rounds"

write.csv(oldPeakPerChrRound4, file="/dcl01/mathias1/data/telomere_mtaub/gwas/results/ALL_GWASResults/conditional/allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop_p_lt_5e-9_peakSNPsByChromosome_round4PostBRAVO.csv", row.names=FALSE, quote=TRUE)


oldPeakPerChrRound6<-read.csv("/dcl01/mathias1/data/telomere_mtaub/gwas/results/ALL_GWASResults/conditional/allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop_p_lt_5e-9_peakSNPsByChromosome_round6.csv", header=TRUE, stringsAsFactors = FALSE)

cond20<-fread(paste0("gunzip -c ", pathToResults, "/conditional/chr20_telomere_adjagesexstudyseqctrbatchPCs_minDP0_cond_20-63678201-C-T_20-63689775-G-A_20-36922795-A-C_20-63695521-G-A_20-63661765-C-T.csv.gz"), header=TRUE, sep=",")
cond20Sub<-cond20 %>% filter(!snpID %in% c(snpIDsToDrop, snpIDsToDrop_lowDepth), Score.pval < 0.1)
cond20Sub %>% arrange(Score.pval) %>% select(snpID, chr, pos, Score.pval, freq, MAC) %>% head()
cond20Sub %>% arrange(Score.pval) %>% select(snpID, chr, pos, Score.pval, freq, MAC) %>% filter(Score.pval < 5e-9)
toInclude20<-(cond20Sub %>% arrange(Score.pval) %>% select(snpID, chr, pos, Score.pval, freq, MAC) %>% filter(snpID == "20:63676585:C:A"))

cond10<-fread(paste0("gunzip -c ", pathToResults, "/conditional/chr10_telomere_adjagesexstudyseqctrbatchPCs_minDP0_cond_10-103916707-C-A_10-103918153-C-A_10-94344908-T-G_10-99514276-A-G_10-103915847-C-T.csv.gz"), header=TRUE, sep=",")
cond10Sub<-cond10 %>% filter(!snpID %in% c(snpIDsToDrop, snpIDsToDrop_lowDepth), Score.pval < 0.1)
cond10Sub %>% arrange(Score.pval) %>% select(snpID, chr, pos, Score.pval, freq, MAC) %>% head()
cond10Sub %>% arrange(Score.pval) %>% select(snpID, chr, pos, Score.pval, freq, MAC) %>% filter(Score.pval < 5e-9)
toInclude10<-(cond10Sub %>% arrange(Score.pval) %>% select(snpID, chr, pos, Score.pval, freq, MAC) %>% filter(snpID == "10:103907794:G:T"))

oldPeakPerChrRound6[oldPeakPerChrRound6$chr == 20, "pos_round6"]<-toInclude20[,"pos"]
oldPeakPerChrRound6[oldPeakPerChrRound6$chr == 20, "snpID_round6"]<-toInclude20[,"snpID"]
oldPeakPerChrRound6[oldPeakPerChrRound6$chr == 20, "Score.pval_round6"]<-toInclude20[,"Score.pval"]
oldPeakPerChrRound6[oldPeakPerChrRound6$chr == 20, "toCondition_round6"]<-paste(oldPeakPerChrRound6[oldPeakPerChrRound6$chr == 20, c("snpID_round1", "snpID_round2","snpID_round3", "snpID_round4", "snpID_round5", "snpID_round6")], collapse=";")

oldPeakPerChrRound6[oldPeakPerChrRound6$chr == 10, "pos_round6"]<-toInclude10[,"pos"]
oldPeakPerChrRound6[oldPeakPerChrRound6$chr == 10, "snpID_round6"]<-toInclude10[,"snpID"]
oldPeakPerChrRound6[oldPeakPerChrRound6$chr == 10, "Score.pval_round6"]<-toInclude10[,"Score.pval"]
oldPeakPerChrRound6[oldPeakPerChrRound6$chr == 10, "toCondition_round6"]<-paste(oldPeakPerChrRound6[oldPeakPerChrRound6$chr == 10, c("snpID_round1", "snpID_round2", "snpID_round3", "snpID_round4", "snpID_round5", "snpID_round6")], collapse=";")

write.csv(oldPeakPerChrRound6, file="/dcl01/mathias1/data/telomere_mtaub/gwas/results/ALL_GWASResults/conditional/allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop_p_lt_5e-9_peakSNPsByChromosome_round6PostBRAVO.csv", row.names=FALSE, quote=TRUE)

## round 7 

oldPeakPerChrRound7<-read.csv("/dcl01/mathias1/data/telomere_mtaub/gwas/results/ALL_GWASResults/conditional/allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop_p_lt_5e-9_peakSNPsByChromosome_round7.csv", header=TRUE, stringsAsFactors = FALSE)

cond20<-fread(paste0("gunzip -c ", pathToResults, "/conditional/chr20_telomere_adjagesexstudyseqctrbatchPCs_minDP0_cond_20-63678201-C-T_20-63689775-G-A_20-36922795-A-C_20-63695521-G-A_20-63661765-C-T_20-63676585-C-A.csv.gz"), header=TRUE, sep=",")
cond20Sub<-cond20 %>% filter(!snpID %in% c(snpIDsToDrop, snpIDsToDrop_lowDepth), Score.pval < 0.1)
cond20Sub %>% arrange(Score.pval) %>% select(snpID, chr, pos, Score.pval, freq, MAC) %>% filter(Score.pval < 5e-9)

cond5<-fread(paste0("gunzip -c ", pathToResults, "/conditional/chr5_telomere_adjagesexstudyseqctrbatchPCs_minDP0_cond_5-1285859-C-A_5-1292331-C-T_5-1287079-G-A_5-1272383-C-T_5-1292843-C-T_5-1280823-A-G.csv.gz"), header=TRUE, sep=",")
cond5Sub<-cond5 %>% filter(!snpID %in% c(snpIDsToDrop, snpIDsToDrop_lowDepth), Score.pval < 0.1)
cond5Sub %>% arrange(Score.pval) %>% select(snpID, chr, pos, Score.pval, freq, MAC) %>% filter(Score.pval < 5e-9)

## typo from above
oldPeakPerChrRound7[oldPeakPerChrRound7$chr == 10, "toCondition_round7"]<-"No more rounds"
## none on chr20 pass 
oldPeakPerChrRound7[oldPeakPerChrRound7$chr == 20, "toCondition_round7"]<-"No more rounds"

write.csv(oldPeakPerChrRound7, file="/dcl01/mathias1/data/telomere_mtaub/gwas/results/ALL_GWASResults/conditional/allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop_p_lt_5e-9_peakSNPsByChromosome_round7PostBRAVO.csv", row.names=FALSE, quote=TRUE)

## check hopefully final round
cond5<-fread(paste0("gunzip -c ", pathToResults, "/conditional/chr5_telomere_adjagesexstudyseqctrbatchPCs_minDP0_cond_5-1285859-C-A_5-1292331-C-T_5-1287079-G-A_5-1272383-C-T_5-1292843-C-T_5-1280823-A-G_5-139637905-T-A.csv.gz"), header=TRUE, sep=",")
cond5Sub<-cond5 %>% filter(!snpID %in% c(snpIDsToDrop, snpIDsToDrop_lowDepth), Score.pval < 0.1)
cond5Sub %>% arrange(Score.pval) %>% select(snpID, chr, pos, Score.pval, freq, MAC) %>% filter(Score.pval < 5e-9)

cond22<-fread("gunzip -c /dcl01/mathias1/data/telomere_mtaub/gwas/results/ALL_GWASResults/conditional/chr22_telomere_adjagesexstudyseqctrbatchPCs_minDP0_cond_22-50532618-T-C.csv.gz", header= TRUE, sep=",")
cond22Sub<-cond22 %>% filter(!snpID %in% c(snpIDsToDrop, snpIDsToDrop_lowDepth), Score.pval < 0.1)
cond22Sub %>% arrange(Score.pval) %>% select(snpID, chr, pos, Score.pval, freq, MAC) %>% filter(Score.pval < 5e-9)

## OLD VERSION OF FUNCTION (kind of)
makePlotAndPullNext <- function(condSNP){
  
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
  
  ## for now this will work, but will need to update once there is more than one file per chr
  ## read in appropriate conditional data files
  for (i in 1:nSNPs){
    currCondRes<-fread(paste0("gunzip -c ", pathToResults, "/conditional/", chr, "_telomere_adjagesexstudyseqctrbatchPCs_minDP0_cond_", paste(SNPTags[1:i], collapse="_"), ".csv.gz"), header=TRUE, sep=",")
    currCondResSub<-currCondRes %>% filter(!snpID %in% c(snpIDsToDrop, snpIDsToDrop_lowDepth), Score.pval < 0.1, !snpID %in% condSNPVec[1:i])
    
  }
  cond1<-fread(paste0("gunzip -c ", pathToResults, "/conditional/", chr, "_*.csv.gz"), header=TRUE, sep=",")
  cond1Sub<-cond1 %>% filter(!snpID %in% c(snpIDsToDrop, snpIDsToDrop_lowDepth), Score.pval < 0.1, !snpID %in% condSNP)
  
  yLimTop<-ceiling(max(max(-log10(nonCondSub$Score.pval)), max(-log10(cond1Sub$Score.pval))))
  
  dtrack_nonCond <- DataTrack(data=-log10(nonCondSub$Score.pval), start=nonCondSub$pos-1, end= nonCondSub$pos, genome="hg38", chromosome=chr, name= paste0(chr, " Non conditional"), type=c("p","g"), col="black", baseline = -log10(5e-9), col.baseline="blue", ylim=c(0, yLimTop))
  
  
  dtrack_cond1 <- DataTrack(data=-log10(cond1Sub$Score.pval), start=cond1Sub$pos-1, end= cond1Sub$pos, genome="hg38", chromosome=chr, name= paste0(chr, " Conditional on ", condSNP), type=c("p","g"), col="black", baseline = -log10(5e-9), col.baseline="blue", ylim=c(0, yLimTop))
  
  ht <- HighlightTrack(trackList=list(dtrack_nonCond, dtrack_cond1, gtrack), start=c(condSNPPos-10000), width=20000, chromosome=chr)
  
  #Create PDF
  png(paste0("/dcl01/mathias1/data/telomere_mtaub/gwas/results/ALL_GWASResults/conditional/", chr, "_cond_", gsub("\\:", "_", condSNP), ".png"),height=1200,width=1600)
  #Plot Data
  plotTracks(list(ht),background.title="darkgray")
  #End
  dev.off()
  
  newPeak <- cond1Sub %>% filter(Score.pval == min(Score.pval)) %>% select(snpID, chr, pos, Score.pval)
  return(newPeak)
}

