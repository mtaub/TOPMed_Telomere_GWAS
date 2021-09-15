## want to start with checking peak SNPs in loci IDed by Rasika using previous version of results

library(readxl)
library(GenomicRanges)
library(data.table)
library(dplyr)



sigSNPs<-read.csv(gzfile("~/Research/OneDrive/telomere/gwas/results/ALL_GWASResults/allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop_p_lt_5e-8.csv.gz"))

sigSNPs$diffPos <- c(0,diff(sigSNPs$pos))
sigSNPs$sameChr <- c(TRUE, sigSNPs$chr[1:(nrow(sigSNPs)-1)] == sigSNPs$chr[2:(nrow(sigSNPs))] )
contigBreaks<-which(!sigSNPs$sameChr  | sigSNPs$diffPos > 200000)
sigSNPs[unique(sort(c(contigBreaks-1, contigBreaks))),1:5]

sigSNPs$locusGroup<-rep(1:(length(contigBreaks)+1), times=c(contigBreaks[1]-1, diff(contigBreaks), nrow(sigSNPs)-contigBreaks[length(contigBreaks)]+1))
groupBreaks<-which(diff(sigSNPs$locusGroup)>0)
sigSNPs[unique(sort(c(groupBreaks, groupBreaks+1))),1:5]

newPeaks<-sigSNPs %>% group_by(locusGroup) %>% filter(Score.pval == min(Score.pval))

write.csv(newPeaks, file="~/Research/OneDrive/telomere/gwas/results/ALL_GWASResults/allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop_p_lt_5e-8_peakSNPsByLocus.csv", row.names=FALSE, quote=FALSE, na="")


## based on Rasika's examination, I need to merge certain groups
## 5-10 should be merged
## 18 and 19
## 25 and 26
## 34 and 35
sigSNPs$locusGroup[sigSNPs$locusGroup %in% c(5:10)]<-"5Merge"
sigSNPs$locusGroup[sigSNPs$locusGroup %in% c("18", "19")]<-"18Merge"
sigSNPs$locusGroup[sigSNPs$locusGroup %in% c("25", "26")]<-"25Merge"
sigSNPs$locusGroup[sigSNPs$locusGroup %in% c("34", "35")]<-"34Merge"


newPeaks<-sigSNPs %>% group_by(locusGroup) %>% filter(Score.pval == min(Score.pval))

newPeaks %>% select(chr, pos, snpID) %>% print(n=Inf)

## pull locus names from some file -- here using most recent LZ plot file that Ras put together for Matt
peakNames<-read.csv("~/Research/OneDrive/telomere/gwas/results/ALL_GWASResults/OASIS QC/Telomere_LZdata_040220_toMatt.csv", header=TRUE, stringsAsFactors = FALSE)
peakNames <- peakNames %>% mutate(pos=Pos38_LZ_Start_500kb + 500000) %>% select(newLOCUS, chr, pos)

## merge in locus names; a couple places Ras had not picked peak SNP when she collapsed
newPeaks <- newPeaks %>% left_join(peakNames)
newPeaks[newPeaks$locusGroup == "5Merge", "newLOCUS"]<-"LINC00901"
newPeaks[newPeaks$locusGroup == "25Merge", "newLOCUS"]<-"LOC101928283"

sigSNPs <- newPeaks %>% select(locusGroup, newLOCUS) %>% right_join(sigSNPs)

write.csv(newPeaks, file="~/Research/OneDrive/telomere/gwas/results/ALL_GWASResults/allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop_p_lt_5e-8_peakSNPsByLocusMerged.csv", row.names=FALSE, quote=FALSE, na="")
write.csv(sigSNPs, file="~/Research/OneDrive/telomere/gwas/results/ALL_GWASResults/allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop_p_lt_5e-8_withLocusGroups.csv", row.names=FALSE, quote=FALSE, na="")

## pull peak SNPs for conditional analysis

peakPerChr<-newPeaks %>% group_by(chr) %>% filter(Score.pval == min(Score.pval)) %>% select(snpID, chr, pos, Score.pval) %>% filter(Score.pval < 5e-9)
write.csv(peakPerChr, file="~/Research/OneDrive/telomere/gwas/results/ALL_GWASResults/allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop_p_lt_5e-9_peakSNPsByChromosome.csv", row.names=FALSE, quote=FALSE, na="")

########################
### OLD CODE -- in case I need for splitting things up again
########################


toCheck<-read.csv("~/Research/OneDrive/telomere/gwas/results/ALL_GWASResults/OASIS QC/Telomere_LZplots_032720_toJim.csv", header=TRUE, stringsAsFactors = FALSE)
newPeaks %>% filter(!snpID %in% toCheck$VARID) %>% select(snpID, chr, pos, Score.pval, filt, low_depth_rate)
toCheck %>% filter(!VARID %in% newPeaks$snpID)

toCheck2<-read.csv("~/Research/OneDrive/telomere/gwas/results/ALL_GWASResults/OLD/allChrs_telomere_adjagesexstudyseqctrbatchPCs_BRAVODrop_p_lt_5e-8_peakSNPsByLocus.csv", header=TRUE, stringsAsFactors = FALSE)
newPeaks %>% filter(!snpID %in% toCheck2$snpID)

newPeaks %>% filter(chr=="3") %>% select(snpID, chr, pos, Score.pval, filt, low_depth_rate)
newPeaks %>% filter(chr=="4") %>% select(snpID, chr, pos, Score.pval, filt, low_depth_rate)
newPeaks %>% filter(chr=="7") %>% select(snpID, chr, pos, Score.pval, filt, low_depth_rate)
newPeaks %>% filter(chr=="X") %>% select(snpID, chr, pos, Score.pval, filt, low_depth_rate)

## in rasika's file but not ours
toCheck %>% filter(!VARID %in% newPeaks$snpID) %>% select(VARID, chr, Pos38, `X.`)
#SNPname         LOCUS        Chr38     Pos38
#<chr>           <chr>        <chr>     <dbl>
#  1 7:124802890:C:T POT1         7     124802890 -- replaced by 7:124812616:C:T
#2 18:660441:GA:G  TYMS         18       660441 -- replaced by 18:666625:C:T
#3 16:69941728:A:T WWP2         16     69941728 -- got included with peak represented by 16:69357811:C:G (TERF2)
#4 9:34116085:T:C  DCAF12       9      34116085 -- replaced by 9:34130437:G:A
#5 10:94264600:C:T PLCE1        10     94264600 -- got included wiht peak represented by 10:94344908:T:G (NOC3L)
#6 16:74226723:A:T LOC101928035 16     74226723 -- no longer significant
#7 10:94700466:G:T CYP2C18      10     94700466 -- no longer significant
#8 6:139024386:A:T ABRACL       6     139024386 -- no longer significant


## in our file but not Rasika's
newPeaks %>% filter(!snpID %in% toCheck$SNPname) %>% select(snpID, chr, pos)
#locusGroup snpID            chr         pos
#<int> <fct>            <fct>     <int>
#  1         19 6:28877976:AAC:A 6      28877976
# 2         24 6:153531000:T:A  6     153531000
#3         25 7:124812616:C:T  7     124812616 -- replaces 7:124802890:C:T for POT1
#4         32 9:34130437:G:A   9      34130437 -- replace 9:34116085:T:C for DCAF12
#5         42 13:72766039:G:A  13     72766039
#6         45 14:91415360:C:T  14     91415360
#7         54 18:666625:C:T    18       666625 -- replaces 18:660441:GA:G for TYMS


toCheck %>% filter(Chr38 == "16") %>% select(SNPname, LOCUS, Chr38, Pos38)
newPeaks %>% filter(chr=="16") %>% select(snpID, chr, pos)

toCheck %>% filter(Chr38 == "10") %>% select(SNPname, LOCUS, Chr38, Pos38)
newPeaks %>% filter(chr=="10") %>% select(snpID, chr, pos)

toCheck %>% filter(Chr38 == "6") %>% select(SNPname, LOCUS, Chr38, Pos38)
newPeaks %>% filter(chr=="6") %>% select(snpID, chr, pos)

## include locus name from Rasika's file where possible
newPeaks <- toCheck %>% select(snpID = SNPname, LOCUS) %>% right_join(newPeaks)
newPeaks$newLOCUS <- newPeaks$LOCUS
newPeaks <- newPeaks %>% select(snpID, LOCUS, newLOCUS, variant.id:locusGroup)
newPeaks[newPeaks$snpID == "7:124812616:C:T", "newLOCUS"] <- "POT1"
newPeaks[newPeaks$snpID == "9:34130437:G:A", "newLOCUS"] <- "DCAF12"
newPeaks[newPeaks$snpID == "18:666625:C:T", "newLOCUS"] <- "TYMS"

sigSNPs <- newPeaks %>% select(peakPos=pos, locusGroup) %>% right_join(sigSNPs, by="locusGroup")
sigSNPs <- sigSNPs %>% mutate(distToPeak = pos - peakPos)
sigSNPs <- sigSNPs %>% select(variant.id:sameChr,peakPos, distToPeak,locusGroup)
sigSNPs %>% group_by(locusGroup) %>% summarize(maxDist = max(abs(distToPeak))) %>% filter(maxDist > 200000)





#toCheck<-read_excel("~/Research/telomere/newManuscript/results/ALL_GWASResults/Chromosome_BRAVO_confirm_loci.xlsx")
#toCheck <- toCheck %>% mutate(Chr38 = gsub("\\:.*", "", SNPname), Pos38=as.numeric(sapply(strsplit(SNPname, split=":"), function(x) x[[2]])))



## split up all loci with distances to peak >200kb
sigSNPs<- sigSNPs %>% mutate(locusGroupSplit = ifelse(distToPeak < -200000, paste0(locusGroup, "Pre"), ifelse(distToPeak > 200000, paste0(locusGroup, "Post"), locusGroup))) 
refinedPeaks<-sigSNPs %>% group_by(locusGroupSplit) %>% filter(Score.pval == min(Score.pval))

sigSNPs <- refinedPeaks %>% select(refinedPeakPos=pos, locusGroupSplit) %>% right_join(sigSNPs, by="locusGroupSplit")
sigSNPs <- sigSNPs %>% mutate(distToRefinedPeak = pos - refinedPeakPos)
sigSNPs <- sigSNPs %>% select(variant.id:sameChr,peakPos, distToPeak,locusGroup, refinedPeakPos, distToRefinedPeak,locusGroupSplit)
sigSNPs %>% group_by(locusGroupSplit) %>% summarize(maxDistRefined = max(abs(distToRefinedPeak))) %>% filter(maxDistRefined > 200000)

## split one more time for last remaining locus with distances > 200kb
sigSNPs<- sigSNPs %>% mutate(locusGroupSplit2 = ifelse(distToRefinedPeak < -200000, paste0(locusGroupSplit, "Pre"), ifelse(distToRefinedPeak > 200000, paste0(locusGroupSplit, "Post"), locusGroupSplit))) 
refinedPeaks2<-sigSNPs %>% group_by(locusGroupSplit2) %>% filter(Score.pval == min(Score.pval))

sigSNPs <- refinedPeaks2 %>% select(refinedPeakPos2=pos, locusGroupSplit2) %>% right_join(sigSNPs, by="locusGroupSplit2")
sigSNPs <- sigSNPs %>% mutate(distToRefinedPeak2 = pos - refinedPeakPos2)
sigSNPs <- sigSNPs %>% select(variant.id:sameChr,peakPos, distToPeak,locusGroup, refinedPeakPos, distToRefinedPeak,locusGroupSplit, refinedPeakPos2, distToRefinedPeak2,locusGroupSplit2)
sigSNPs %>% group_by(locusGroupSplit2) %>% summarize(maxDistRefined2 = max(abs(distToRefinedPeak2))) %>% filter(maxDistRefined2 > 200000)


refinedPeaks2 <- toCheck %>% select(snpID = SNPname, LOCUS) %>% right_join(refinedPeaks2)
refinedPeaks2$newLOCUS <- refinedPeaks2$LOCUS
refinedPeaks2 <- refinedPeaks2 %>% select(snpID, LOCUS, newLOCUS, variant.id:locusGroupSplit2)
refinedPeaks2[refinedPeaks2$snpID == "7:124812616:C:T", "newLOCUS"] <- "POT1"
refinedPeaks2[refinedPeaks2$snpID == "9:34130437:G:A", "newLOCUS"] <- "DCAF12"
refinedPeaks2[refinedPeaks2$snpID == "18:666625:C:T", "newLOCUS"] <- "TYMS"
refinedPeaks2[refinedPeaks2$snpID == "16:69972064:A:G", "newLOCUS"] <- "WWP2"
