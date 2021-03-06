################################################
#### Code for creating supplementary figure with three parts:
### A. Computel vs TelSeq results on JHS data
### B. Computel vs TelSeq computing time boxplot on JHS data
### C. Computel and TelSeq vs Southern Blot on JHS data
### D. Batch-corrected TelSeq vs Southern Blot on JHS data
### E. Batch-corrected TelSeq vs FlowFISH on GeneSTAR data
################################################

library(ggplot2)
library(tidyverse)
library(readxl)
library(gridExtra)
library(ggpubr)


## need data set with uncorrected length estimates
## not sure I am actually going to use this here though...
load("~/Research/OneDrive/telomere/gwas/results/allResMerge_allSamples_PCAdj.rda")
load("~/Research/OneDrive/telomere/gwas/results/allResMerge_forAnalysis_031520.rda")
forAnalysisNoAdj <- allResMerge %>% select(NWDID, LENGTH_ESTIMATE) %>% right_join(forAnalysis)
#save(forAnalysisNoAdj, file = "~/Research/OneDrive/telomere/gwas/results/allResMerge_forAnalysis_031520_wUnadjTL.rda")




## Computel results -- generated by UM on GeneSTAR and JHS samples
load("~/Research/OneDrive/telomere/gwas/inputFiles/allResComp.rda")
allResComp<-data.frame(NWDID=names(allResComp), Computel=allResComp)
dim(allResComp) ## n = 5257 (includes GeneSTAR and JHS)
allResComp <- allResComp %>% mutate(Computel=Computel/1000)



## Southern blot data with accession numbers (N=2470)
ltpJHS<-read.csv("~/Research/OneDrive/telomere/gwas/inputFiles/mean_TRF_dbgap.csv")
ltpJHS<-ltpJHS %>% select(Accession = SUBJECT_ID, Southern=Mean_TRF)
dim(ltpJHS) ## 2470

## age and sex with accession numbers (may be more samples than in phenoAll)
covars<-read.table("~/Research/OneDrive/telomere/gwas/inputFiles/JHS_telomere_data_09142017.txt", header=TRUE, stringsAsFactors=FALSE, sep="\t") %>% select(Accession = SUBJECT_ID, wgs_age = wgs_age, telomere_age = telomere_age, male = male)
dim(covars)

ltpJHS<-full_join(covars, ltpJHS)
dim(ltpJHS)

## need this to map id from file with age to telSeq data (only 3406 samples)
# downloaded from dbGaP PhenoGenotypeFiles/RootStudyConsentSet_phs000964.TOPMed_WGS_JHS.v2.p1.c1.HMB-IRB-NPU/PhenotypeFiles/
mapJHS<-read.table("~/Research/OneDrive/telomere/gwas/inputFiles/phs000964.v2.pht004839.v2.p1.TOPMed_WGS_JHS_Sample.MULTI.txt", comment.char="#", header=TRUE, stringsAsFactors = FALSE) %>% select(Accession = Accession, NWDID = SUBJECT_ID)
dim(mapJHS)

ltpJHS<-ltpJHS %>% full_join(mapJHS)
dim(ltpJHS)
## remove two duplicates
ltpJHS <- ltpJHS %>% filter(!duplicated(NWDID))
dim(ltpJHS)
sum(!is.na(ltpJHS$Southern))

## all the samples that I have plate information for
JHSBatch<-read_excel("~/Research/OneDrive/telomere/gwas/inputFiles/JHS_sample_prep_information.xls") %>% select(NWDID = nwd_ids, Plate=plate_ids, Flowcells = Flowecells, DmuxDate = 'Dmux date', Lane = 'Lane #') %>% select(NWDID, Plate) %>% filter(!is.na(NWDID))
dim(JHSBatch) ## n=3418


## want to maximize samples that have TelSeq, Batch, Computel, Southern blot, i.e., don't filter if possible
forAnalysisNoAdj<-forAnalysisNoAdj %>% left_join(JHSBatch)
forAnalysisNoAdj<-forAnalysisNoAdj %>% left_join(allResComp)
forAnalysisNoAdj<-forAnalysisNoAdj %>% left_join(ltpJHS)


## this contains way more data than I need, but want to be sure I have it all
forPlots<-forAnalysisNoAdj %>% select(NWDID, TelSeq=LENGTH_ESTIMATE, TelSeqAdj = LENGTH_ESTIMATE_BATCHADJ, Computel, Southern, Plate, study, wgs_age, telomere_age, sex, seq_center_new) %>% mutate(TelSeqAdjRecenter = TelSeqAdj + 3.311832)
dim(forPlots)

forSBPlots<-forPlots %>% filter(!is.na(Computel), !is.na(Southern), !is.na(Plate))
dim(forSBPlots) ## 2398 


allJHSLong2<-forSBPlots %>% gather(key="Method", value="Estimate", TelSeq, Computel) 
allJHSLong2$Method<-factor(allJHSLong2$Method, levels=c("TelSeq", "Computel"))
oindTel <- order(tapply(forSBPlots$TelSeq, forSBPlots$Plate, median, na.rm=TRUE))
allJHSLong2$Plate <- factor(allJHSLong2$Plate, levels=levels(factor(forSBPlots$Plate))[oindTel])

getCorText<-function(x,y){
  r <- cor(x, y, use="complete.obs")
  txt <- format(c(r, 0.123456789), digits = 2)[1]
  txt <- paste("r= ", txt, sep = "")
  return(txt)
}

corrJHS = allJHSLong2 %>% 
  group_by(Method) %>% 
  summarize(cor = getCorText(Estimate, Southern))

plotLabs<-data.frame(labText=c(" ", "c"), x=3.7, y=4.8, Method=c("TelSeq", "Computel"), plate_ids=NA)

SBPlot<-allJHSLong2 %>% 
  ggplot(aes(x = Estimate, y = Southern, col=Plate)) + 
  geom_point() + 
  facet_wrap( ~ Method, scales = "free") +
  theme(legend.position="none") +
  geom_text(aes(x=1.3, y=10, label = cor), 
            data = filter(corrJHS, Method == "TelSeq"), colour = "black") +
  geom_text(aes(x=1.3, y=10, label = cor), 
            data = filter(corrJHS, Method == "Computel"), colour = "black") +
  xlab("WGS Estimate (kb)") +
  ylab("Southern Blot (kb)") +
  geom_text(data=plotLabs, mapping=aes(x=x, y=y, label=labText), col="black", cex=6) +
  theme_bw() +
  theme(legend.position="none")
print(SBPlot)

## add plot with batch adjustment
SBPlotAdj<-forSBPlots %>% 
  ggplot(aes(x = TelSeqAdjRecenter, y = Southern, col=Plate)) + 
  geom_point() + 
  theme(legend.position="none") +
  xlab("bPC-adjusted TelSeq (kb)") +
  ylab("Southern Blot (kb)") +
  annotate("text", x = 2.1, y = 10, label = getCorText(forSBPlots$TelSeqAdjRecenter, forSBPlots$Southern)) +
  annotate("text", x=5.4, y=4.95, label="d", cex=6) +
  theme_bw() +
  theme(legend.position="none")
print(SBPlotAdj)

# models comparing percent variablity explained
summary(lm(TelSeq ~ wgs_age, data = forSBPlots))
summary(lm(TelSeqAdj ~ wgs_age, data = forSBPlots))

summary(lm(TelSeq ~ Plate, data = forSBPlots))
summary(lm(TelSeqAdj ~ Plate, data = forSBPlots))

# this time around, just sticking with the same samples that I have SB data for
#forCompPlots<-forPlots %>% filter(!is.na(Computel), !is.na(Plate))
#dim(forCompPlots) ## 3362
forCompPlots <- forSBPlots

compPlot<-forCompPlots %>%
  ggplot(aes(x = TelSeq, y = Computel, col=Plate)) + 
  geom_point() + 
  theme(legend.position="none") + 
  xlab("TelSeq (kb)") + 
  ylab("Computel (kb)") +
  annotate("text", x = 1.0, y = 4, label = getCorText(forCompPlots$TelSeq, forCompPlots$Computel)) +
  annotate("text", x=4.2, y=0.5, label="a", cex=6) +
  theme_bw() + 
  theme(legend.position="none")
print(compPlot)

## want to check on compute times for telseq and computel

telSeqDir<-"~/Research/OneDrive/telomere/gwas/inputFiles/call-run_telseq"
allDirsTel<-system(paste0("ls ", telSeqDir), intern=TRUE)

telTimes<-sapply(allDirsTel, function(x){
  if (length(system(paste0("ls ", telSeqDir, "/", x, "/*.telseq.out"), intern=TRUE))>0){
    basename<-gsub(".*\\/", "", gsub(".telseq.out", "", system(paste0("ls ", telSeqDir, "/", x, "/*.telseq.out"), intern=TRUE)))
    shard<-gsub(".*\\-", "", x)
    res<-system(paste0("tail -1 ", telSeqDir, "/", x, "/run_telseq-", shard, "-stderr.log"), intern=TRUE)
    res<-as.numeric(gsub("s", "", gsub(".*\\ ", "", res)))
    names(res)<-basename
    return(res)
  }
  maxAttempt<-max(as.numeric(gsub(".*\\-", "", system(paste0("ls -d ", telSeqDir, "/", x, "/attempt*"), intern=TRUE))))
  if (length(system(paste0("ls ", telSeqDir, "/", x, "/attempt-", maxAttempt, "/*.telseq.out"), intern=TRUE))>0){
    basename<-gsub(".*\\/", "", gsub(".telseq.out", "", system(paste0("ls ", telSeqDir, "/", x, "/attempt-", maxAttempt, "/*.telseq.out"), intern=TRUE)))
    shard<-gsub(".*\\-", "", x)
    res<-system(paste0("tail -1 ", telSeqDir, "/", x, "/attempt-", maxAttempt, "/run_telseq-", shard, "-stderr.log"), intern=TRUE)
    res<-as.numeric(gsub("s", "", gsub(".*\\ ", "", res)))
    names(res)<-basename
    return(res)
  }
  return("No results")
})

names(telTimes)<-sub(".*\\.", "", names(telTimes))
telTimesClean<-telTimes[!telTimes < 200]
telTimesCleanDF<-data.frame(NWDID=names(telTimesClean), TelSeq=telTimesClean)

compDir<-"~/Research/OneDrive/telomere/gwas/inputFiles/call-run_computel"
allDirsComp<-system("ls ~/Research/OneDrive/telomere/gwas/inputFiles/call-run_computel", intern=TRUE)

compTimes<-sapply(allDirsComp, function(x){
  if (length(system(paste0("ls ", compDir, "/", x, "/*.computel.tsv"), intern=TRUE))>0){
    basename<-gsub(".*\\/", "", gsub(".computel.tsv", "", system(paste0("ls ", compDir, "/", x, "/*.computel.tsv"), intern=TRUE)))
    shard<-gsub(".*\\-", "", x)
    res1<-system(paste0("head -2 ", compDir, "/", x, "/run_computel-", shard, "-stdout.log | tail -1"), intern=TRUE)
    res1<-strptime(res1,  "%a %b %d %H:%M:%S UTC %Y")
    res2<-system(paste0("tail -1 ", compDir, "/", x, "/run_computel-", shard, "-stdout.log"), intern=TRUE)
    res2<-strptime(res2,  "%a %b %d %H:%M:%S UTC %Y")
    res<-as.numeric(difftime(res2, res1, units="secs"))
    names(res)<-basename
    return(res)
  }
  maxAttempt<-max(as.numeric(gsub(".*\\-", "", system(paste0("ls -d ", compDir, "/", x, "/attempt*"), intern=TRUE))))
  if (length(system(paste0("ls ", compDir, "/", x, "/attempt-", maxAttempt, "/*.computel.tsv"), intern=TRUE))>0){
    basename<-gsub(".*\\/", "", gsub(".computel.tsv", "", system(paste0("ls ", compDir, "/", x, "/attempt-", maxAttempt, "/*.computel.tsv"), intern=TRUE)))
    shard<-gsub(".*\\-", "", x)
    res1<-system(paste0("head -2 ", compDir, "/", x, "/attempt-", maxAttempt,"/run_computel-", shard, "-stdout.log | tail -1"), intern=TRUE)
    res1<-strptime(res1,  "%a %b %d %H:%M:%S UTC %Y")
    res2<-system(paste0("tail -1 ", compDir, "/", x, "/attempt-", maxAttempt,"/run_computel-", shard, "-stdout.log"), intern=TRUE)
    res2<-strptime(res2,  "%a %b %d %H:%M:%S UTC %Y")
    res<-as.numeric(difftime(res2, res1, units="secs"))
    names(res)<-basename
    return(res)
  }
  return("No results")
})

names(compTimes)<-sub(".*\\.", "", names(compTimes))
compTimesDF<-data.frame(NWDID=names(compTimes), Computel=compTimes)

bothTimesDF<-forCompPlots %>% select(NWDID) %>% left_join(telTimesCleanDF)
bothTimesDF<-bothTimesDF %>% left_join(compTimesDF)

all(bothTimesDF$NWDID %in% forSBPlots$NWDID)

bothTimesDFLong<-bothTimesDF %>% gather(key="Method", value="Time", TelSeq, Computel)
bothTimesDFLong$Method<-factor(bothTimesDFLong$Method, levels=c("TelSeq", "Computel"))


timePlot<-bothTimesDFLong %>% 
  ggplot(aes(x = Method, y = Time/3600)) + 
  geom_boxplot() + 
  ylab("Time (hours)")+
  annotate("text", x=2.5, y=0.7, label="b", cex=6) +
  #ylim(c(0.5, 21))+
  theme_bw()
print(timePlot)

## flow fish plot

## need to use TelSeq results from first round of UM running because it is appropriate for the shorter read lenghts
## get telRes19 object from getResFromUM.R (run lines 1-74 in getResFromUM.R)

telRes19Long<-telRes19 %>% dplyr::select(NWD_ID, TelSeq=UMTel, Computel=UMComp, TL_Lymphs, TL_Grans) %>% gather(key="Method", value="Estimate", TelSeq, Computel) %>% gather(key="CellType", value="FlowFish", TL_Lymphs, TL_Grans)
telRes19Long$Method<-factor(telRes19Long$Method, levels=c("TelSeq", "Computel"))

getCorText<-function(x,y){
  r <- cor(x, y, use="complete.obs")
  txt <- format(c(r, 0.123456789), digits = 2)[1]
  txt <- paste("r= ", txt, sep = "")
  return(txt)
}

corr2 = telRes19Long %>% 
  group_by(Method, CellType) %>% 
  summarize(cor = getCorText(Estimate, FlowFish))




ffPlot<-telRes19Long %>% filter(Method == "TelSeq") %>% 
  ggplot(aes(x = Estimate, y = FlowFish, col=CellType)) + 
  geom_point() + 
  geom_text(aes(x = 5.5, y = 7, label = cor), 
            data = filter(corr2, CellType=="TL_Grans", Method=="TelSeq"), colour = "blue") +
  geom_text(aes(x = 5.5 , y = 6.85, label = cor), 
            data = filter(corr2, CellType=="TL_Lymphs", Method=="TelSeq"), colour = "orange") +
  scale_colour_manual(values = c("blue", "orange"), labels=c(TL_Grans = "Gran", TL_Lymphs="Lymph")) +
  theme(legend.justification=c(0,1), legend.position=c(0,1))+
  xlab("TelSeq (kb)") +
  ylab("FlowFISH (kb)") +
  geom_text(aes(x=7.4, y=3.65, label="e"), col="black", cex=6) + 
  theme_bw()
print(ffPlot)

pdf("~/Research/OneDrive/telomere/gwas/manuscript/figPDFs/SuppFig1_ComputelTelSeq.pdf", height=8, width=12)
grid.arrange(compPlot, timePlot, SBPlot, SBPlotAdj, ffPlot, nrow=2)
dev.off()

############################
####### Suppl Fig 2
############################

novelLociRep<-read.csv("~/Research/OneDrive/telomere/gwas/manuscript/tableCSVs/SupplTable5_replication_novel_TOPMedResults.csv", header= TRUE, stringsAsFactors = FALSE, na.strings = "-")
novelLociRepLong <- novelLociRep %>% mutate(isSig = !threshold == ">0.05") %>% select(TOPMed = Est_joint, Li = li_BETA, Dorajoo = dorajoo_beta, isSig) %>% pivot_longer(cols = Li:Dorajoo, names_to = "study", values_to = "estimate") %>% filter(!is.na(estimate))


getCorPvalText <- function(x,y){
  r <- cor.test(x, y, use="complete.obs")
  rtxt <- format(c(r$estimate, 0.123456789), digits = 2)[1]
  ptxt <- format(c(r$p.value, 0.123456789), digits = 2)[1]
  txt <- paste("r = ", rtxt, "\n p = ", ptxt, sep = "")
  return(txt)
}
corrNovelAll = novelLociRepLong %>% 
  group_by(study) %>% 
  summarize(cor = getCorPvalText(estimate, TOPMed)) 
corrNovelAll$lab <- c("a", "b")

corrNovelSig = novelLociRepLong %>% 
  group_by(study, isSig) %>% 
  summarize(cor = getCorPvalText(estimate, TOPMed)) %>%
  filter(isSig)
corrNovelSig$lab <- c("c", "d")

novelPlotAll<-novelLociRepLong %>% 
  ggplot(aes(x = TOPMed, y = estimate)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap( ~ study, scales = "free") +
  theme(legend.position="none") +
  geom_text(aes(x=-12, y=0.058, label = cor), 
            data = filter(corrNovelAll, study == "Dorajoo"), colour = "black") +
  geom_text(aes(x=-72, y=0.045, label = cor), 
            data = filter(corrNovelAll, study == "Li"), colour = "black") +
  geom_text(aes(x=28, y=-0.05, label = lab), 
            data = filter(corrNovelAll, study == "Dorajoo"), colour = "black") +
  geom_text(aes(x=27, y=-0.06, label = lab), 
            data = filter(corrNovelAll, study == "Li"), colour = "black") +
  geom_hline(yintercept=0, linetype="dashed", color = "red") +
  geom_vline(xintercept=0, linetype="dashed", color = "red") +
  xlab("TOPMed estimate (bp)") +
  ylab("Prior study estimate (T/S)") +
  theme_bw() +
  theme(legend.position="none")
print(novelPlotAll)

pdf("~/Research/OneDrive/telomere/gwas/manuscript/figPDFs/SuppFig2_NovelPlotAll.pdf", height=4, width=6)
print(novelPlotAll)
dev.off()


novelPlotSig<-novelLociRepLong %>% filter(isSig) %>%  
  ggplot(aes(x = TOPMed, y = estimate)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap( ~ study, scales = "free") +
  theme(legend.position="none") +
  geom_text(aes(x=-12, y=0.058, label = cor), 
            data = filter(corrNovelSig, study == "Dorajoo"), colour = "black") +
  geom_text(aes(x=-35, y=0.045, label = cor), 
            data = filter(corrNovelSig, study == "Li"), colour = "black") +
  geom_text(aes(x=28, y=-0.05, label = lab), 
            data = filter(corrNovelSig, study == "Dorajoo"), colour = "black") +
  geom_text(aes(x=27, y=-0.05, label = lab), 
            data = filter(corrNovelSig, study == "Li"), colour = "black") +
  geom_hline(yintercept=0, linetype="dashed", color = "red") +
  geom_vline(xintercept=0, linetype="dashed", color = "red") +
  xlab("TOPMed estimate (bp)") +
  ylab("Prior study estimate (T/S)") +
  theme_bw() +
  theme(legend.position="none")
print(novelPlotSig)

pdf("~/Research/OneDrive/telomere/gwas/manuscript/figPDFs/SuppFig2_NovelPlotSig.pdf", height=4, width=6)
print(novelPlotSig)
dev.off()

pdf("~/Research/OneDrive/telomere/gwas/manuscript/figPDFs/SuppFig2_NovelPlotBoth.pdf", height=8, width=6)
ggarrange(novelPlotAll, novelPlotSig, nrow = 2)
dev.off()


### For LocusZoom plots
### need to convert all to pngs to make a smaller supplementary figure
## pdf files were copied from ~/Research/OneDrive/telomere/gwas/results/LocusZooms/individual_plots

allFiles<-system("ls *.pdf /dcl01/mathias1/data/telomere_mtaub/gwas/results/LocusZooms/individual_plots", intern=TRUE)

for (currFile in allFiles){
  fileBase <- sub(".pdf", "", currFile)
  system(paste0("convert /dcl01/mathias1/data/telomere_mtaub/gwas/results/LocusZooms/individual_plots/", fileBase, ".pdf /dcl01/mathias1/data/telomere_mtaub/gwas/results/LocusZooms/individual_plots/", fileBase, ".png"))
}

