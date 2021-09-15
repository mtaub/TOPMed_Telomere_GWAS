

library(qqman)
library(GenomicRanges)
library(readxl)
library(dplyr)
library(data.table)


makeManhattanNew<-function(pathToFile, fileBasename, rangesToPull, outputTag="", pathToOutput=pathToFile, annotateTable=NULL, yLim=c(0,55), makeQQ=FALSE, popTag="ALL"){
  allLoci<-fread(paste0("gunzip -c ", pathToFile, fileBasename, ".csv.gz"), header=TRUE, sep=",")


  if(makeQQ){
    ## make a qqplot
    source("~/Research/telomere/newManuscript/code/func.R")
    allLociRare<-allLoci %>% filter(freq < 0.01 | freq > 0.99)
    allLociCommon <- allLoci %>% filter(freq >= 0.01 & freq <= 0.99)
    bmp(filename=paste0(pathToOutput, fileBasename,  "_QQPlot_Rare01.bmp"), width=3.5, height=3.5, units='in', res = 200,  bg="white", type="cairo")
    fooRare<-my.qq.func(allLociRare$Score.pval, hm=2000, gc=TRUE) ## gc inflation factor
    my.qq.plot(fooRare, mt = paste0(popTag, ", MAF < 0.01 "), mt.cex = 0.8, lambda.cex = 0.5, cex.lab = 0.8, cex.axis = 0.8)
    dev.off()
    bmp(filename=paste0(pathToOutput, fileBasename,  "_QQPlot_Common01.bmp"), width=3.5, height=3.5, units='in', res = 200,  bg="white", type="cairo")
    fooCommon<-my.qq.func(allLociCommon$Score.pval, hm=2000, gc=TRUE) ## gc inflation factor
    my.qq.plot(fooCommon, mt = paste0(popTag, ", MAF >= 0.01 "), mt.cex = 0.8, lambda.cex = 0.5, cex.lab = 0.8, cex.axis = 0.8)
    dev.off()
  }
    
  
    
  ## select needed columns and filter by p-value
  allLociSub<-allLoci %>% filter(Score.pval < 0.01)
  allLociSub<-allLociSub %>% select(CHR=chr, BP=pos, P='Score.pval', SNP=snpID) %>% mutate(CHR = as.numeric(ifelse(CHR == "X", 23,CHR)))
  allLociRanges<-GRanges(seqnames=allLociSub$CHR, ranges=IRanges(start=allLociSub$BP, width=1))
  
  sigInPrev<-findOverlaps(allLociRanges, rangesToPull, select="arbitrary", ignore.strand=TRUE)
  
  ## Ras has requested only red for all prev loci
  ## will need to add blue for suggestive novel loci
  highlight.red<-allLociSub$SNP[which(!is.na(findOverlaps(allLociRanges, rangesToPull[rangesToPull$color == "Red"], select="arbitrary", ignore.strand=TRUE)))]
  highlight.blue<-allLociSub$SNP[which(!is.na(findOverlaps(allLociRanges, rangesToPull[rangesToPull$color == "Blue"], select="arbitrary", ignore.strand=TRUE)))]
  #highlight.green<-allLociSub$SNP[which(!is.na(findOverlaps(allLociRanges, rangesToPull[rangesToPull$color == "green"], select="arbitrary", ignore.strand=TRUE)))]
  
  
  bmp(filename=paste0(pathToOutput, fileBasename, outputTag,  "_Manhattan.bmp"), width=7.3, height=3.5, units='in', res = 200,  bg="white", type="cairo")
  par(mai=c(0.62,0.62,0.05,0.05),
      cex=0.50) 
  manhattanCMAT(allLociSub,col = c("gray50", "gray80"), 
                highlight.red=highlight.red,
                highlight.blue=highlight.blue,
                #highlight.purple=highlight.purple,
                #ylim=yLim,
                #highlight.green=highlight.green,
                cexPoints=0.4, alphaLevel=1, annotateTable=annotateTable,
                suggestiveline=FALSE, genomewideline=-log10(5e-9), main=popTag, chrlabs=c(1:22, "X"))
  dev.off()
}

resForPlot<-read_excel("/dcl01/mathias1/data/telomere_mtaub/gwas/results/forRasikaOASISLookup-withnovelty_v2_20200512.xlsx")

annotateTable<-resForPlot  %>% mutate(chr = gsub("\\:.*", "", snpID), chr = as.numeric(ifelse(chr == "X", 23,chr)), pos = as.numeric(sapply(strsplit(snpID, split=":"), function(x) x[2]))) %>% dplyr::select(SNP=snpID, label=LocusName.Final, chr, pos, negKb = `"-Kb"`, posKb = `"+Kb"`, novelty) %>% filter(!is.na(negKb))
rangesToPull<-GRanges(seqnames=annotateTable$chr, ranges=IRanges(start=annotateTable$pos - 1000*annotateTable$negKb, end=annotateTable$pos + 1000*annotateTable$posKb))
rangesToPull$color<-ifelse(is.na(annotateTable$novelty), "Red", "Blue")
## need to widen just the region on chr3
start(ranges(rangesToPull)[5]) <- start(ranges(rangesToPull)[5]) - 101000
end(ranges(rangesToPull)[5]) <- end(ranges(rangesToPull)[5]) + 2950000
annotateTable<-annotateTable %>% dplyr::select(SNP, label)


source("~/Research/telomere/newManuscript/code/myManhattan.R")

basePath<-"/dcl01/mathias1/data/telomere_mtaub/gwas/results/"

## ALL samples
makeManhattanNew(paste0(basePath, "ALL_GWASResults/"), "allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop", outputTag="_Labels", rangesToPull=rangesToPull, annotateTable=annotateTable)
makeManhattanNew(paste0(basePath, "White_GWASResults/"), "White_allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop", outputTag="_Labels", rangesToPull=rangesToPull, annotateTable=annotateTable, popTag="White")
makeManhattanNew(paste0(basePath, "Black_GWASResults/"), "Black_allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop", outputTag="_Labels", rangesToPull=rangesToPull, annotateTable=annotateTable, popTag="Black")
makeManhattanNew(paste0(basePath, "Asian_GWASResults/"), "Asian_allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop", outputTag="_Labels", rangesToPull=rangesToPull, annotateTable=annotateTable, popTag="Asian")
makeManhattanNew(paste0(basePath, "HispanicLatino_GWASResults/"), "HL_allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop", outputTag="_Labels", rangesToPull=rangesToPull, annotateTable=annotateTable, popTag="Hispanic/Latino")
makeManhattanNew(paste0(basePath, "Brazilian_GWASResults/"), "Brazilian_allChrs_telomere_adjagesexbatchPCs_minDP0_BRAVODepthDrop", outputTag="_Labels", rangesToPull=rangesToPull, annotateTable=annotateTable, popTag="Brazilian")
makeManhattanNew(paste0(basePath, "Samoan_GWASResults/"), "Samoan_allChrs_telomere_adjagesexseqctrbatchPCs_minDP0_BRAVODepthDrop", outputTag="_Labels", rangesToPull=rangesToPull, annotateTable=annotateTable, popTag="Samoan")


makeManhattanNew(paste0(basePath, "ALL_GWASResults/"), "allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop", outputTag="_noLabels", rangesToPull=rangesToPull, makeQQ=TRUE)
gc()
## White samples
makeManhattanNew(paste0(basePath, "White_GWASResults/"), "White_allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop", outputTag="_noLabels", rangesToPull=rangesToPull, makeQQ=TRUE, popTag="White")
gc()

## Black samples
makeManhattanNew(paste0(basePath, "Black_GWASResults/"), "Black_allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop", outputTag="_noLabels", rangesToPull=rangesToPull, makeQQ=TRUE, popTag="Black")
gc()

## Asian samples
makeManhattanNew(paste0(basePath, "Asian_GWASResults/"), "Asian_allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop", outputTag="_noLabels", rangesToPull=rangesToPull, makeQQ=TRUE, popTag="Asian")
gc()

## HispanicLatino samples
makeManhattanNew(paste0(basePath, "HispanicLatino_GWASResults/"), "HL_allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop", outputTag="_noLabels", rangesToPull=rangesToPull, makeQQ=TRUE, popTag="Hispanic/Latino")
gc()

## Brazilian samples
makeManhattanNew(paste0(basePath, "Brazilian_GWASResults/"), "Brazilian_allChrs_telomere_adjagesexbatchPCs_minDP0_BRAVODepthDrop", outputTag="_noLabels", rangesToPull=rangesToPull, makeQQ=TRUE, popTag="Brazilian")
gc()

## Samoan samples
makeManhattanNew(paste0(basePath, "Samoan_GWASResults/"), "Samoan_allChrs_telomere_adjagesexseqctrbatchPCs_minDP0_BRAVODepthDrop", outputTag="_noLabels", rangesToPull=rangesToPull, makeQQ=TRUE, popTag="Samoan")

