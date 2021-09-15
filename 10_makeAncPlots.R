library(tidyverse)
library(readxl)
library(GenomicRanges)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(rmeta)
library(meta)


### FOREST PLOTS
## use the multivariate regression results from Matt and Kruthika's forest plot code

allRes<-read.csv("~/Research/OneDrive/telomere/gwas/results/MERGED_allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop_p_lt_5e-9_conditionalAllSNPResSub_OASISAnno_jointResults_20200512.csv", header = TRUE, stringsAsFactors = FALSE)

colsToPull <- c("snpID", "LocusName.Final", "rsNum", "freq", "Qp", "Est_joint", "SE_joint", paste0(c("Est_joint_", "SE_joint_"), rep(c("white", "black", "hl", "asian"), each = 2)), paste0(c("freq"), c(".White", ".Black", ".HispanicLatino", ".Asian")))
allRes <- allRes %>% select(colsToPull) %>% dplyr::rename(Est_joint_all = Est_joint, SE_joint_all = SE_joint, freq_all = freq, freq_white = freq.White, freq_black = freq.Black, freq_hl = freq.HispanicLatino, freq_asian = freq.Asian)

allResLong<-left_join(left_join(allRes %>% select(snpID, LocusName.Final, rsNum, Qp, Est_joint_all, Est_joint_white, Est_joint_black, Est_joint_hl, Est_joint_asian) %>% pivot_longer(cols = Est_joint_all:Est_joint_asian, names_to = "ancestry", values_to = "Est") %>% mutate(ancestry = sub("Est_joint_", "", ancestry)), allRes %>% select(snpID, LocusName.Final, rsNum, Qp, SE_joint_all, SE_joint_white, SE_joint_black, SE_joint_hl, SE_joint_asian) %>% pivot_longer(cols = SE_joint_all:SE_joint_asian, names_to = "ancestry", values_to = "SE") %>% mutate(ancestry = sub("SE_joint_", "", ancestry))), allRes %>% select(snpID, LocusName.Final, rsNum, Qp, freq_all, freq_white, freq_black, freq_hl, freq_asian) %>% pivot_longer(cols = freq_all:freq_asian, names_to = "ancestry", values_to = "AAF") %>% mutate(ancestry = sub("freq_", "", ancestry)))
allResLong <- allResLong %>% mutate(Est = 1000*Est, SE = 1000*SE, CIlow = Est - 1.96*SE, CIhigh = Est + 1.96*SE, AAFPct = ifelse(!is.na(AAF), paste0(format(round(100*AAF,1), digits = 2), "%"), NA), AAFPct = ifelse(AAFPct == " 0.0%", "<0.1%", AAFPct))

allResLong <- allResLong %>% mutate(ancestry = factor(ancestry, levels=c("all", "white", "black", "hl", "asian"), labels = c("Pooled", "European", "African", "Hispanic/Latino", "Asian")), plotLab = paste(rsNum, snpID, LocusName.Final), I2Lab = paste0("Cochran's Q p-value: ", round(Qp, 2)))


lapply(split(allResLong, allResLong$plotLab), function(tmpSub){
  print(unique(tmpSub$plotLab))
  tableText <- tmpSub %>% mutate(Est = ifelse(is.na(Est), "", format(Est, digits=3))) %>% select(ancestry, Est, AAFPct) %>% as.matrix()
  colnames(tableText) <- NULL
  tableText <- rbind(rbind(rbind(c(unique(tmpSub$plotLab), NA, NA), c("Ancestry", "Est (bp)", "AAF")), tableText), c(unique(tmpSub$I2Lab), NA, NA))
  
#  forestplot(labeltext = tableText,mean = c(NA, NA, tmpSub$Est),lower = c(NA, NA, tmpSub$CIlow),upper = c(NA, NA, tmpSub$CIhigh),zero=0,is.summary=c(TRUE,TRUE,TRUE, rep(FALSE,4)), xlog=F, col=meta.colors(box=c("red", "royalblue", "red", "royalblue", "green"),line=c("red", "darkblue", "red", "darkblue"), summary="royalblue"))

  pdf(file = paste0("~/Research/OneDrive/telomere/gwas/results/figures/forestPlots/", unique(tmpSub$LocusName.Final), "_", unique(tmpSub$rsNum) ,".pdf"), height = 4, width=8)
  forestplot(labeltext = tableText,mean = c(NA, NA, tmpSub$Est, NA),lower = c(NA, NA, tmpSub$CIlow, NA),upper = c(NA, NA, tmpSub$CIhigh, NA),zero=0,is.summary=c(TRUE,TRUE,TRUE, rep(FALSE,4), TRUE), xlog=F, col=meta.colors(box="royalblue",line="darkblue", summary="royalblue"), graphwidth = unit(3, "inches"), align = c("l", "l", "l"))
  dev.off()
})

## write text for input to latex file
plotTags <- paste0(allRes$LocusName.Final, "_", allRes$rsNum)
write(sapply(split(plotTags, rep(1:(length(plotTags)/2), each = 2)), function(x) paste0("\\layoutTwo{", x[1], "}{", x[2],"}")), file="/Users/mtaub/Research/OneDrive/telomere/gwas/code/ForestPlotLatexText.txt")

############################
####### Heterogeneity tests
############################


multiVarRes <- read.table("~/Research/OneDrive/telomere/gwas/results/hits_summary_with_ancestry_with_joint_20200512.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
multiVarRes <- multiVarRes %>% dplyr::rename(Est_joint_all = Est_joint, SE_joint_all = SE_joint)

multiVarResHet <- multiVarRes %>% mutate(w_all = 1/SE_joint_all^2, w_white = 1/SE_joint_white^2, w_black = 1/SE_joint_black^2, w_hl = 1/SE_joint_hl^2, w_asian = 1/SE_joint_asian^2)

w_mat <- multiVarResHet %>% select(w_all:w_asian)
theta_mat <- multiVarResHet %>% select(Est_joint_all, Est_joint_white, Est_joint_black, Est_joint_hl, Est_joint_asian)
se_mat <- multiVarResHet %>% select(SE_joint_all, SE_joint_white, SE_joint_black, SE_joint_hl, SE_joint_asian)
Q_all <- rowSums(w_mat * ((theta_mat - (rowSums(w_mat * theta_mat, na.rm=TRUE)/rowSums(w_mat, na.rm=TRUE)))^2), na.rm=TRUE)
I2_all <- pmax(0, (Q_all - (rowSums(!is.na(w_mat)) - 1))/Q_all)

multiVarResHet %>% filter(I2_all > 0.5)

Q_sub <- rowSums(w_mat[,-1] * ((theta_mat[,-1] - (rowSums(w_mat[,-1] * theta_mat[,-1], na.rm=TRUE)/rowSums(w_mat[,-1], na.rm=TRUE)))^2), na.rm=TRUE)
Q_sub_p <- 1-pchisq(Q_sub, df = rowSums(!is.na(w_mat[,-1])) - 1)
I2_sub <- pmax(0, (Q_sub - (rowSums(!is.na(w_mat[,-1])) - 1))/Q_sub)


multiVarResHet$Q <- Q_sub
multiVarResHet$Qp <- Q_sub_p
multiVarResHet$I2 <- I2_sub


multiVarResHet %>% filter(I2_sub > 0.5)
multiVarResHet %>% filter(I2_sub > 0.75)

forMetagen <- left_join(multiVarRes %>% select(snpID, Est_joint_white, Est_joint_black, Est_joint_hl, Est_joint_asian,  SE_joint_white, SE_joint_black, SE_joint_hl, SE_joint_asian) %>% pivot_longer(cols = Est_joint_white:Est_joint_asian, names_to = "ancestry", values_to = "Est") %>% mutate(ancestry = sub("Est_joint_", "", ancestry)), multiVarRes %>% select(snpID, Est_joint_white, Est_joint_black, Est_joint_hl, Est_joint_asian,  SE_joint_white, SE_joint_black, SE_joint_hl, SE_joint_asian) %>% pivot_longer(cols = SE_joint_white:SE_joint_asian, names_to = "ancestry", values_to = "SE") %>% mutate(ancestry = sub("SE_joint_", "", ancestry))) %>% select(snpID, ancestry, Est, SE)

I2CI<-sapply(split(forMetagen, forMetagen$snpID), function(x){
  testOut<-metagen(TE = x$Est[!is.na(x$Est)], seTE = x$SE[!is.na(x$Est)], studlab = x$ancestry[!is.na(x$Est)], comb.fixed = TRUE)
  return(paste0(round(testOut$I2, 2), " (", round(testOut$lower.I2,2), ", ", round(testOut$upper.I2,2), ")"))
  
}
)

Qp<-sapply(split(forMetagen, forMetagen$snpID), function(x){
  testOut<-metagen(TE = x$Est[!is.na(x$Est)], seTE = x$SE[!is.na(x$Est)], studlab = x$ancestry[!is.na(x$Est)], comb.fixed = TRUE)
  return(testOut$pval.Q)
  
}
)

I2CIDf<-data.frame(snpID = names(I2CI), I2CI = I2CI, stringsAsFactors = FALSE)

multiVarResHet <- left_join(multiVarResHet, I2CIDf)

write.csv(multiVarResHet, file="~/Research/OneDrive/telomere/gwas/results/hits_summary_with_ancestry_with_joint_20200512_wHetTest.csv", row.names=FALSE, quote=TRUE)


########################
### OLD CODE
########################

colnames(condALL)[colnames(condALL) == "Est_Primary"] <- "Est.ALL"
colnames(condALL)[colnames(condALL) == "Est.SE_Primary"] <- "Est.SE.ALL"

condALL <- condALL %>% mutate(plotLab = factor(plotLab, levels = condALL$plotLab, labels=condALL$plotLab))


condAllLong<-left_join(condALL %>% select(snpID, plotLab, Est.ALL, Est.White, Est.Black, Est.HispanicLatino, Est.Asian) %>% pivot_longer(cols = Est.ALL:Est.Asian, names_to = "ancestry", values_to = "Est") %>% mutate(ancestry = sub("Est.", "", ancestry)), condALL %>% select(snpID, plotLab, Est.SE.ALL, Est.SE.White, Est.SE.Black, Est.SE.HispanicLatino, Est.SE.Asian) %>% pivot_longer(cols = Est.SE.ALL:Est.SE.Asian, names_to = "ancestry", values_to = "Est.SE") %>% mutate(ancestry = sub("Est.SE.", "", ancestry)))
condAllLong <- condAllLong %>% mutate(CIlow = Est - 1.96*Est.SE, CIhigh = Est + 1.96*Est.SE)

condAllLong <- condAllLong %>% mutate(ancestry = factor(ancestry, levels=c("ALL", "White", "Black", "HispanicLatino", "Asian")))

png(filename = paste0("~/Research/OneDrive/telomere/gwas/results/MERGED_allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop_p_lt_1e-5_forestPlots.png"), height = 1000, width = 1000)
condAllLong %>% ggplot(aes(x = ancestry, y = Est, ymin = CIlow, ymax = CIhigh)) +
  geom_pointrange(aes(col = ancestry)) + 
  geom_hline(yintercept = 0, linetype = 2) +
  facet_wrap(~plotLab,strip.position="top") +
  #  facet_wrap(~plotLab,strip.position="top",scales = "free_y") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  coord_cartesian(ylim=c(-0.25, 0.25))
dev.off()

### HEATMAPS
## will create ancestry specific result heatmaps for p-values and effect sizes

locusDef <- read.csv("~/Research/OneDrive/telomere/gwas/results/old/locusDefinitions.csv", header=TRUE, stringsAsFactors = FALSE)
locusDefGR <- GRanges(seqnames = locusDef$chr, range = IRanges(start = locusDef$start, end = locusDef$end))

mergedRes <- read.csv(gzfile("~/Research/OneDrive/telomere/gwas/results/MERGED_allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop_p_lt_1e-5.csv.gz"), header=TRUE, stringsAsFactors = FALSE)
mergedResGR <- GRanges(seqnames = mergedRes$chr, range = IRanges(start = mergedRes$pos, width=1))

resOL <- findOverlaps(locusDefGR, mergedResGR)

mergedResSub <- mergedRes[subjectHits(resOL),]
mergedResSub$locusIndex <- queryHits(resOL)
mergedResSub$locusFINAL <- locusDef$locusFINAL[queryHits(resOL)]
mergedResSub <- mergedResSub %>% mutate(locusFINAL = factor(locusFINAL, levels = mergedResSub$locusFINAL[!duplicated(mergedResSub$locusFINAL)], labels=mergedResSub$locusFINAL[!duplicated(mergedResSub$locusFINAL)]))

pvalCols <- paste("Score.pval", c("ALL", "White", "Black", "HispanicLatino", "Asian"), sep=".")
pvalBreaks <- c(5e-9, 5e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 0.1)


col_fun_pvals = colorRamp2(pvalBreaks, rev(brewer.pal(n = 9, name = 'YlOrRd')))

lgd = Legend(col_fun = col_fun_pvals, title = "p-value", at = c(5e-9, 5e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 0.1))

groupings <- list(1:12, 13:18, 19:30, 31:37)

for (i in 1:length(groupings)){
  testDatDF <- mergedResSub  %>% filter(locusIndex %in% groupings[[i]])
  testDat <- testDatDF %>% select(pvalCols)
  colnames(testDat) <- sub("Score.pval.", "", colnames(testDat))
  testDat<- t(as.matrix(testDat))
  
  png(filename = paste0("~/Research/OneDrive/telomere/gwas/results/MERGED_allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop_p_lt_1e-5_pvalHeatmap_group", i, ".png"), height = 200, width = 1000)
  p<-Heatmap(testDat, cluster_rows = FALSE, cluster_columns = FALSE, col = col_fun_pvals, heatmap_legend_param = list(title = "p-value", at = c(5e-9, 5e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 0.1), color_bar = "discrete"), column_split = testDatDF$locusFINAL, column_title = "",na_col="black", gap = unit(3, "mm"))
  print(p)
  dev.off()
}


BetaCutoffs<-rev(c(-0.25, -0.1, -0.05, -0.02, 0.02, 0.05, 0.1, 0.25))


col_fun_pvals = colorRamp2(BetaCutoffs, rev(brewer.pal(n = 8, name = 'RdBu')))
lgd_beta = Legend(col_fun = col_fun_pvals, title = "beta", at = c(-0.25, -0.1, -0.05, -0.02, 0.02, 0.05, 0.1, 0.25))

betaCols<-paste("Est", c("ALL", "White", "Black", "HispanicLatino", "Asian"), sep=".")


for (i in 1:length(groupings)){
  testDatDF <- mergedResSub  %>% filter(locusIndex %in% groupings[[i]])
  testDat <- testDatDF %>% select(betaCols)
  colnames(testDat) <- sub("Est.", "", colnames(testDat))
  testDat<- t(as.matrix(testDat))
  
  png(filename = paste0("~/Research/OneDrive/telomere/gwas/results/MERGED_allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop_p_lt_1e-5_betaHeatmap_group", i, ".png"), height = 200, width = 1000)
  Heatmap(testDat, cluster_rows = FALSE, cluster_columns = FALSE, col = col_fun_pvals, heatmap_legend_param = list(title = "beta", at = c(-0.25, -0.1, -0.05, -0.02, 0.02, 0.05, 0.1, 0.25), color_bar = "discrete"), column_split = testDatDF$locusFINAL, column_title_rot = 90)
  dev.off()
  
  
}

