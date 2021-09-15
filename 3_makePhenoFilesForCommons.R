###################################
####### Create files for analysis on Analysis Commons
###################################


library(dplyr)
library(readxl)
library(tidyr)
library(tidyverse)

## all DCC harmonized age files are available in 

setwd("/Users/mtaub/Research/OneDrive - Johns Hopkins/telomere/gwas")

## most recent results from 031520: have corrected sequencing center for SAFS (seq_center_new variable)
load("./results/allResMerge_filteredCompleteData.rda")


## select only the columns that I need for the genetic analysis
forAnalysis <- allResMerge %>% select(NWDID, study, seq_center_new, sex, age_at_dna_blood_draw_wgs,studyAncAdjust, LENGTH_ESTIMATE_BATCHADJ, pcAir1:pcAir11)
dim(forAnalysis) # 109122 samples

save(forAnalysis, file="./results/allResMerge_forAnalysis_031520.rda")
write.csv(forAnalysis, file="./results/allResMerge_forAnalysis_031520.csv", row.names=FALSE, quote=FALSE)

## want to create ancestry-specific files using hareCoarse_wProb values

## need to load full HARE probabilities again since I want to check the distributions in the different groups
## contains an object named out that has the updated hare groups
load("inputFiles/hareGroups_5xOutlier_20200219.RData")
hareAnc<-out %>% rename(NWDID = IID) %>% select(NWDID, hareFine_wProb, hareCoarse_wProb)
hareProb<-read.table("inputFiles/HARE_probs_wIDs_input_v2.txt", header=TRUE, stringsAsFactors = FALSE)
hareProb <-hareProb %>% rename(NWDID = IID)
hareAncProb<-left_join(hareAnc, hareProb, by="NWDID")
hareAncProb <- hareAncProb %>% filter(NWDID %in% allResMerge$NWDID)
hareAncProb <- hareAncProb %>% mutate(HispanicLatino = Dominican + PuertoRican + CentralAmerican + Mexican + SouthAmerican + Cuban)
hareAncProb <- hareAncProb %>% mutate(analysisGroup = fct_collapse(hareFine_wProb, White = c("White", "Amish"), Asian = c("Asian", "Taiwanese", "HanChinese"), Black = c("Black"), HispanicLatino = c("CentralAmerican", "CostaRican", "Cuban", "Dominican", "Mexican", "PuertoRican", "SouthAmerican"), Brazilian = c("Brazilian"), Samoan = c("Samoan")))

hareAncProbLong <- hareAncProb %>% gather(key=probGroup, value=prob, Dominican:HispanicLatino)

hareAncProbLong %>% filter(analysisGroup == probGroup) %>% group_by(analysisGroup) %>% summarise(n=n(), count90 = sum(prob < 0.9, na.rm=TRUE), frac90 = mean(prob < 0.9, na.rm=TRUE), count70 = sum(prob < 0.7, na.rm=TRUE), frac70 = mean(prob < 0.7, na.rm=TRUE), count50 = sum(prob < 0.5, na.rm=TRUE), frac50 = mean(prob < 0.5, na.rm=TRUE))

## total to drop is 612 samples using cutoff of 0.7
hareAncProbLong %>% filter(analysisGroup == probGroup) %>% group_by(analysisGroup) %>% summarise(n=n(), count90 = sum(prob < 0.9, na.rm=TRUE), frac90 = mean(prob < 0.9, na.rm=TRUE), count70 = sum(prob < 0.7, na.rm=TRUE), frac70 = mean(prob < 0.7, na.rm=TRUE), count50 = sum(prob < 0.5, na.rm=TRUE), frac50 = mean(prob < 0.5, na.rm=TRUE)) %>% summarise(totDrop=sum(count70))
IDsToDrop<-hareAncProbLong %>% filter(analysisGroup == probGroup, prob < 0.7) %>% pull(NWDID) 

## 
forAnalysisByAncestry <- allResMerge %>% ungroup(study) %>% filter(!NWDID %in% IDsToDrop) %>% mutate(analysisGroup = fct_collapse(hareFine_wProb, White = c("White", "Amish"), Asian = c("Asian", "Taiwanese", "HanChinese"), Black = c("Black"), HispanicLatino = c("CentralAmerican", "CostaRican", "Cuban", "Dominican", "Mexican", "PuertoRican", "SouthAmerican"), Brazilian = c("Brazilian"), Samoan = c("Samoan")))

## White
forAnalysisWhiteAmish<-forAnalysisByAncestry %>% filter(analysisGroup =="White")  %>% select(NWDID, study, seq_center_new, sex, age_at_dna_blood_draw_wgs, LENGTH_ESTIMATE_BATCHADJ, pcAir1:pcAir11)
studyToJoinWhiteAmish<-forAnalysisWhiteAmish %>% count(study) %>% filter(n < 30) %>% pull(study)
forAnalysisWhiteAmish <- forAnalysisWhiteAmish %>% ungroup(study) %>% mutate(study = ifelse(study %in% studyToJoinWhiteAmish, "MergedStudy", study))
table(forAnalysisWhiteAmish$study, useNA="always")
## also need to reassign sequencing center for GeneSTAR to be BROAD instead of MACROGEN
forAnalysisWhiteAmish <- forAnalysisWhiteAmish %>% mutate(seq_center_new = ifelse(seq_center_new == "MACROGEN", "BROAD", seq_center_new))
save(forAnalysisWhiteAmish, file="./results/White_forAnalysis_031920.rda")
write.csv(forAnalysisWhiteAmish, file="./results/White_forAnalysis_031920.csv", row.names=FALSE, quote=FALSE)
## 85 samples in merged study

forAnalysisBlack<-forAnalysisByAncestry %>% filter(analysisGroup =="Black")   %>% select(NWDID, study, seq_center_new, sex, age_at_dna_blood_draw_wgs, LENGTH_ESTIMATE_BATCHADJ, pcAir1:pcAir11)
studyToJoinBlack<-forAnalysisBlack %>% count(study) %>% filter(n < 30) %>% pull(study)
forAnalysisBlack <- forAnalysisBlack %>% ungroup(study) %>% mutate(study = ifelse(study %in% studyToJoinBlack, "MergedStudy", study))
table(forAnalysisBlack$study, useNA="always")
save(forAnalysisBlack, file="./results/Black_forAnalysis_031920.rda")
write.csv(forAnalysisBlack, file="./results/Black_forAnalysis_031920.csv", row.names=FALSE, quote=FALSE)
## 59 samples in merged study

forAnalysisAsian<-forAnalysisByAncestry %>% filter(analysisGroup =="Asian")  %>% select(NWDID, study, seq_center_new, sex, age_at_dna_blood_draw_wgs, LENGTH_ESTIMATE_BATCHADJ, pcAir1:pcAir11)
studyToJoinAsian<-forAnalysisAsian %>% count(study) %>% filter(n < 30) %>% pull(study)
forAnalysisAsian <- forAnalysisAsian %>% ungroup(study) %>% mutate(study = ifelse(study %in% studyToJoinAsian, "MergedStudy", study))
table(forAnalysisAsian$study, useNA="always")
save(forAnalysisAsian, file="./results/Asian_forAnalysis_031920.rda")
write.csv(forAnalysisAsian, file="./results/Asian_forAnalysis_031920.csv", row.names=FALSE, quote=FALSE)
## 57 samples in merged study

forAnalysisBrazilian<-forAnalysisByAncestry %>% filter(analysisGroup =="Brazilian")  %>% select(NWDID, sex, age_at_dna_blood_draw_wgs, LENGTH_ESTIMATE_BATCHADJ, pcAir1:pcAir11)
save(forAnalysisBrazilian, file="./results/Brazilian_forAnalysis_031920.rda")
write.csv(forAnalysisBrazilian, file="./results/Brazilian_forAnalysis_031920.csv", row.names=FALSE, quote=FALSE)

forAnalysisSamoan<-forAnalysisByAncestry %>% filter(analysisGroup =="Samoan")  %>% select(NWDID, seq_center_new, sex, age_at_dna_blood_draw_wgs, LENGTH_ESTIMATE_BATCHADJ, pcAir1:pcAir11)
save(forAnalysisSamoan, file="./results/Samoan_forAnalysis_031920.rda")
write.csv(forAnalysisSamoan, file="./results/Samoan_forAnalysis_031920.csv", row.names=FALSE, quote=FALSE)

## more involved to deal with HL samples
forAnalysisHL<-forAnalysisByAncestry %>% filter(analysisGroup =="HispanicLatino")  
studyToJoinHL<-forAnalysisHL %>% count(study) %>% filter(n < 30) %>% pull(study)
forAnalysisHL <- forAnalysisHL %>% ungroup(study) %>% mutate(studyHL = ifelse(study %in% studyToJoinHL, "MergedStudy", study))
table(forAnalysisHL$study, useNA="always")
table(forAnalysisHL$studyHL, useNA="always")
## 212 samples in merged study

ancGroups<-c("Dominican", "PuertoRican", "CentralAmerican", "Mexican", "SouthAmerican", "Cuban", "Black", "White", "Asian")
forAnalysisHL <- forAnalysisHL %>% left_join(hareAncProb)
hareAncOrderedHL<-t(apply(forAnalysisHL[,ancGroups], 1, function(x) ancGroups[order(x, decreasing=TRUE)]))
colnames(hareAncOrderedHL)<-c("First", "Second", "Third", "Fourth", "Fifth", "Sixth", "Seventh", "Eighth", "Ninth")
forAnalysisHL <- cbind(forAnalysisHL, hareAncOrderedHL)

forAnalysisHL %>% group_by(study, studyAncAdjust) %>% count() %>% print(n=Inf)
forAnalysisHL %>% filter(study %in% studyToJoinHL) %>% group_by(study, studyAncAdjust) %>% count() %>% print(n=Inf)
forAnalysisHL %>% filter(study %in% studyToJoinHL) %>% group_by(study, hareFine_wProb) %>% count() %>% print(n=Inf)

## going to start from studyAncAdjust value for het variances, and combine for specific studies as needed
forAnalysisHL$studyAncAdjustHL<-forAnalysisHL$studyAncAdjust

## CAMP: 
forAnalysisHL[forAnalysisHL$study == "CAMP", "studyAncAdjustHL"]<-"CAMP_All"
## ChildrensHS_IGERA: 
forAnalysisHL[forAnalysisHL$study == "ChildrensHS_IGERA", "studyAncAdjustHL"]<-"ChildrensHS_IGERA_All"

## MergedStudy: Not sure yet...
forAnalysisHL %>% filter(studyHL == "MergedStudy") %>% group_by(hareFine_wProb) %>% count()
forAnalysisHL[forAnalysisHL$studyHL %in% c("MergedStudy"), "studyAncAdjustHL"]<-"MergedStudy_All"

## MLOF: combine black, white, asian
forAnalysisHL[forAnalysisHL$studyAncAdjust %in% c("MLOF_Asian", "MLOF_Black", "MLOF_White"), "studyAncAdjustHL"]<-"MLOF_NonHispanic"

## PCGC_CHD: Put Black and White with Mexican since biggest group; no clear easy alternative
forAnalysisHL %>% filter(studyAncAdjust %in% c("PCGC_CHD_Black", "PCGC_CHD_White")) %>% select(First:Ninth)
forAnalysisHL[forAnalysisHL$studyAncAdjust %in% c("PCGC_CHD_Black", "PCGC_CHD_White"), "studyAncAdjustHL"] <- "PCGC_CHD_Mexican"

## SAPPHIRE_asthma: put white with Mexican (biggest group)
forAnalysisHL %>% filter(studyAncAdjust %in% c("SAPPHIRE_asthma_White")) %>% select(First:Ninth)
forAnalysisHL[forAnalysisHL$studyAncAdjust %in% c("SAPPHIRE_asthma_White"), "studyAncAdjustHL"] <- "SAPPHIRE_asthma_Mexican"

## WHI: put white and black with Mexican (biggest group)
forAnalysisHL %>% filter(studyAncAdjust %in% c("WHI_Black", "WHI_White")) %>% select(First:Ninth)
forAnalysisHL[forAnalysisHL$studyAncAdjust %in% c("WHI_Black", "WHI_White"), "studyAncAdjustHL"] <- "WHI_Mexican"

forAnalysisHL %>% count(studyHL, studyAncAdjustHL) %>% print(n=Inf)

forAnalysisHL<- forAnalysisHL  %>% select(NWDID, study = studyHL, seq_center_new, sex, age_at_dna_blood_draw_wgs,studyAncAdjust = studyAncAdjustHL, LENGTH_ESTIMATE_BATCHADJ, pcAir1:pcAir11)
forAnalysisHL %>% count(study, seq_center_new) %>% print(n=Inf)
forAnalysisHL %>% count(study, studyAncAdjust) %>% print(n=Inf)

save(forAnalysisHL, file="./results/HL_forAnalysis_032320.rda")
write.csv(forAnalysisHL, file="./results/HL_forAnalysis_032320.csv", row.names=FALSE, quote=FALSE)


## extreme value analysis

lowerBound <- 
forAnalysisExtremes<- allResMerge %>% filter(!(LENGTH_ESTIMATE_BATCHADJ > quantile(LENGTH_ESTIMATE_BATCHADJ, 0.01) & LENGTH_ESTIMATE_BATCHADJ < quantile(LENGTH_ESTIMATE_BATCHADJ, 0.99)))
