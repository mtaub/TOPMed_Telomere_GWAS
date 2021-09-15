## Using DCC harmonized age files where available


library(dplyr)
library(readxl)
library(tidyr)

## all DCC harmonized age files are available in 

setwd("/Users/mtaub/Research/OneDrive - Johns Hopkins/telomere/gwas")

colsToKeep<-c("SUBJECT_ID", "unique_subject_key", "topmed_study", "age_at_dna_blood_draw_wgs")
filesToRead<-list.files(path= paste0("inputFiles/phase2SourceFiles/"), pattern="unofficial.*")

allAge<-do.call(rbind, lapply(filesToRead, function(currFileName){
  currFile<-read.table(paste0("inputFiles/phase2SourceFiles/", currFileName), header=TRUE, stringsAsFactors = FALSE)
  currFile<-currFile[,colsToKeep]
}))
sort(unique(allAge$topmed_study))


## Need to add
## AFLMU : AFLMUIDs_MFS.csv
AFLMU = read_excel(paste0("inputFiles/phase2SourceFiles/aflmu_age-at-bleed-all.xlsx"))
## need to get SUBJECT_ID values to create unique_subject_key column


## Amish
  
Amish = read.delim(paste0("inputFiles/phase2SourceFiles/phs000956.v3.pht005002.v1.p1.c2.TOPMed_WGS_Amish_Subject_Phenotypes.HMB-IRB-MDS.txt"),stringsAsFactors=FALSE, skip = 10)
Amish <- Amish %>% mutate(unique_subject_key = paste0("Amish_", SUBJECT_ID), topmed_study="Amish") %>% select(SUBJECT_ID, unique_subject_key, topmed_study, age_at_dna_blood_draw_wgs=age_at_dna)
Amish$age_at_dna_blood_draw_wgs[Amish$age_at_dna_blood_draw_wgs == "90+"] <- "91"
Amish$age_at_dna_blood_draw_wgs<-as.numeric(Amish$age_at_dna_blood_draw_wgs)

  
## Australian Familial AF
## BioVU_AF
## Camp
Camp = read.csv(paste0("inputFiles/phase2SourceFiles/TOPMed_CAMP_age_at_draw.csv"),stringsAsFactors=FALSE)
Camp <- Camp %>% mutate(unique_subject_key = paste0("CAMP_", SUBJECT_ID), topmed_study="CAMP") %>% select(SUBJECT_ID, unique_subject_key, topmed_study, age_at_dna_blood_draw_wgs=age_at_dna_blood_draw)


## CARDIA
CARDIA<-read.table("inputFiles/phase2SourceFiles/unofficial_topmed_dcc_age_at_dna_blood_draw_wgs_withSampleID_freeze8_20200122_phs001612.txt", header = TRUE, stringsAsFactors = FALSE)
CARDIA<-CARDIA %>% rename(NWDID=SAMPLE_ID)
## want to drop existing CARDIA data
allAge <- allAge %>% filter(topmed_study != "CARDIA")



## CARE, PIMA, PCGC, PUSH

CAREetc<-read.csv("inputFiles/phase2SourceFiles/age_dna_CARE_PIMA_PCGC_PUSH.csv", header=TRUE, stringsAsFactors = FALSE)
CAREetc <- CAREetc %>% select(SUBJECT_ID, unique_subject_key, topmed_study = study, age_at_dna_blood_draw_wgs)

## CATHGEN
## CCAF

## CFS: new data uploaded: 
CFS<-read.table("inputFiles/phase2SourceFiles/CFS_Phase3_DNA_Sample_Ages_BEC_20200117.txt", header = TRUE, stringsAsFactors = FALSE)
CFS<-CFS %>% mutate(topmed_study="CFS") %>% select(SUBJECT_ID, unique_subject_key, topmed_study, age_at_dna_blood_draw_wgs = age_at_dna_blood_draw)

## CRA: new data from Cecelia
CRA<-read.csv(paste0("inputFiles/phase2SourceFiles/TOPMed_CRA_age_at_draw.csv"), stringsAsFactors = FALSE)
dim(CRA)
allAgeCRA<-allAge %>% filter(topmed_study == "CRA")
sum(allAgeCRA$SUBJECT_ID %in% CRA$SUBJECT_ID)
CRA<-CRA %>% left_join(allAgeCRA)

## in the new data set, not equal to what was in the old data set
CRA %>% filter(!is.na(age_at_dna_blood_draw), age_at_dna_blood_draw != age_at_dna_blood_draw_wgs)
table(is.na(CRA$age_at_dna_blood_draw), is.na(CRA$age_at_dna_blood_draw_wgs))

## use new file from Cecelia by default
CRA$age_at_dna_blood_draw_merge <- ifelse(!is.na(CRA$age_at_dna_blood_draw), CRA$age_at_dna_blood_draw, ifelse(!is.na(CRA$age_at_dna_blood_draw_wgs), CRA$age_at_dna_blood_draw_wgs, NA))

CRA %>% filter(!is.na(age_at_dna_blood_draw), !is.na(age_at_dna_blood_draw_wgs)) %>% mutate(ageDiff=age_at_dna_blood_draw - age_at_dna_blood_draw_wgs) %>% summarise(fracSmall=mean(abs(ageDiff)<1))

CRA <- CRA %>% mutate(unique_subject_key = paste0("CRA_", SUBJECT_ID))
  
CRA<-CRA %>% select(SUBJECT_ID, unique_subject_key, topmed_study, age_at_dna_blood_draw_wgs = age_at_dna_blood_draw_merge)
allAge <- allAge %>% filter(topmed_study != "CRA")

## EGCUT

## GALA/SAGE: new data from Angel, just the missing samples
GALASAGE<-read.table("inputFiles/phase2SourceFiles/phenotype_freeze_2020_0113_atgcv1.1_2020_0113_AM_topmedf8telomere.563.txt", header=TRUE, stringsAsFactors = FALSE, sep="\t")
allAgeGALASAGE<-allAge %>% filter(topmed_study %in% c("GALAI", "GALAII", "SAGE"))
## no repeats in there, so should just be able to select and merge new data with old
sum(GALASAGE$unique_subject_key %in% allAgeGALASAGE$unique_subject_key)

GALASAGE<-GALASAGE %>% select(SUBJECT_ID = NWDID, unique_subject_key, topmed_study = study, age_at_dna_blood_draw_wgs = age)

## GENAF
## GGAF
## INSPIRE_AF

## IPF
# issues with ID swaps not being reflected in DCC file
IPFAge<-read.table("inputFiles/phase2SourceFiles/IPF_HHWG_TELOMERES_PHENOTYPE_20200115_JHC.txt", header=TRUE, stringsAsFactors = FALSE, sep="\t")
IPFAge <- IPFAge %>% mutate(unique_subject_key = paste0("IPF_", ID), topmed_study="IPF") %>% select(SUBJECT_ID=ID, unique_subject_key, topmed_study, age_at_dna_blood_draw_wgs=AgeAtBloodDraw)
allAge <- allAge %>% filter(!topmed_study == "IPF")

## JHU_AF
## LTRC

LTRC<-read.delim(paste0("inputFiles/phase2SourceFiles/LTRC_R040_age_dna_draw.txt"),stringsAsFactors=FALSE)
LTRC <- LTRC %>% mutate(unique_subject_key = paste0("LTRC_", SUBJECT_ID), topmed_study="LTRC") %>% select(SUBJECT_ID, unique_subject_key, topmed_study, age_at_dna_blood_draw_wgs=age_at_dna_blood_draw)

## MGH_AF
## miRhythm
## MPP
## Partners
## PCGC_CHD
## PIMA
## PMBB_AF
## PUSH_SCD

## SAFS: new data uploaded, just the missing samples
SAFS<-read.table("inputFiles/phase2SourceFiles/age_at_blood_draw_for_margaret.tsv", header=TRUE, stringsAsFactors = FALSE, sep="\t")
allAgeSAFS<-allAge %>% filter(topmed_study == "SAFS")
# only missing samples uploaded, so can just merge in below
sum(SAFS$deidentified_subject %in% allAgeSAFS$SUBJECT_ID)

SAFS <- SAFS %>% mutate(unique_subject_key = paste0("SAFS_", deidentified_subject), topmed_study="SAFS")
SAFS <- SAFS %>% select(SUBJECT_ID = deidentified_subject, unique_subject_key, topmed_study, age_at_dna_blood_draw_wgs = age_at_blood_draw)

## wait on AFLMU and CARDIA until below
allAge<-rbind(allAge, Amish, Camp,LTRC, CRA, GALASAGE, SAFS, IPFAge, CFS, CAREetc)



load("./results/allResMerge.rda")

## need to map AFLMU NWDIDs to SUBJECT_IDs and unique_subject_key values
AFLMU<-allResMerge %>% filter(study=="AFLMU") %>% select(NWDID, unique_subject_key, SUBJECT_ID=subject_id) %>% right_join(AFLMU)
AFLMU<-AFLMU %>% mutate(topmed_study = "AFLMU") %>% select(SUBJECT_ID, unique_subject_key, topmed_study, age_at_dna_blood_draw_wgs = age)
allAge<-rbind(allAge, AFLMU)


CARDIAMerge<-allResMerge %>% filter(study=="CARDIA") %>% select(NWDID, unique_subject_key, SUBJECT_ID=subject_id) %>% right_join(CARDIA, by="NWDID")
CARDIAMerge<-CARDIAMerge %>% mutate(topmed_study = "CARDIA") %>% select(SUBJECT_ID = SUBJECT_ID.x, unique_subject_key=unique_subject_key.x, topmed_study, age_at_dna_blood_draw_wgs)
allAge<-rbind(allAge, CARDIAMerge)

allAge[allAge$unique_subject_key %in% allAge$unique_subject_key[duplicated(allAge$unique_subject_key)],]
sum(duplicated(allAge$unique_subject_key))
dim(allAge)
allAge <- allAge[!duplicated(allAge$unique_subject_key),]
dim(allAge)


## no more issues with duplicated unique_subject_key values except for control samples before merge
sum(duplicated(allResMerge$unique_subject_key))
allResMerge %>% filter(unique_subject_key %in% unique_subject_key[duplicated(unique_subject_key)]) %>% group_by(study) %>% summarise(n=n()) %>% print(n=Inf)

dim(allResMerge)
allResMerge <- allResMerge %>% left_join(allAge)
dim(allResMerge)
sum(duplicated(allResMerge$unique_subject_key))


allResMerge %>% group_by(study) %>% summarize(nAge=sum(!is.na(age_at_dna_blood_draw_wgs)), nNoAge=sum(is.na(age_at_dna_blood_draw_wgs))) %>% print(n=Inf)
dim(allResMerge)

## ancestry from HARE

## contains an object named out that has the updated hare groups
load("inputFiles/hareGroups_5xOutlier_20200219.RData")

hareAnc<-out %>% rename(NWDID = IID) %>% select(NWDID, hareFine_wProb, hareCoarse_wProb)

#hareAnc<-read.table("inputFiles/het_resid_var_group_v2.txt", header=TRUE, stringsAsFactors = FALSE)
#hareAnc <- hareAnc %>% rename(NWDID=IID)
#hareAnc <- hareAnc %>% mutate(justAnc = sub(".*\\_", "", het_resid_var_group))
#allResMerge<-allResMerge %>% left_join(hareAnc)

#hareFull<-read.table("inputFiles/hareGroups_18feb2020.txt",header = TRUE, stringsAsFactors = FALSE)
#hareFull <-hareFull %>% rename(NWDID = IID) 

hareProb<-read.table("inputFiles/HARE_probs_wIDs_input_v2.txt", header=TRUE, stringsAsFactors = FALSE)
hareProb <-hareProb %>% rename(NWDID = IID)

hareAncProb<-left_join(hareAnc, hareProb, by="NWDID")
#hareAncProb <- hareAncProb %>% left_join(hareFull)

## ancestry from DCC file
#dccDemo<-read.table("inputFiles/freeze8_demographic_annot_2020-01-17.txt", stringsAsFactors = FALSE, header=TRUE, sep="\t")
#dccDemo <- dccDemo %>% rename(NWDID = sample.id)
#allResMerge<-allResMerge %>% left_join(dccDemo)

## studies that are not in the DCC file
allResMerge %>% filter(!NWDID %in% hareAnc$NWDID) %>% group_by(study) %>% summarise(n=n()) %>% print(n=Inf)
hareAncProbStudy<-allResMerge %>% select(NWDID, study) %>% left_join(hareAncProb)
hareAncProbStudy <- hareAncProbStudy %>% mutate(studyAnc = paste(study, hareFine_wProb, sep="_"))
hareAncProbStudy[hareAncProbStudy$studyAnc == "FHS_NA", "studyAnc"] <- NA

## starting point of counts using HARE assignments
hareAncProbStudy %>% group_by(study) %>% count(studyAnc) %>% print(n=Inf)
## 787 samples are in groups < 30
hareAncProbStudy %>% group_by(study) %>% count(studyAnc) %>% filter(n<30) %>% ungroup(study) %>% summarise(total=sum(n))

ancGroups<-c("Dominican", "PuertoRican", "CentralAmerican", "Mexican", "SouthAmerican", "Cuban", "Black", "White", "Asian")
all(ancGroups == colnames(hareAncProbStudy)[5:13])
hareAncOrdered<-t(apply(hareAncProbStudy[,ancGroups], 1, function(x) ancGroups[order(x, decreasing=TRUE)]))
colnames(hareAncOrdered)<-c("First", "Second", "Third", "Fourth", "Fifth", "Sixth", "Seventh", "Eighth", "Ninth")

hareAncProbStudy <- cbind(hareAncProbStudy, hareAncOrdered, stringsAsFactors=FALSE)
hareAncProbStudy <- hareAncProbStudy %>% group_by(study, studyAnc) %>% summarise(nStart=n()) %>% right_join(hareAncProbStudy) 
as.data.frame(hareAncProbStudy %>% filter(nStart < 30, hareFine_wProb != First) )
hareAncProbStudy <- hareAncProbStudy %>% mutate(hareFineSecond = ifelse(nStart >= 30, hareFine_wProb, ifelse(hareFine_wProb != First, First, Second)))
hareAncProbStudy <- hareAncProbStudy %>% mutate(studyAncSecond = paste(study, hareFineSecond, sep="_"))
hareAncProbStudy %>% group_by(study) %>% count(studyAncSecond) %>% print(n=Inf)

## where are there still problems using second choice?
write.csv(hareAncProbStudy %>% group_by(study, studyAnc) %>% summarise(nInit=n()) %>% print(n=Inf), file="results/initHARECts.csv", row.names=FALSE, quote=FALSE)

write.csv(hareAncProbStudy %>% group_by(study, studyAncSecond) %>% summarise(nSecond=n()) %>% print(n=Inf), file="results/secondHARECts.csv", row.names=FALSE, quote=FALSE)

hareAncProbStudy %>% group_by(study, studyAncSecond) %>% summarise(nSecond=n()) %>% print(n=Inf)
hareAncProbStudy$studyAncAdjust<-hareAncProbStudy$studyAncSecond

## CAMP: can be fixed by making one hispanic group
hareAncProbStudy[hareAncProbStudy$studyAncSecond %in% c("CAMP_CentralAmerican", "CAMP_Cuban", "CAMP_Dominican", "CAMP_PuertoRican", "CAMP_SouthAmerican", "CAMP_Mexican"), "studyAncAdjust"]<-"CAMP_Hispanic"
## CARDIA: group lone Mexican with White (slightly higher probability)
hareAncProbStudy[hareAncProbStudy$studyAncSecond %in% c("CARDIA_Mexican"), "studyAncAdjust"]<-"CARDIA_White"
## CARE_BADGER: make one group
hareAncProbStudy[hareAncProbStudy$study %in% c("CARE_BADGER"), "studyAncAdjust"]<-"CARE_BADGER_All"
## CARE_CLIC: only 12 total people -- propose to merge with CARE_PACT
hareAncProbStudy[hareAncProbStudy$study %in% c("CARE_CLIC"), "studyAncAdjust"]<-"CARE_CLIC_PACT_All"
## CARE_PACT: on 22 total people -- propose to merge with CARE_CLIC
hareAncProbStudy[hareAncProbStudy$study %in% c("CARE_PACT"), "studyAncAdjust"]<-"CARE_CLIC_PACT_All"
## CARE_TREXA: make one group
hareAncProbStudy[hareAncProbStudy$study %in% c("CARE_TREXA"), "studyAncAdjust"]<-"CARE_TREXA_All"
## CCAF: group lone Puerto Rican with White
hareAncProbStudy[hareAncProbStudy$studyAncSecond %in% c("CCAF_PuertoRican"), "studyAncAdjust"]<-"CCAF_White"

## ChildrensHS_GAP: only 7 total people: group with ChildrensHS_IGERA Mexican or White as appropriate (5 mexican, 2 white)
hareAncProbStudy[hareAncProbStudy$studyAnc %in% c("ChildrensHS_GAP_Mexican", "ChildrensHS_GAP_CentralAmerican"), "studyAncAdjust"]<-"ChildrensHS_IGERA_GAP_Mexican"
hareAncProbStudy[hareAncProbStudy$studyAnc %in% c("ChildrensHS_GAP_White"), "studyAncAdjust"]<-"ChildrensHS_IGERA_GAP_White"
## ChildrensHS_IGERA: group all (4) non-Mexican with White
hareAncProbStudy[hareAncProbStudy$studyAncSecond %in% c("ChildrensHS_IGERA_White","ChildrensHS_IGERA_Black", "ChildrensHS_IGERA_CentralAmerican", "ChildrensHS_IGERA_SouthAmerican"), "studyAncAdjust"]<-"ChildrensHS_IGERA_GAP_White"
hareAncProbStudy[hareAncProbStudy$studyAncSecond %in% c("ChildrensHS_IGERA_Mexican"), "studyAncAdjust"]<-"ChildrensHS_IGERA_GAP_Mexican"

## ChildrensHS_MetaAir: make one group
hareAncProbStudy[hareAncProbStudy$study %in% c("ChildrensHS_MetaAir"), "studyAncAdjust"]<-"ChildrensHS_MetaAir_All"

## COPDGene: group two outliers with White
hareAncProbStudy[hareAncProbStudy$studyAncSecond %in% c("COPDGene_Dominican"), "studyAncAdjust"]<-"COPDGene_White"

## DECAF: only six people total: drop? **********
## ECLIPSE: group six outliers with White
hareAncProbStudy[hareAncProbStudy$studyAncSecond %in% c("ECLIPSE_CentralAmerican", "ECLIPSE_PuertoRican", "ECLIPSE_SouthAmerican"), "studyAncAdjust"]<-"ECLIPSE_White"

## FHS: two missing drop? one outlier group with White (i.e., make one group)
hareAncProbStudy[hareAncProbStudy$studyAncSecond %in% c("FHS_SouthAmerican"), "studyAncAdjust"]<-"FHS_White"

## GALAI: six outliers; group with Mexican or Puerto Rican, whichever has higher probability?
hareAncProbStudy[hareAncProbStudy$studyAncSecond %in% c("GALAI_Black", "GALAI_CentralAmerican",  "GALAI_SouthAmerican"), "studyAncAdjust"]<-"GALAI_Mexican"
hareAncProbStudy[hareAncProbStudy$studyAncSecond %in% c("GALAI_Cuban"), "studyAncAdjust"]<-"GALAI_PuertoRican"

## GALAII: eight outliers; group with whichever other group has highest probability?
hareAncProbStudy[hareAncProbStudy$studyAncSecond %in% c("GALAII_Black", "GALAII_Cuban", "GALAII_White"), "studyAncAdjust"]<-paste0("GALAII_", ancGroups[1:5][apply(hareAncProbStudy[hareAncProbStudy$studyAncSecond %in% c("GALAII_Black", "GALAII_Cuban", "GALAII_White"),ancGroups[1:5]], 1, which.max)])

## HVH: make one group
hareAncProbStudy[hareAncProbStudy$study %in% c("HVH"), "studyAncAdjust"] <- "HVH_All"

## HyperGEN: make one group
hareAncProbStudy[hareAncProbStudy$study %in% c("HyperGEN"), "studyAncAdjust"] <- "HyperGEN_Black"
## INSPIRE_AF: make one group
hareAncProbStudy[hareAncProbStudy$study %in% c("INSPIRE_AF"), "studyAncAdjust"] <- "INSPIRE_AF_White"
## IPF: make one group
hareAncProbStudy[hareAncProbStudy$study %in% c("IPF"), "studyAncAdjust"] <- "IPF_All"

## LTRC: two outliers: group one with white and one with black (by higher probability)
hareAncProbStudy[hareAncProbStudy$studyAncSecond == "LTRC_CentralAmerican", "studyAncAdjust"] <- "LTRC_Black"
hareAncProbStudy[hareAncProbStudy$studyAncSecond == "LTRC_Dominican", "studyAncAdjust"] <- "LTRC_White"

## Mayo_VTE: make one group
hareAncProbStudy[hareAncProbStudy$study %in% c("Mayo_VTE"), "studyAncAdjust"] <- "Mayo_VTE_All"

## MGH_AF: make one group
hareAncProbStudy[hareAncProbStudy$study %in% c("MGH_AF"), "studyAncAdjust"] <- "MGH_AF_All"

## MLOF: 19 Cuban and Dominican: put in next highest group (Third for all)
hareAncProbStudy[hareAncProbStudy$studyAncSecond %in% c("MLOF_Cuban", "MLOF_Dominican"), "studyAncAdjust"]<-paste0("MLOF_", hareAncProbStudy %>% filter(studyAncSecond %in% c("MLOF_Cuban", "MLOF_Dominican")) %>% pull(Third))


## OMG_SCD: make one group
hareAncProbStudy[hareAncProbStudy$study %in% c("OMG_SCD"), "studyAncAdjust"] <- "OMG_SCD_Black"

## PCGC_CHD: 13 outliers: put in next highest group among those with enough people (i.e. not Cuban)
hareAncProbStudy[hareAncProbStudy$studyAncSecond %in% c("PCGC_CHD_CentralAmerican", "PCGC_CHD_SouthAmerican"), "studyAncAdjust"]<-paste0("PCGC_CHD_", ancGroups[c(1,2,4,7:9)][apply(hareAncProbStudy[hareAncProbStudy$studyAncSecond %in% c("PCGC_CHD_CentralAmerican", "PCGC_CHD_SouthAmerican"),ancGroups[c(1,2,4,7:9)]], 1, which.max)])

## PharmHU: group all non-black for the second choice together
hareAncProbStudy[hareAncProbStudy$study %in% c("PharmHU") & !hareAncProbStudy$studyAncSecond %in% c("PharmHU_Black"), "studyAncAdjust"]<-"PharmHU_NonBlack"


## PIMA: make one group
hareAncProbStudy[hareAncProbStudy$study %in% c("PIMA"), "studyAncAdjust"] <- "PIMA_All"
## PMBB_AF: make one group
hareAncProbStudy[hareAncProbStudy$study %in% c("PMBB_AF"), "studyAncAdjust"] <- "PMBB_AF_White"
## PUSH_SCD: make one group
hareAncProbStudy[hareAncProbStudy$study %in% c("PUSH_SCD"), "studyAncAdjust"] <- "PUSH_SCD_All"
## SAFS: group one outlier with Mexican
hareAncProbStudy[hareAncProbStudy$studyAncSecond == "SAFS_Cuban", "studyAncAdjust"] <- "SAFS_Mexican"
## SAGE: make one group
hareAncProbStudy[hareAncProbStudy$study %in%  c("SAGE"), "studyAncAdjust"] <- "SAGE_All"
## SAPPHIRE_asthma: group 6 outliers with white (5) or Puerto Rican (1) (third choice)
hareAncProbStudy[hareAncProbStudy$studyAncSecond %in% c("SAPPHIRE_asthma_CentralAmerican", "SAPPHIRE_asthma_SouthAmerican"), "studyAncAdjust"]<-paste0("SAPPHIRE_asthma_", hareAncProbStudy %>% filter(studyAncSecond %in% c("SAPPHIRE_asthma_CentralAmerican", "SAPPHIRE_asthma_SouthAmerican")) %>% pull(Third))

## SARP: make one hispanic group from FIRST choice (studyAnc) ****
hareAncProbStudy[hareAncProbStudy$studyAnc %in% c("SARP_CentralAmerican", "SARP_Cuban", "SARP_Dominican", "SARP_PuertoRican", "SARP_SouthAmerican", "SARP_Mexican"), "studyAncAdjust"]<-"SARP_Hispanic"


## VU_AF: group 5 outliers with White
hareAncProbStudy[hareAncProbStudy$studyAncSecond %in% c("VU_AF_CentralAmerican", "VU_AF_Mexican", "VU_AF_SouthAmerican"), "studyAncAdjust"] <- "VU_AF_White"
## walk_PHaSST: make one group
hareAncProbStudy[hareAncProbStudy$study %in% c("walk_PHaSST"), "studyAncAdjust"] <- "walk_PHaSST_All"

## WHI: group 16 with third highest group
hareAncProbStudy[hareAncProbStudy$studyAncSecond %in% c("WHI_CentralAmerican", "WHI_SouthAmerican"), "studyAncAdjust"]<-paste0("WHI_", hareAncProbStudy %>% filter(studyAncSecond %in% c("WHI_CentralAmerican", "WHI_SouthAmerican")) %>% pull(Third))

hareAncProbStudy[hareAncProbStudy$study %in% c("DECAF"), "studyAncAdjust"] <- NA
hareAncProbStudy[hareAncProbStudy$studyAncSecond %in% c("FHS_NA"), "studyAncAdjust"] <- NA

write.csv(hareAncProbStudy %>% group_by(study, studyAncAdjust) %>% summarise(nFinal=n()) %>% print(n=Inf), file="results/finalHARECts.csv", row.names=FALSE, quote=FALSE)


allResMerge<-hareAncProbStudy %>% select(NWDID, hareFine_wProb, hareCoarse_wProb, studyAncAdjust) %>% right_join(allResMerge)
dim(allResMerge)


allResMerge %>% group_by(study) %>% summarize(nAnc = sum(!is.na(studyAncAdjust)), nNoAnc=sum(is.na(studyAncAdjust)), nAge=sum(!is.na(age_at_dna_blood_draw_wgs)), nNoAge=sum(is.na(age_at_dna_blood_draw_wgs))) %>% print(n=Inf)
dim(allResMerge)


write.csv(allResMerge %>% group_by(study) %>% summarize(nAnc = sum(!is.na(studyAncAdjust)), nNoAnc=sum(is.na(studyAncAdjust)), nAge=sum(!is.na(age_at_dna_blood_draw_wgs)), nNoAge=sum(is.na(age_at_dna_blood_draw_wgs))) %>% print(n=Inf), file="results/nAncAgeSummary_031519.csv", row.names=FALSE)

save(allResMerge, file="./results/allResMerge_AgeAnc.rda")



###################################
####### PCAir results
###################################

pcAir<-read.table("inputFiles/pcair_results.txt", stringsAsFactors = FALSE)
colnames(pcAir)<-c("NWDID", paste0("pcAir", 1:11))

sum(allResMerge$NWDID %in% pcAir$NWDID)

allResMerge <- allResMerge %>% left_join(pcAir)
## 9 from FHS, others are NA for study
allResMerge %>% filter(is.na(pcAir1)) %>% group_by(study) %>% summarize(n=n())


save(allResMerge, file="./results/allResMerge_AgeAncPCs.rda")

###################################
####### Filter and count by study/center/ancestry
###################################

dim(allResMerge)
## only keep complete observations
## for now, ignoring those missing ancestry
##allResMerge <- allResMerge %>% filter(!is.na(ancestry.group), !is.na(study), !(is.na(age_at_dna_blood_draw_wgs)), !is.na(PC1), !is.na(pcAir1), !is.na(sex))
allResMerge <- allResMerge %>% filter(!is.na(topmed_study), !(is.na(age_at_dna_blood_draw_wgs)), !is.na(studyAncAdjust), !is.na(PC1), !is.na(pcAir1), !is.na(sex))
dim(allResMerge) # 109122

allResMerge %>% count(study) %>% print(n=Inf)

## create a new seq_center_new variable to reset ILLUMINA to BROAD
allResMerge <- allResMerge %>% mutate(seq_center_new = ifelse(seq_center != "ILLUMINA", seq_center, "BROAD"))

save(allResMerge, file="./results/allResMerge_filteredCompleteData.rda")



###################################
####### Stop here for actual analysis
###################################


hareAncCoarse <- out %>% select(NWDID = IID, hareCoarse_wProb)
allResMerge <- allResMerge %>% left_join(hareAncCoarse)
allResMerge %>%  count(hareFine_wProb)
allResMerge %>%  count(hareCoarse_wProb)
allResMerge %>% ungroup(study) %>% count(hareFine_wProb)
allResMerge %>% ungroup(study) %>%  count(hareCoarse_wProb)


write.csv(allResMerge %>%  count(hareFine_wProb)%>%  pivot_wider(names_from=hareFine_wProb, values_from=n), file="./results/hareFine_ByStudy.csv", row.names=FALSE, quote=FALSE, na="-")

write.csv(allResMerge %>%  count(hareCoarse_wProb)%>%  pivot_wider(names_from=hareCoarse_wProb, values_from=n), file="./results/hareCoarse_ByStudy.csv", row.names=FALSE, quote=FALSE, na="-")




###################################
####### Stop here for actual analysis
###################################


allResMerge %>% group_by(study, ancestry.group) %>% summarise(n=n()) %>% spread(ancestry.group, n) %>% print(n=Inf)
allResMerge %>% group_by(study, seq_center, ancestry.group) %>% summarise(n=n()) %>% spread(ancestry.group, n) %>% print(n=Inf)

## for now, ignoring ancestry groups
##forAnalysis<-allResMerge %>% group_by(study, seq_center, ancestry.group) %>% summarise(n=n()) %>% filter(n > 50) %>% inner_join(allResMerge, by=c("study", "seq_center", "ancestry.group"))
forAnalysis<-allResMerge %>% group_by(study, seq_center) %>% summarise(n=n()) %>% filter(n > 50) %>% inner_join(allResMerge, by=c("study", "seq_center"))
dim(forAnalysis) # 100551
forAnalysis %>%  group_by(study, ancestry.group) %>% summarise(n=n()) %>% spread(ancestry.group, n) %>% print(n=Inf)
forAnalysis %>%  group_by(study, seq_center) %>% summarise(n=n()) %>% spread(seq_center, n) %>% print(n=Inf)

forAnalysis %>%  group_by(study, seq_center) %>% summarise(n=n()) %>% spread(seq_center, n) %>% mutate(total = sum(BAYLOR,BROAD, ILLUMINA, MACROGEN, NYGC, UW, WASHU, na.rm=TRUE)) %>% print(n=Inf)





### write out MESA results for Rebecca/Nathan/John
write.csv(allResMerge %>% filter(study == "MESA"), file="./results/allResMerge_MESA.csv", row.names=FALSE, quote=FALSE)





## check PharmHU
newPharmHU<-read.table("./inputFiles/allSourceFiles/PharmHU_HHWG_TELOMERES_PHENOTYPE_20181024_VS.txt", header=TRUE, stringsAsFactors = FALSE, sep="\t")
newPharmHU <- rename(newPharmHU, SUBJECT_ID=ID)
newPharmHUMerge<-allResMerge %>% filter(study=="PharmHU") %>% select(SUBJECT_ID=subject_id, age_at_dna_blood_draw_wgs) %>% right_join(newPharmHU)
## there are 8 more samples in the previous file, including all from the DCC file as well, but they are not in our data set
newPharmHUMerge %>% filter(is.na(age_at_dna_blood_draw_wgs), !is.na(AgeAtBloodDraw))
newPharmHUMerge %>% filter(!is.na(age_at_dna_blood_draw_wgs), is.na(AgeAtBloodDraw))

## check SAFS -- nothing new here
newSAFS<-read.table("./inputFiles/allSourceFiles/SAFS.txt", header=TRUE, stringsAsFactors = FALSE, sep="\t")
newSAFS <- rename(newSAFS, SUBJECT_ID=Deidentified.Subject)
newSAFSMerge<-allResMerge %>% filter(study=="SAFS") %>% select(SUBJECT_ID=subject_id, age_at_dna_blood_draw_wgs, topmed_phase) %>% right_join(newSAFS)
newSAFSMerge %>% filter(is.na(age_at_dna_blood_draw_wgs), !is.na(topmed_phase), !is.na(Age.At.Blood.Draw))

newSAFS2<-read.csv("~/Downloads/SAFS_have_age.csv", stringsAsFactors = FALSE)
newSAFS2<-newSAFS2 %>% rename(subject_id = SUBJECT_ID)
newSAFS2Merge<-allResMerge %>% filter(study == "SAFS") %>% right_join(newSAFS2)

## want to create data with SAFS indivdiuals missing age
## want to write out ID info for samples that are missing age for GALAI, GALAII and SAGE
SAFSMissing<-allResMerge %>% filter(study %in% c("SAFS"), is.na(age_at_dna_blood_draw_wgs)) %>% select(NWDID:ancestry.group)
write.csv(SAFSMissing, file="./results/SAFS_MissingAge.csv", row.names=FALSE, quote=FALSE)

CFSMissing<-allResMerge %>% filter(study %in% c("CFS"), is.na(age_at_dna_blood_draw_wgs)) %>% select(NWDID:ancestry.group)
write.csv(SAFSMissing, file="./results/CFS_MissingAge.csv", row.names=FALSE, quote=FALSE)

PharmHUMissing<-allResMerge %>% filter(study %in% c("PharmHU"), is.na(age_at_dna_blood_draw_wgs)) %>% select(NWDID:ancestry.group)
write.csv(SAFSMissing, file="./results/PharmHU_MissingAge.csv", row.names=FALSE, quote=FALSE)
