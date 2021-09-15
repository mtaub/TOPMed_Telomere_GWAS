## using files from /telomere/inputFiles
## This script will 
##  1. merge TL data with DCC manifest
##  2. Drop samples with read lengths < 150
##  3. Clean study name for UCSF_AF samples
##  4. Drop samples due to non-consent for this analysis
##  5. Drop samples that are duplicates or trios
##  6. Merge in genotype PCs (11) (LATER)
##  7. Merge in batch PCs (200) (LATER)



library(tidyverse)
library(readxl)


setwd("/Users/mtaub/Research/OneDrive - Johns Hopkins/telomere/gwas")

## updated with 12/20 file
dccDat<-read.table("inputFiles/freeze8_sample_annot_2020-01-17.txt", header=TRUE, sep="\t", stringsAsFactors = FALSE)
dccDat<-rename(dccDat, NWDID=sample.id)

## for now need to replace IPF data from this file with the old version
#dccDat <- dccDat %>% filter(study != "IPF")
## I also need to fix the unique_subject_id and subject_id values for the IPF samples
## will use old version of dcc data file
#dccDatOld<-read.table("inputFiles/freeze8_sample_annot_2019-10-08.txt", header=TRUE, sep="\t", stringsAsFactors = FALSE)
#dccDatOld <- dccDatOld %>% filter(study == "IPF") %>% rename(NWDID=sample.id)

#dccDat<-rbind(dccDat, dccDatOld)


## read in and combine full TL data from Josh
allRes1<-read.csv("inputFiles/parsed_telseq_output.csv", stringsAsFactors=FALSE)
dim(allRes1)  # 93219 (as of 11/19/18, no duplicate NWDIDs)
allRes1 <- allRes1 %>% dplyr::select(NWDID, LENGTH_ESTIMATE, base_coord, studyname)
allRes1$phase<-"Phase1"
allRes2<-read.csv("inputFiles/telseq_output_2019_12_03.csv", stringsAsFactors = FALSE)
dim(allRes2)
allRes2 <- allRes2 %>% dplyr::select(NWDID, LENGTH_ESTIMATE, base_coord, studyname)
allRes2$phase<-"Phase2"

sum(allRes2$NWDID %in% allRes1$NWDID)

allRes<-rbind(allRes1, allRes2)

## 885 samples from Josh's results that are not in manifest; DROP
nrow(allRes) - sum(allRes$NWDID %in% dccDat$NWDID)

allResMerge<-dccDat %>% left_join(allRes)

## drop samples that don't have TL estimate
sum(is.na(allResMerge$LENGTH_ESTIMATE))
allResMerge<-allResMerge %>% filter(!is.na(LENGTH_ESTIMATE))
dim(allResMerge)

## regress out batch PCs from this full set of data with 200 PCs; this is also what will be posted
## to dbGaP

load("inputFiles/nofail_svd.pcs.200.rda")
sum(allResMerge$NWDID %in% batchPCs$NWDID)


allResMerge <- allResMerge %>% left_join(batchPCs)
allResMerge %>% filter(is.na(PC1)) %>% group_by(study) %>% summarize(n=n()) %>% print(n=Inf)
allResMerge %>% filter(is.na(PC1)) %>% group_by(funding) %>% summarize(n=n()) %>% print(n=Inf)

## for all samples with batch PCs, regress out 200 of them
## doing this now to include all possible samples; these will eventually be uploaded to dbGaP
batchLMOut<-lm(LENGTH_ESTIMATE ~ ., data=select(allResMerge, LENGTH_ESTIMATE, PC1:PC200), na.action=na.exclude)
allResMerge$LENGTH_ESTIMATE_BATCHADJ<-resid(batchLMOut)

save(allResMerge, file="results/allResMerge_allSamples_PCAdj.rda")


## drop all samples with read length < 150
## read length
## will lose 908 + 1090 = 1998 samples with read length < 150
table(allResMerge$base_coord, useNA="always")
allResMerge %>% filter(base_coord < 150) %>% group_by(study) %>% summarise(n=n())
allResMerge<-allResMerge %>% filter(base_coord > 150) 
dim(allResMerge)

## according to Stephanie Gogarten's email of 10/30/19, 118 samples with topmed_project = "AFGen" and study = NA are from study "UCSF_AF" 
allResMerge %>% filter(is.na(study)) %>% group_by(topmed_project) %>% summarise(n=n())
allResMerge[is.na(allResMerge$study) & allResMerge$topmed_project == "AFGen", "study"]<-"UCSF_AF"
allResMerge %>% filter(is.na(study)) %>% group_by(topmed_project) %>% summarise(n=n())

table(allResMerge$study, allResMerge$exclude, useNA="always")

## CONSENT???
## In phase1, I dropped c("COPDGene.DS-CS-RD", "COPDGene.NA" )
## Should I do that here? I assume yet
table(allResMerge$consent, useNA="always")
allResMerge$studyConsent<-paste(allResMerge$study, allResMerge$consent, sep=".")
table(allResMerge$studyConsent)
consentDF<-data.frame(count=table(allResMerge$studyConsent, useNA="always")) 
#write.csv(consentDF, file="inputFiles/consentDF.csv", row.names=FALSE)

## 291 dropped for consent here
studyConsentToDrop<-c("COPDGene.DS-CS-RD", "COPDGene.NA" )
allResMerge$isConsented<-!allResMerge$studyConsent %in% studyConsentToDrop





allResMerge <- allResMerge %>% filter(isConsented, !exclude)
dim(allResMerge)


## duplicates from DCC
dccDups<-read.table("inputFiles/freeze8_duplicates_2019-12-20.txt", stringsAsFactors = FALSE, header=TRUE)


## controls are pairs with each other
# confirm controls only appear with other controls
table(dccDups$study1, dccDups$study2, useNA="always")[,"Control"]
table(dccDups$study1, dccDups$study2, useNA="always")["Control",]
dupsSub<-filter(dccDups, !(study1 == "Control" & study2 == "Control")) ## left with 2365 duplicates
dim(dupsSub)

## merge in TL estimates: want to preferentially keep subjects who actually have length estimates
dupsSub<-allResMerge %>% select(NWDID, LENGTH_ESTIMATE, seq_center, phase) %>% right_join(dupsSub, by=c("NWDID" = "ID1"))
dupsSub<-allResMerge %>% select(NWDID, LENGTH_ESTIMATE, seq_center, phase) %>% right_join(dupsSub, by=c("NWDID" = "ID2"))
dupsSub<-dupsSub %>% select(ID1=NWDID.y, est1=LENGTH_ESTIMATE.y, center1=seq_center.y, consentGroup1=phase.y, ID2=NWDID, est2=LENGTH_ESTIMATE.x, center2=seq_center.x, consentGroup2=phase.x, study1, study2, MZtwinID)



## there are some triplicates that I might need to deal with differently
isTripID1<-dupsSub$ID1[dupsSub$ID1 %in% dupsSub$ID2 | duplicated(dupsSub$ID1)]
isTripID2<-dupsSub$ID2[dupsSub$ID2 %in% dupsSub$ID1 | duplicated(dupsSub$ID2)]
trips<-dupsSub %>% filter(ID1 %in% c(isTripID1, isTripID2) | ID2 %in% c(isTripID1, isTripID2))
dim(trips) ## 96

groups<-list(c(trips[1,"ID1"], trips[1, "ID2"]))
for(i in 2:nrow(trips)){
  whichMember<-which(sapply(groups, function(x) trips[i,"ID1"] %in% x | trips[i,"ID2"] %in% x))
  if(length(whichMember) == 0) groups<-c(groups, list(c(trips[i, "ID1"], trips[i, "ID2"])))
  if (!length(whichMember) == 0) groups[[whichMember]] <- unique(c(groups[[whichMember]], c(trips[i, "ID1"], trips[i, "ID2"])))
  if(length(whichMember) > 1) print("there is a problem") # if there are two different groups that the IDs are in because they should be merged
}
## want to drop all but one from each group, so will select one to keep and drop the rest
## complicated routine for selecting exactly one of a set of more than duplicates
set.seed(1239849)
tripsToDrop<-unlist(sapply(groups, function(x){
  tripsSub<-subset(trips, ID1 %in% x | ID2 %in% x)
  candsToKeep<-c(tripsSub$ID1[!is.na(tripsSub$est1)], tripsSub$ID2[!is.na(tripsSub$est2)]) # keep one that at least has an estimate
  toKeep<-NULL
  if(length(candsToKeep) >= 1) toKeep<-sample(candsToKeep, 1)
  toDrop<-setdiff(x, toKeep)
  return(toDrop)
}))
tripsToKeep<-sapply(groups, function(x) setdiff(x, tripsToDrop))
tripsToKeep<-unlist(tripsToKeep)

length(tripsToDrop) # 63 of which 43 are in data set
sum(tripsToDrop %in% allResMerge$NWDID)

## only want to drop duplicates where both of them have length estimates (otherwise, keep the one that has one and the other one is already excluded)
## if neither duplicate has an estimate, no need to keep them around
## actually I am going to change this -- these samples could link other samples that are in our data set, so will filter them out after making trips 
## but I have preferentially kept one of a group that has an estimate, so it does not matter now if I drop these
dupsSub<-dupsSub %>% filter(!(is.na(est1) & is.na(est2)))  ## 272 entries meet this criteria; now have 2093 duplicates
dim(dupsSub %>% filter((is.na(est1) | is.na(est2))))  # 832 of these now that I am not dropping double hits

# now deal with the remaining duplicates: 2008 entries
# first, restrict to only those not involved in the triples/groups
dupsSubNoTrips<-dupsSub %>% filter(!ID1 %in% c(trips$ID1, trips$ID2) & !ID2 %in% c(trips$ID1, trips$ID2))
dim(dupsSubNoTrips)  # 2008
## then, drop any from the pair that do not have TL data
## have already dropped those where both have no TL data
sum(is.na(dupsSubNoTrips$est1) & is.na(dupsSubNoTrips$est2))
idsToDrop<-c(dupsSubNoTrips$ID1[is.na(dupsSubNoTrips$est1)], dupsSubNoTrips$ID2[is.na(dupsSubNoTrips$est2)])
length(idsToDrop)  # 811
idsToKeep<-c(dupsSubNoTrips$ID2[is.na(dupsSubNoTrips$est1)], dupsSubNoTrips$ID1[is.na(dupsSubNoTrips$est2)])
length(idsToKeep)  # 811

dupsSubNoTripsNoNAs<-dupsSubNoTrips %>% filter(!is.na(est1) & !is.na(est2))
dim(dupsSubNoTripsNoNAs)  # 1197 (1197 + 811 = 2008)

# then randomly drop one from each remaining pair
set.seed(92387489)
idx<-sample(c(1,2), size=nrow(dupsSubNoTripsNoNAs), replace=TRUE)
idsToDrop<-c(idsToDrop, dupsSubNoTripsNoNAs$ID1[idx == 1], dupsSubNoTripsNoNAs$ID2[idx==2])
length(idsToDrop) # 2008
idsToKeep<-c(idsToKeep, dupsSubNoTripsNoNAs$ID2[idx == 1], dupsSubNoTripsNoNAs$ID1[idx==2])
length(idsToKeep) # 2008

all(c(dupsSub$ID1, dupsSub$ID2) %in% c(idsToKeep, idsToDrop, tripsToKeep, tripsToDrop))
length(unique(c(idsToKeep, idsToDrop, tripsToKeep, tripsToDrop)))
length(unique(c(idsToKeep, tripsToKeep)))
length(unique(c(idsToDrop, tripsToDrop)))
length(unique(c(dupsSub$ID1, dupsSub$ID2)))
# verify we are keeping one from each group
sapply(groups, function(x) sum(x %in% tripsToKeep))
all(c(idsToKeep, tripsToKeep) %in% allResMerge$NWDID)



## want to change phase in the full file  to DupRemoved
## not all ids in idsToDrop and tripsToDrop are in this file since we preferentially selected ones with missing data to remove: 855 total removed
allResMerge[allResMerge$NWDID %in% c(idsToDrop, tripsToDrop), "phase"]<-"DupRemoved"

table(allResMerge$phase, useNA="always")

dupCheck<-subset(allResMerge, phase != "DupRemoved")
sum(dupsSub$ID1 %in% dupCheck$NWDID & dupsSub$ID2 %in% dupCheck$NWDID)
sum(dupsSub$ID1 %in% allResMerge$NWDID & dupsSub$ID2 %in% allResMerge$NWDID)
idsRemoved<-subset(allResMerge, phase == "DupRemoved")$NWDID

allResMerge %>% group_by(study, phase) %>% summarise(n=n()) %>% spread(phase,n) %>% print(n=Inf)
allResMerge %>% group_by(study, seq_center,phase) %>% summarise(n=n()) %>% spread(phase,n) %>% print(n=Inf)

## want to check that there are no samples slated for analysis that still have duplicates
allResMerge[!is.na(allResMerge$unique_subject_key)  & allResMerge$phase %in% c("Phase1", "Phase2"),] %>% filter(unique_subject_key %in% unique_subject_key[duplicated(unique_subject_key)]) %>% arrange(unique_subject_key)

## no more duplicates remaining, based on unique_subject_key

## there are still some duplicate samples based on unique_subject_key, which were not in dups file
## As per 11/28 email with Stephanie Gogarten, I will drop both subjects since there is a data mapping issue
idsToDealWith<-allResMerge[!is.na(allResMerge$unique_subject_key)  & allResMerge$phase %in% c("Phase1", "Phase2"),] %>% filter(unique_subject_key %in% unique_subject_key[duplicated(unique_subject_key)]) %>% select(unique_subject_key) %>% unique(.) %>% as.vector(.)
idsToDealWith<-idsToDealWith[,1]
length(idsToDealWith)


## go ahead and drop the duplicates
allResMerge <- allResMerge %>% filter(phase != "DupRemoved")
dim(allResMerge)

save(allResMerge, file="results/allResMerge.rda")








