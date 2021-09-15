library(readxl)
library(GenomicRanges)
library(data.table)
library(dplyr)
library(tidyverse)

############################################
### Table 1: Peak SNPs with p < 5e-9, independent loci as identified from conditional analysis
############################################

## File created by Matt, merging all primary and conditional results for 59 variants
multiVarRes<-read.table("~/Research/OneDrive/telomere/gwas/results/hits_summary_with_ancestry_with_joint_20200512.txt", header = TRUE, stringsAsFactors = FALSE, sep="\t") %>% select(-novelty)

## File from me with all primary results (all + population subgroups) and some conditional results, with OASIS annotation
condRes<-read.csv(gzfile("~/Research/OneDrive/telomere/gwas/results/MERGED_allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop_p_lt_5e-9_conditionalAllSNPResSub_OASISAnno.csv.gz"), header=TRUE, stringsAsFactors = FALSE)

allRes<-left_join(multiVarRes, condRes)

## File from Matt with results from tests for heterogeneity added
hetRes <- read.csv("~/Research/OneDrive/telomere/gwas/results/hits_summary_with_ancestry_with_joint_20200512_wHetTest.csv", header = TRUE, stringsAsFactors = FALSE)

hetRes <- hetRes %>% select(snpID, Q, Qp, I2, I2CI)

allRes <- left_join(allRes, hetRes)

## File with some locus info, including novelty (do not use column from Matt's results)
moreInfo <- read_excel("~/Research/OneDrive/telomere/gwas/results/forRasikaOASISLookup-withnovelty_v2_20200512.xlsx")
moreInfo <- moreInfo %>% mutate(asterisk = ifelse(is.na(asterisk) , "", asterisk)) %>% select(snpID, asterisk, novelty)

allRes <- left_join(allRes, moreInfo)

allRes <- allRes %>% mutate(annotation = ifelse(Type %in% c("exonic", "intronic") & !(Function == ""), Function, RglmDB)) %>% mutate(annotation = ifelse(annotation == "", NA, annotation), Qp = ifelse(Qp == 1, NA, Qp))

allResForTable <- allRes %>% select(chr, pos, LocusName.Final, rsNum, asterisk, novelty, annotation,  MAC, pval_primary, PVE_Primary, pval_joint, pval_joint_white, pval_joint_black, pval_joint_hl, pval_joint_asian, Est_joint, Est_joint_white, Est_joint_black, Est_joint_hl, Est_joint_asian, Qp)

allResForTable <- allResForTable %>% mutate(asterisk = ifelse(pval_primary > 5e-9, "*", asterisk), Est_joint = 1000*Est_joint, Est_joint_white = 1000*Est_joint_white, Est_joint_black = 1000*Est_joint_black, Est_joint_hl = 1000*Est_joint_hl, Est_joint_asian = 1000* Est_joint_asian) %>% select(chr:rsNum, asterisk, novelty:Qp)

allResForTable <- allResForTable %>% arrange(as.numeric(chr), pos)

write.csv(allResForTable, file = "~/Research/OneDrive/telomere/gwas/manuscript/tableCSVs/Table1_20200519.csv", row.names=FALSE, quote = FALSE, na = "-")


write.csv(allRes %>% arrange(as.numeric(chr), pos), file = "~/Research/OneDrive/telomere/gwas/results/MERGED_allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop_p_lt_5e-9_conditionalAllSNPResSub_OASISAnno_jointResults_20200512.csv", row.names = FALSE, quote = TRUE)

#allResForNancy <- allRes %>% select(chr, pos, snpID, LocusName.Final, rsNum, ref, alt, Est_joint_pooled = Est_joint, SE_joint_pooled = SE_joint, pval_joint_pooled = pval_joint,  Est_joint_white, SE_joint_white, pval_joint_white, Est_joint_black, SE_joint_black, pval_joint_black, Est_joint_hl, SE_joint_hl, pval_joint_hl, Est_joint_asian, SE_joint_asian, pval_joint_asian ) %>% mutate_at(vars(Est_joint_pooled, SE_joint_pooled, Est_joint_white, SE_joint_white, Est_joint_black, SE_joint_black, Est_joint_hl, SE_joint_hl, Est_joint_asian, SE_joint_asian), .funs = funs(. * 1000)) %>% arrange(as.numeric(chr), pos)

#write.csv(allResForNancy, file = "~/Research/OneDrive/telomere/gwas/results/jointResultsForNancy_20200512.csv", row.names = FALSE, quote = FALSE)

#######################################################
### Table S1: Counts by study and ancestry group
#######################################################

load("~/Research/OneDrive/telomere/gwas/results/allResMerge_forAnalysis_031520.rda")
load("~/Research/OneDrive/telomere/gwas/results/White_forAnalysis_031920.rda")
load("~/Research/OneDrive/telomere/gwas/results/Black_forAnalysis_031920.rda")
load("~/Research/OneDrive/telomere/gwas/results/Asian_forAnalysis_031920.rda")
load("~/Research/OneDrive/telomere/gwas/results/HL_forAnalysis_032320.rda")

allSubGroups<-rbind(forAnalysisWhiteAmish %>% mutate(ancestryGroup = "White") %>% select(NWDID, ancestryGroup), rbind(forAnalysisBlack %>% mutate(ancestryGroup = "Black") %>% select(NWDID, ancestryGroup), rbind(forAnalysisHL %>% mutate(ancestryGroup = "Hispanic/Latino") %>% select(NWDID, ancestryGroup), forAnalysisAsian %>% mutate(ancestryGroup = "Asian") %>% select(NWDID, ancestryGroup))))

forAnalysis <- left_join(forAnalysis, allSubGroups)
forAnalysis <- forAnalysis %>% mutate(ancestryGroup = ifelse(is.na(ancestryGroup), "Other", ancestryGroup), ancestryGroup = factor(ancestryGroup, levels=c("White", "Black", "Hispanic/Latino", "Asian", "Other")))


demoTable<-left_join(forAnalysis %>% group_by(Study=study) %>% summarise('Total Count'=n(), 'Male Count (Pct)'=paste0(sum(sex == "M"), " (", round(mean(sex=="M"),2), ")"), 'Mean Age (SD, Range)'=paste0(round(mean(age_at_dna_blood_draw_wgs),0), " (", round(sd(age_at_dna_blood_draw_wgs),1), ", ", round(min(age_at_dna_blood_draw_wgs),0), "-", round(max(age_at_dna_blood_draw_wgs),0), ")")), forAnalysis %>% group_by(Study=study, ancestryGroup) %>% summarise(n=n()) %>% spread(ancestryGroup, n))
demoTable<-left_join(demoTable,forAnalysis %>% group_by(Study=study, seq_center_new) %>% summarise(n=n()) %>% spread(seq_center_new, n) )


## add row at bottom with totals
Total<-c("Total", ungroup(forAnalysis) %>% summarise('Total Count'=n(), 'Male Count (Pct)'=paste0(sum(sex == "M"), " (", round(mean(sex=="M"),2), ")"), 'Mean Age (SD, Range)'=paste0(round(mean(age_at_dna_blood_draw_wgs),0), " (", round(sd(age_at_dna_blood_draw_wgs),1), ", ", round(min(age_at_dna_blood_draw_wgs),0), "-", round(max(age_at_dna_blood_draw_wgs),0), ")")), ungroup(forAnalysis) %>% group_by(ancestryGroup) %>% summarise(n=n()) %>% spread(ancestryGroup, n), ungroup(forAnalysis) %>% group_by(seq_center_new) %>% summarise(n=n()) %>% spread(seq_center_new, n))

demoTable<-rbind(demoTable, unlist(Total)) %>% print(n=Inf)
write.csv(demoTable, "~/Research/OneDrive/telomere/gwas/manuscript/tableCSVs/SupplTable1_Demographics.csv", row.names=FALSE, quote=TRUE, na="-")


## want to pull studies used in previous paper
#load("~/Research/telomere/manuscript/results/samplesForReplication.rda")
#load("~/Research/telomere/manuscript/results/samplesForDiscovery_112918.rda")
#prevStudies <- unique(c(forDiscovery$study, forReplication$study))

#currStudies <- unique(forAnalysis$study)

#all(prevStudies %in% currStudies)
#setdiff(currStudies, prevStudies)

#studyInfo<- read_excel("~/Research/OneDrive/telomere/gwas/inputFiles/Access to study on dbGAP.xlsx")
#all(forAnalysis$study %in% studyInfo$Study)

#ctsByStudy <- forAnalysis %>% count(study) %>% mutate(newStudy = !study %in% prevStudies)
#ctsByStudy <- studyInfo %>% select(study = Study, Topmed_project, topmed_phs) %>% right_join(ctsByStudy)
#write.csv(ctsByStudy, file = "~/Research/OneDrive/telomere/gwas/manuscript/basicStudyInfo.csv", row.names=FALSE, quote= FALSE)

##############################
### Table S2: Conditional results
##############################

condResTable <- allRes %>% select("snpID", "chr", "pos", "LocusName.Final", "rsNum", "round_Cond", "Est_primary", "pval_primary", paste(c("Est_cond", "pval_cond"), rep(1:6, each = 2), sep = "_")) %>% mutate_at(vars(contains("Est_")), ~(.*1000)) %>% mutate(round_Cond = ifelse(is.na(round_Cond), "Primary", sub("cond", "Cond_", round_Cond)))

ldTable<-read.csv("~/Research/OneDrive/telomere/gwas/results/sentinels_pairwise_ld_all.csv", header=TRUE, stringsAsFactors = FALSE)
ldTable <- ldTable %>% dplyr::rename(snpID = other_snp)

condResTable <- condResTable %>% left_join(ldTable)
condResTable <- condResTable %>% mutate(r2 = ifelse(round_Cond == "Primary", NA, r2), dprime = ifelse(round_Cond == "Primary", NA, dprime)) %>% select(chr:round_Cond, r2, dprime, Est_primary:pval_cond_6)


write.csv(condResTable, "~/Research/OneDrive/telomere/gwas/manuscript/tableCSVs/SupplTable2_ConditionalRes.csv", row.names=FALSE, quote=TRUE, na="")


##############################
### Table S3: SNPs from previous studies
##############################

## code for pulling variants from primary results
#gunzip -c /dcl01/mathias1/data/telomere_mtaub/gwas/results/ALL_GWASResults/allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop.csv.gz | awk  'NR==FNR{a[$1,$2]="foo";next}; a[$2,$3]=="foo" {print}' FS=',' /dcl01/mathias1/data/telomere_mtaub/gwas/results/novel_loci/published_sentinels_unique_chrPos.csv FS=',' - > /dcl01/mathias1/data/telomere_mtaub/gwas/results/novel_loci/allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop_publishedSentinels.csv
#gunzip -c /dcl01/mathias1/data/telomere_mtaub/gwas/results/White_GWASResults/White_allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop.csv.gz | awk  'NR==FNR{a[$1,$2]="foo";next}; a[$2,$3]=="foo" {print}' FS=',' /dcl01/mathias1/data/telomere_mtaub/gwas/results/novel_loci/published_sentinels_unique_chrPos.csv FS=',' - > /dcl01/mathias1/data/telomere_mtaub/gwas/results/novel_loci/White_allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop_publishedSentinels.csv
#gunzip -c /dcl01/mathias1/data/telomere_mtaub/gwas/results/Black_GWASResults/Black_allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop.csv.gz | awk  'NR==FNR{a[$1,$2]="foo";next}; a[$2,$3]=="foo" {print}' FS=',' /dcl01/mathias1/data/telomere_mtaub/gwas/results/novel_loci/published_sentinels_unique_chrPos.csv FS=',' - > /dcl01/mathias1/data/telomere_mtaub/gwas/results/novel_loci/Black_allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop_publishedSentinels.csv
#gunzip -c /dcl01/mathias1/data/telomere_mtaub/gwas/results/HispanicLatino_GWASResults/HL_allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop.csv.gz | awk  'NR==FNR{a[$1,$2]="foo";next}; a[$2,$3]=="foo" {print}' FS=',' /dcl01/mathias1/data/telomere_mtaub/gwas/results/novel_loci/published_sentinels_unique_chrPos.csv FS=',' - > /dcl01/mathias1/data/telomere_mtaub/gwas/results/novel_loci/HL_allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop_publishedSentinels.csv
#gunzip -c /dcl01/mathias1/data/telomere_mtaub/gwas/results/Asian_GWASResults/Asian_allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop.csv.gz | awk  'NR==FNR{a[$1,$2]="foo";next}; a[$2,$3]=="foo" {print}' FS=',' /dcl01/mathias1/data/telomere_mtaub/gwas/results/novel_loci/published_sentinels_unique_chrPos.csv FS=',' - > /dcl01/mathias1/data/telomere_mtaub/gwas/results/novel_loci/Asian_allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop_publishedSentinels.csv

#gunzip -c /dcl01/mathias1/data/telomere_mtaub/gwas/results/Asian_GWASResults/Asian_allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop.csv.gz | head -1 > /dcl01/mathias1/data/telomere_mtaub/gwas/results/novel_loci/allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop_colNames.csv

prevSentinels <- read_excel("~/Research/OneDrive/telomere/gwas/results/ALL_GWASResults/novel_loci/published_sentinels-betas_PMCID.xlsx")  %>% rename(rsNum = SNP) %>% mutate(P = sub("−", "-", P))
prevSentinelsPos <- read.table("~/Research/OneDrive/telomere/gwas/results/ALL_GWASResults/novel_loci/pairwise-sentinel-LD-reqd_files/published_sentinels_positions.txt", header=TRUE, stringsAsFactors = FALSE, sep="\t") %>% rename(rsNum = Variation.ID) 

moreInfo <- read_excel("~/Research/OneDrive/telomere/gwas/results/forRasikaOASISLookup-withnovelty_v2_20200512.xlsx")
annotateTable<-moreInfo  %>% mutate(chr = gsub("\\:.*", "", snpID), chr = as.numeric(ifelse(chr == "X", 23,chr)), pos = as.numeric(sapply(strsplit(snpID, split=":"), function(x) x[2]))) %>% dplyr::select(SNP=snpID, label=LocusName.Final, chr, pos, negKb = `"-Kb"`, posKb = `"+Kb"`, novelty) %>% filter(!is.na(negKb))
rangesToPull<-GRanges(seqnames=annotateTable$chr, ranges=IRanges(start=annotateTable$pos - 1000*annotateTable$negKb, end=annotateTable$pos + 1000*annotateTable$posKb))
start(ranges(rangesToPull)[5]) <- start(ranges(rangesToPull)[5]) - 101000
end(ranges(rangesToPull)[5]) <- end(ranges(rangesToPull)[5]) + 2950000



#write.csv(prevSentinelsPos %>% mutate(chrPos = paste(Chromosome, Position, sep=":")) %>% filter(!duplicated(chrPos)) %>% select(Chromosome, Position), file = "~/Research/OneDrive/telomere/gwas/results/ALL_GWASResults/novel_loci/pairwise-sentinel-LD-reqd_files/published_sentinels_unique_chrPos.csv", row.names = FALSE, quote = FALSE)

# confirming only 57 unique entries, so will drop  any 
length(unique(prevSentinelsPos$rsNum))
length(unique(apply(prevSentinelsPos, 1, function(x) paste(x, collapse = ":"))))
prevSentinelsPos <- prevSentinelsPos %>% filter(!duplicated(rsNum))

prevSentinels <- left_join(prevSentinels, prevSentinelsPos)
prevSentinels <- prevSentinels %>% mutate(chrPos = paste(Chromosome, Position, sep=":"))

colNames <- read.csv("~/Research/OneDrive/telomere/gwas/results/ALL_GWASResults/novel_loci/pairwise-sentinel-LD-reqd_files/allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop_colNames.csv")
colNames <- names(colNames)
prevAll<-read.csv("~/Research/OneDrive/telomere/gwas/results/ALL_GWASResults/novel_loci/pairwise-sentinel-LD-reqd_files/allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop_publishedSentinels.csv", header=FALSE, stringsAsFactors = FALSE)
names(prevAll)<-c(colNames[1:4], paste0(colNames[5:14], ".ALL"), colNames[15:length(colNames)])
for (ancTag in c("White", "Black", "HL", "Asian")){
  prevCurr<-read.csv(paste0("~/Research/OneDrive/telomere/gwas/results/ALL_GWASResults/novel_loci/pairwise-sentinel-LD-reqd_files/", ancTag, "_allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop_publishedSentinels.csv"), header=FALSE, stringsAsFactors = FALSE)
  names(prevCurr)<-c(colNames[1:4], paste0(colNames[5:14], ".", ancTag), colNames[15:length(colNames)])
  prevAll <- left_join(prevAll, prevCurr)
}
prevAll <- prevAll %>% mutate(chrPos = paste(chr, pos, sep=":"))
prevAll <- prevAll %>% mutate(Est.ALL = 1000*Est.ALL, Est.White = 1000*Est.White, Est.Black = 1000*Est.Black, Est.HL = 1000*Est.HL, Est.Asian = 1000*Est.Asian)

prevSentinels <- left_join(prevSentinels, prevAll)
# there is one with multi alleles, but only one matches ref/alt from Dorajoo
prevSentinels %>% filter(chrPos == "8:73008648")
prevSentinels <- prevSentinels %>% filter(!snpID == "8:73008648:T:A")

## need to check on signs for betas and remove betas where EA is missing
prevSentinels <- prevSentinels %>% mutate(EA = toupper(EA))
prevSentinels %>% filter(EA != alt) %>% select(rsNum, EA, chr, pos, ref, alt, Beta, Est.ALL)
prevSentinels <- prevSentinels %>% mutate(Beta = ifelse(EA != alt, -Beta, ifelse(is.na(EA), NA, Beta)), Sign = ifelse(is.na(Beta), "N/A", ifelse(Beta < 0, "-", "+")), PMCID = substr(PMCID, 1, 10))

prevSentinelsGR <- GRanges(seqnames=prevSentinels$Chromosome, ranges=IRanges(start=prevSentinels$Position, width = 1))
table1PrevOL <- findOverlaps(prevSentinelsGR, rangesToPull)

prevSentinels$inTable1 <- "Yes"
prevSentinels$inTable1[setdiff(1:nrow(prevSentinels), queryHits(table1PrevOL))] <- "No"

# changed to only include trans-ethnic results
prevSentinelsForTable <- prevSentinels %>% select("chr", "pos", "rsNum", "ref", "alt", "GENE", "Author", "Year", "PMCID", "P", "Sign", "inTable1", paste(c("Score.pval", "Est", "freq"), rep(c("ALL"), each = 3), sep=".")) %>% arrange(chr, pos, Author)

prevSentinelsForTable <- prevSentinelsForTable %>% filter(!as.numeric(P) > 5e-8)

write.csv(prevSentinelsForTable, file = "~/Research/OneDrive/telomere/gwas/manuscript/tableCSVs/SupplTable3_published_sentinels_TOPMedResults.csv", row.names=FALSE, quote = TRUE, na = "-")

## ADD COLUMN: is this in Table 1 or not? Add to the ALL section (and check -- if it is exactly those with p < 5x10-9 or from loci from table 1. Question is really do we identify all these loci?? Or which do we not identify?)

## For replication table summary : correlation of effect sizes for replication of novel loci??

##############################
### Table S4: Annotation of Table 1 variants
##############################


annoColsToPull <- c("chr", "pos", "LocusName.Final", "rsNum", "ref", "alt", "Type", "Function", "AAchg", "GENE", "RglmDB", "eigenPC", "Dnase", "Reg", "SIFT",	"PP2_HDIV","PP2_HVAR",	"LRT",	"MT",	"MA",	"FATHMM",	"metaSVM",	"metaLR",	"PhyloP",	"SiPhy",	"GERP..",	"CADD",	"ClinVar",	"Phast")

write.csv(allRes %>% arrange(as.numeric(chr), pos) %>% select(annoColsToPull), file = "~/Research/OneDrive/telomere/gwas/manuscript/tableCSVs/SupplTable4_Annotation.csv", row.names=FALSE, quote = TRUE, na = "")

##############################
### Table S5: Replication of novel findings
##############################

# dorajoo
dorajoo<-fread("~/Research/OneDrive/telomere/gwas/results/ALL_GWASResults/novel_loci/replicated_loci/replication_dependent_files/dorajoo_replication.txt")
colnames(dorajoo)<-c("V1", "chr", "hg37pos", "other_allele", "test_allele", "pvalue", "beta", "se", "p_het")
#dorajoo has a bunch of SNPs where it's rsID:pos:ref:alt
#I need to separate it to just the rsID
library(stringr)
dorajoo$count<-str_count(dorajoo$V1, ':')
zero<-subset(dorajoo, dorajoo$count == 0)
colnames(zero)[1]<-"SNP"
three<-subset(dorajoo, dorajoo$count == 3)
three[,c("SNP", "pos", "ref", "alt"):=tstrsplit(V1, ":", fixed=T)]
three<-three[,c(11,2:10)]
dorajoo<-rbind(zero, three)
#there are some rsID's which say GSA-rsid# get rid of the prefix
dorajoo$SNP<-sub("GSA-", "", dorajoo$SNP)
dor<-dorajoo[,c(1:8)]
colnames(dor)<-c("rsNum", "dorajoo_chr", "dorajoo_hg37_pos", "dorajoo_NEA", "dorajoo_EA", "dorajoo_pvalue",
                 "dorajoo_beta", "dorajoo_se")

# li
li<-fread("~/Research/OneDrive/telomere/gwas/results/ALL_GWASResults/novel_loci/replicated_loci/replication_dependent_files/li_replication.txt")
colnames(li)<-c("rsNum", "li_chr", "li_hg37pos", "li_EA", "li_NEA", "li_EAF", "li_BETA", "li_SE", "li_pvalue",
                "n", "fdr")
li<-li[,c(1:9)]

repTable <- left_join(allResForTable, dor) %>% left_join(li) %>% filter(novelty == "Novel") %>% select( chr, pos, LocusName.Final, rsNum,Est_joint, li_pvalue, li_BETA, li_SE, dorajoo_pvalue, dorajoo_beta, dorajoo_se)
nTests<-repTable %>% filter(!is.na(li_pvalue) | !is.na(dorajoo_pvalue)) %>% nrow()
repTable <- repTable %>% mutate(threshold = ifelse(rowSums(cbind(li_pvalue, dorajoo_pvalue) < 0.05/nTests, na.rm = TRUE)>0, paste0("<", round(0.05/nTests,4)), ifelse(rowSums(cbind(li_pvalue, dorajoo_pvalue) < 0.05, na.rm = TRUE)>0, "<0.05", ifelse(rowSums(is.na(cbind(li_pvalue, dorajoo_pvalue)))==2, NA, ">0.05")))) 

write.csv(repTable, file = "~/Research/OneDrive/telomere/gwas/manuscript/tableCSVs/SupplTable5_replication_novel_TOPMedResults.csv", row.names=FALSE, quote = TRUE, na = "-")



##############################
### Table S6: Expanded population level results
##############################

colsForPrimaryTable <- c("chr", "pos", "LocusName.Final", "rsNum", paste0("Est_primary", c("", "_white", "_black", "_hl", "_asian")), paste0("pval_primary", c("", "_white", "_black", "_hl", "_asian")), paste0("freq", c("", ".White", ".Black", ".HispanicLatino", ".Asian")), paste0("PVE", c("_Primary", ".White", ".Black", ".HispanicLatino", ".Asian")))

primaryTable <- allRes %>% select(colsForPrimaryTable) %>% mutate_at(vars(Est_primary, Est_primary_white, Est_primary_black, Est_primary_hl, Est_primary_asian), .funs = funs(. * 1000)) %>% arrange(as.numeric(chr), pos)

write.csv(primaryTable, file = "~/Research/OneDrive/telomere/gwas/manuscript/tableCSVs/SupplTable6_primaryResults.csv", row.names=FALSE, quote = TRUE, na = "-")

##############################
### Table S9: BioVU PheWAS Results
##############################

phewasAAOld <- read_excel("~/Research/OneDrive/telomere/BioVU/AA_PheWAS_combined.xlsx")
phewasAA <- read_excel("~/Research/OneDrive/telomere/BioVU/re/AA_PheWAS_combined.xlsx")
phewasAA <- phewasAA %>% group_by(snp) %>% mutate(nTests = sum(!is.na(OR)), newBonf = ifelse(p < 0.05/nTests, TRUE, FALSE))
phewasAA <- phewasAA %>% group_by(phenotype) %>% mutate(nSnps = sum(!is.na(OR)))
phewasAASub<-phewasAA %>% filter(!is.na(OR)) 
# they match
table(phewasAASub$bonferroni, phewasAASub$newBonf)
phewasAASig <- phewasAASub %>% filter(newBonf) %>%  mutate(nSamps = paste(n_cases, n_controls, sep = " / "), ancGroup = "AA") %>% select(rsNum = snp, phenotype, group, description, OR, p, nSamps, ancGroup)


phewasBySNP<-phewasAASig %>% ungroup() %>% left_join(allRes %>% select(chr, pos, rsNum, LocusName.Final, ref, alt, Est_joint)) %>% mutate(Est_joint = 1000*Est_joint) %>% arrange(as.numeric(chr), pos, p) %>% select(chr, pos, LocusName.Final, rsNum, ref, alt, Est_joint, phenotype, group, description, OR, p, nSamps, ancGroup)

phewasEA <- read_excel("~/Research/OneDrive/telomere/BioVU/re/EA_PheWAS_combined.xlsx")
phewasEA <- phewasEA %>% group_by(snp) %>% mutate(nTests = sum(!is.na(OR)), newBonf = ifelse(p < 0.05/nTests, TRUE, FALSE))
phewasEA <- phewasEA %>% group_by(phenotype) %>% mutate(nSnps = sum(!is.na(OR)))
phewasEASub<-phewasEA %>% filter(!is.na(OR)) 
# they match
table(phewasEASub$bonferroni, phewasEASub$newBonf)
phewasEASig <- phewasEASub %>% filter(newBonf) %>%  mutate(nSamps = paste(n_cases, n_controls, sep = " / "), ancGroup = "EA") %>% select(rsNum = snp, phenotype, group, description, OR, p, nSamps, ancGroup)

phewasBySNPEA<-phewasEASig %>% ungroup() %>% left_join(allRes %>% select(chr, pos, rsNum, LocusName.Final, ref, alt, Est_joint)) %>% mutate(Est_joint = 1000*Est_joint) %>% arrange(as.numeric(chr), pos, p) %>% select(chr, pos, LocusName.Final, rsNum, ref, alt, Est_joint, phenotype, group, description, OR, p, nSamps, ancGroup)

#phewasEA <- phewasEA %>% mutate(OR = exp(beta), nSamps = NA, ancGroup = "EA") %>% select(rsNum, phenotype, group, description, OR, p, nSamps, ancGroup)
#phewasEA <- phewasEA %>% left_join(allRes %>% select(chr, pos, rsNum, LocusName.Final, ref, alt, Est_joint)) %>% mutate(Est_joint = 1000*Est_joint) %>% arrange(as.numeric(chr), pos, p) %>% select(chr, pos, LocusName.Final, rsNum, ref, alt, Est_joint, phenotype, group, description, OR, p, nSamps, ancGroup)

phewasBySNPBoth <- rbind(phewasBySNP, phewasBySNPEA) %>% arrange(as.numeric(chr), pos, p)

write.csv(phewasBySNPBoth, file = "~/Research/OneDrive/telomere/gwas/manuscript/tableCSVs/SupplTable9_Part2_PheWASBySNP.csv", row.names=FALSE, quote = TRUE, na = "-")


#phewasEAPRS <- read_delim("~/Research/OneDrive/telomere/BioVU/Telomere_EA_PRS_Phewas.txt", delim = "\t")
phewasEAPRS <- read_excel("~/Research/OneDrive/telomere/BioVU/re/EA_PRS_Scaled_PheWAS_Results.xlsx")
phewasEAPRSSub <- phewasEAPRS %>% filter(!is.na(OR)) 
phewasEAPRSSub <- phewasEAPRSSub %>% mutate(newBonf = ifelse(p < 0.05/nrow(phewasEAPRSSub), TRUE, FALSE))
# they match
table(phewasEAPRSSub$bonferroni, phewasEAPRSSub$newBonf)
phewasEAPRSSig <- phewasEAPRSSub %>% filter(newBonf) %>%  mutate(nSamps = paste(n_cases, n_controls, sep = " / "), ancGroup = "EA") %>% select(phenotype, group, description, OR, p, nSamps, ancGroup) %>% arrange(phenotype)

phewasAAPRS <- read_excel("~/Research/OneDrive/telomere/BioVU/re/AA_PRS_Scaled_PheWAS_Results.xlsx")
phewasAAPRSSub <- phewasAAPRS %>% filter(!is.na(OR)) 
phewasAAPRSSub <- phewasAAPRSSub %>% mutate(newBonf = ifelse(p < 0.05/nrow(phewasAAPRSSub), TRUE, FALSE))
# they match
table(phewasAAPRSSub$bonferroni, phewasAAPRSSub$newBonf)
phewasAAPRSSig <- phewasAAPRSSub %>% filter(newBonf) %>%  mutate(nSamps = paste(n_cases, n_controls, sep = " / "), ancGroup = "AA") %>% select(phenotype, group, description, OR, p, nSamps, ancGroup) 

## no significant AA PRS results
write.csv(phewasEAPRSSig, file = "~/Research/OneDrive/telomere/gwas/manuscript/tableCSVs/SupplTable9_Part1_PheWASByPRS.csv", row.names=FALSE, quote = TRUE, na = "-")

##############################
### Table SXX: BioVU LabWAS Results -- LEAVE OUT FOR NOW
##############################

## my bonferroni correction does not match hers here or for EAs

labwasAA <- read_csv("~/Research/OneDrive/telomere/BioVU/AA_LabWAS_Results_Combined.csv")
labwasAA <- labwasAA %>% group_by(Predictor) %>% mutate(nTests = sum(!is.na(OR)), newBonf = ifelse(p < 0.05/nTests, TRUE, FALSE))
table(labwasAA$Bonferroni, labwasAA$newBonf)
labwasAASig <- labwasAA %>% filter(newBonf) %>% mutate(ancGroup = "AA") %>% select(Predictor, Full_name, Group, N, p, OR, newBonf, ancGroup)


labwasEA <- read_csv("~/Research/OneDrive/telomere/BioVU/EA_LabWAS_Results_Combined.csv")
labwasEA <- labwasEA %>% group_by(Predictor) %>% mutate(nTests = sum(!is.na(OR)), newBonf = ifelse(p < 0.05/nTests, TRUE, FALSE))
table(labwasEA$Bonferroni, labwasEA$newBonf)
labwasEASig <- labwasEA %>% filter(newBonf) %>% mutate(ancGroup = "EA") %>% select(Predictor, Full_name, Group, N, p, OR, newBonf, ancGroup)

labwasAAPRS <- read_csv("~/Research/OneDrive/telomere/BioVU/Telomere_AA_PRS_LabWAS_Results.csv")
labwasAAPRS <- labwasAAPRS %>% mutate(newBonf = ifelse(p < 0.05/nrow(labwasAAPRS), TRUE, FALSE))
table(labwasAAPRS$Bonferroni, labwasAAPRS$newBonf) # nothing significant for AA PRS

labwasEAPRS <- read_csv("~/Research/OneDrive/telomere/BioVU/Telomere_EA_PRS_LabWAS_Results.csv")
labwasEAPRS <- labwasEAPRS %>% mutate(newBonf = ifelse(p < 0.05/nrow(labwasEAPRS), TRUE, FALSE))
table(labwasEAPRS$Bonferroni, labwasEAPRS$newBonf) 
labwasEAPRSSig <- labwasEAPRS %>% filter(newBonf) %>% mutate(ancGroup = "EA") %>% select(Full_name, Group, OR, p, N, ancGroup) %>% arrange(p)
write.csv(labwasEAPRSSig, file = "~/Research/OneDrive/telomere/gwas/manuscript/tableCSVs/SupplTable10_Part1_LabWASByPRS.csv", row.names=FALSE, quote = TRUE, na = "-")


labwasBySNP<-rbind(labwasEASig %>% ungroup(), labwasAASig %>% ungroup())  %>% rename(rsNum = Predictor) %>% left_join(allRes %>% select(chr, pos, rsNum, LocusName.Final, ref, alt, Est_joint)) %>% mutate(Est_joint = 1000*Est_joint) %>% arrange(as.numeric(chr), pos, p) %>% select(chr, pos, LocusName.Final, rsNum, ref, alt, Est_joint, Full_name, Group, OR, p, N, ancGroup)
write.csv(labwasBySNP, file = "~/Research/OneDrive/telomere/gwas/manuscript/tableCSVs/SupplTable10_Part2_LabWASBySNP.csv", row.names=FALSE, quote = TRUE, na = "-")



##############################
### Table S10: UKBiobank Results
##############################

## these results use a per-SNP bonferroni cutoff since that is what was done for BioVU
UKBBRes <- read.csv("~/Research/OneDrive/telomere/gwas/results/UKBB/UKBB_SAIGEResSigBySNP.csv", header= TRUE, stringsAsFactors = FALSE)

UKBBRes <- UKBBRes %>% rename(rsNum = ID) %>% left_join(allResForTable %>% select(rsNum, LocusName.Final, Est_joint)) %>% mutate(OR = exp(beta))

UKBBResTable <- UKBBRes %>% mutate(nSamps = paste(Ncases, Ncontrols, sep = " / ")) %>% select(X.CHROM, POS, LocusName.Final, rsNum, REF, ALT, Est_joint, PheCode, Group, Description, OR, pval, nSamps) %>% arrange(as.numeric(X.CHROM), POS, pval)

write.csv(UKBBResTable, file = "~/Research/OneDrive/telomere/gwas/manuscript/tableCSVs/SupplTable11_UKBBResults.csv", row.names=FALSE, quote = TRUE, na = "-")




##############################
### Additional checks of reported paper results
##############################

# count of significant loci
length(unique(allResForTable$LocusName.Final))

# count of novel loci
allResForTable %>% filter(!duplicated(LocusName.Final)) %>% count(novelty) %>% filter(novelty == "Novel")

# count of novel variants
allResForTable %>%  count(novelty) %>% filter(novelty == "Novel")


# how many loci replicated
repTable %>% group_by(LocusName.Final) %>% summarise(isSig = any(threshold %in% c("<0.0026", "<0.05"))) %>% count(isSig)

# count of loci with more than one signal
allResForTable %>% count(LocusName.Final) %>% filter(n > 1)

# count that directly match prior signals
# this is based on things labeled ** which was done by hand -- SHOULD BE DOUBLE CHECKED
allResForTable %>% count(asterisk)

# counts of different sample groups
nrow(forAnalysis)
nrow(forAnalysisWhiteAmish)
nrow(forAnalysisBlack)
nrow(forAnalysisHL)
nrow(forAnalysisAsian)
nrow(forAnalysis) - (nrow(forAnalysisWhiteAmish) + nrow(forAnalysisBlack) + nrow(forAnalysisHL) + nrow(forAnalysisAsian))

# percent male
mean(forAnalysis$sex == "M")

# age ranges
range(forAnalysis$age_at_dna_blood_draw_wgs)

# total number of variants computed as
#[mtaub@compute-060 ALL_GWASResults]$ gunzip -c allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0.csv.gz | wc -l
#162583130

# count of replicated prior loci and their names
# to get 25 prior loci group ("TERC" "LRRC34 (TERC)") ("OBFC1" "STN1 (OBFC1)") ("ZNF208" "ZNF676") ("RTEL1/STMN3" "RTEL1" "RTEL1/ZBTB46")
# 30 - 5 extras = 25
# so there are 16 total of the 25 prior that are replicated by us
unique(prevSentinelsForTable$GENE)
length(allResForTable %>% filter(!duplicated(LocusName.Final), is.na(novelty)) %>% pull(LocusName.Final))
sub("\\_", "\\/", paste(allResForTable %>% filter(!duplicated(LocusName.Final), is.na(novelty)) %>% pull(LocusName.Final), collapse = ", ") )
# this second one includes PRRC2A which we think is actually not necessarily included in the nearby HSPA1A signal (due to low LD between our sentinel and the previously reported SNV) and is missing ZNF257/ZNF676 just due to naming differences
unique(prevSentinelsForTable %>% filter(Score.pval.ALL<5e-9 | GENE %in% allResForTable$LocusName.Final) %>% pull(GENE))

# close to replicated loci; only SENP7 and CTC1 were not previoulsy mentioned
prevSentinelsForTable %>% filter(Score.pval.ALL < 0.05, Score.pval.ALL > 5e-9) %>% select(GENE, rsNum, Score.pval.ALL)

# distance between PRRC2A and HSPA1A
allResForTable %>% filter(LocusName.Final == "HSPA1A") %>% pull(pos) - prevSentinelsForTable %>% filter(GENE == "PRRC2A") %>% pull(pos)

# not replicated
sub("\\_", "\\/", paste(prevSentinelsForTable %>% filter(inTable1 == "No", Score.pval.ALL > 0.05) %>% pull(GENE), collapse = ", "))

# for known loci, how many have multiple independent SNVs?
allResForTable %>% filter(is.na(novelty)) %>% count(LocusName.Final) %>% filter(n > 1)
allResForTable %>% filter(is.na(novelty)) %>% count(asterisk) %>% filter(asterisk == "**") 

# number/list of loci with missense coding variants
allResForTable %>% filter(is.na(novelty), annotation == "missense") %>% pull(LocusName.Final)

# number of variants in known loci with RegulomeDB score < 7
allResForTable %>% filter(is.na(novelty)) %>% count(annotation)


# count of variants falling in novel loci
allResForTable %>% filter(novelty == "Novel") %>% nrow()

# replication results
# note: TYMS is getting counted twice
repTable %>% filter(!is.na(li_pvalue) | !is.na(dorajoo_pvalue)) %>% nrow()
repTable %>% filter(threshold == "<0.0026")
paste(repTable %>% filter(threshold == "<0.0026") %>% pull(LocusName.Final), collapse = ", ")
repTable %>% filter(threshold == "<0.05")
paste(repTable %>% filter(threshold == "<0.05") %>% pull(LocusName.Final), collapse = ", ")

# p-value for SAMHD1
allResForTable %>% filter(LocusName.Final == "SAMHD1") %>% pull(pval_joint)

# correlation numbers for replication results
repTable %>% summarise(corTOPMedLi = cor.test(Est_joint, li_BETA)$estimate, corPTOPMedLi = cor.test(Est_joint, li_BETA)$p.value, corTOPMedDorajoo = cor.test(Est_joint, dorajoo_beta)$estimate, corPTOPMedDorajoo = cor.test(Est_joint, dorajoo_beta)$p.value)

repTable %>% filter(!is.na(threshold), !threshold == ">0.05") %>%  summarise(corTOPMedLi = cor.test(Est_joint, li_BETA)$estimate, corPTOPMedLi = cor.test(Est_joint, li_BETA)$p.value, corTOPMedDorajoo = cor.test(Est_joint, dorajoo_beta)$estimate, corPTOPMedDorajoo = cor.test(Est_joint, dorajoo_beta)$p.value)

# min and max effect sizes of common and rare variants
allRes %>% filter(freq >= 0.05 & freq <= 0.95) %>% summarise(minPrimary = min(abs(Est_primary)), maxPrimary = max(abs(Est_primary))) 
allRes %>% filter(!(freq >= 0.05 & freq <= 0.95)) %>% summarise(minPrimary = min(abs(Est_primary)), maxPrimary = max(abs(Est_primary))) 

# which SNVs have p-values for Cochrane's Q < 0.05
allResForTable %>% filter(Qp < 0.05)

# group-specific p-values for TINF2
allResForTable %>% filter(rsNum == "rs28372734")
allResForTable %>% filter(rsNum == "rs8016076")

phewasEAPRSSig %>% count(group)
phewasEAPRSSig %>% filter(group == "neoplasms") %>% pull(OR)



##############################
### Additional data selections
##############################

## subset JHS, WHI and GeneSTAR
load("~/Research/OneDrive/telomere/gwas/results/allResMerge_forAnalysis_031520.rda")
forAnalysisSubStudies<-forAnalysis %>% filter(study %in% c("JHS", "WHI", "GeneSTAR"))
write.csv(forAnalysisSubStudies, file = "~/Research/OneDrive/telomere/gwas/results/allResMerge_forAnalysis_031520_JHSWHIGeneSTAR.csv", row.names = FALSE, quote=FALSE)

### Table S3: All SNPs with p < 5e-9 in at least one group??

mergedRes<-read.csv(gzfile("~/Research/OneDrive/telomere/gwas/results/MERGED_allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop_p_lt_5e-8_OASISAnno.csv.gz"), header=TRUE, stringsAsFactors = FALSE)

allWithLocusGroup<-read.csv("~/Research/OneDrive/telomere/gwas/results/ALL_GWASResults/allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop_p_lt_5e-8_withLocusGroups.csv", header = TRUE, stringsAsFactors = FALSE)

mergedRes <- allWithLocusGroup %>% select(locusGroup, newLOCUS, snpID) %>% right_join(mergedRes)
mergedRes$diffPos <- c(0,diff(mergedRes$pos))
mergedRes$sameChr <- c(TRUE, mergedRes$chr[1:(nrow(mergedRes)-1)] == mergedRes$chr[2:(nrow(mergedRes))] )
contigBreaks<-which(!mergedRes$sameChr  | mergedRes$diffPos > 200000)
mergedRes[unique(sort(c(contigBreaks-1, contigBreaks))),1:5]

mergedRes$locusGroupMerged<-rep(1:(length(contigBreaks)+1), times=c(contigBreaks[1]-1, diff(contigBreaks), nrow(mergedRes)-contigBreaks[length(contigBreaks)]+1))
new8OutsideLoci<-mergedRes %>% filter(locusGroupMerged %in% (mergedRes %>% group_by(locusGroupMerged) %>% summarise(allNA = all(is.na(newLOCUS))) %>% filter(allNA==TRUE) %>% pull(locusGroupMerged))) %>% select("snpID", paste(c("freq", "MAC", "Score.pval"), rep(c("Black", "White", "Asian", "HispanicLatino", "Samoan", "Brazilian"), each=3), sep="."))

write.csv(new8OutsideLoci, file="~/Research/OneDrive/telomere/gwas/results/MERGED_allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop_p_lt_5e-8_forBRAVOCheck_notInLoci.csv", row.names=FALSE, quote=FALSE)

## check against what Rasika pulled
library(readxl)
rasResult<-read_excel("~/Research/OneDrive/telomere/gwas/results/ALL_GWASResults/OASIS QC/Chromosome_BRAVO_confirm.xlsx", sheet="Stratified Anlaysis Lookups")
colnames(rasResult)<-c("snpID", "snpIDMod", "Outcome", "Path")
rasResultOdd<- rasResult %>% filter(Outcome == "odd") %>% pull(snpID)

mergedResOdd<-mergedRes %>% filter(snpID %in% rasResultOdd) %>% select(paste(c("freq", "MAC", "Score.pval"), rep(c("Black", "White", "Asian", "HispanicLatino", "Samoan", "Brazilian"), each=3), sep="."))
sum(mergedResOdd$Score.pval.Black < 5e-8, na.rm=TRUE)
sum(mergedResOdd$Score.pval.White < 5e-8, na.rm=TRUE)
sum(mergedResOdd$Score.pval.Asian < 5e-8, na.rm=TRUE)
sum(mergedResOdd$Score.pval.HispanicLatino < 5e-8, na.rm=TRUE)
sum(mergedResOdd$Score.pval.Brazilian < 5e-8, na.rm=TRUE)
sum(mergedResOdd$Score.pval.Samoan < 5e-8, na.rm=TRUE)
mergedResOdd %>% select("freq.White", "MAC.White", "Score.pval.White")

sum(mergedResOdd$Score.pval.White < 5e-9, na.rm=TRUE)

## want to pull new positions for Rasika to check
new8<-mergedRes %>% filter(Score.pval.ALL > 5e-8) %>% pull(snpID)
write.csv(new8, file="~/Research/OneDrive/telomere/gwas/results/MERGED_allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop_p_lt_5e-8_forBRAVOCheck.csv", row.names=FALSE, quote=FALSE)

justPvals<-mergedRes[,grep("pval", colnames(mergedRes))]
nrow(justPvals) - sum(justPvals$Score.pval.ALL < 5e-8) ## 247
sum(rowSums(justPvals < 5e-9, na.rm=TRUE)>0)
idx9<-which(rowSums(justPvals < 5e-9, na.rm=TRUE)>0)

new9<-mergedRes[idx9,] %>% filter(Score.pval.ALL > 5e-9) %>% pull(snpID)
write.csv(new9, file="~/Research/OneDrive/telomere/gwas/results/MERGED_allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop_p_lt_5e-9_forBRAVOCheck.csv", row.names=FALSE, quote=FALSE)

justPvals <- justPvals %>% filter(rowSums(justPvals < 5e-9, na.rm=TRUE)>0)
nrow(justPvals) - sum(justPvals$Score.pval.ALL < 5e-9) ## 149
notSigAll <- justPvals %>% filter(justPvals$Score.pval.ALL > 5e-9)
## of ones not sig in pooled analysis, all are unique to one other group
table(rowSums(notSigAll < 5e-9, na.rm=TRUE))
sum(notSigAll$Score.pval.White < 5e-9, na.rm=TRUE) ## 111
sum(notSigAll$Score.pval.Black < 5e-9, na.rm=TRUE) ## 28
sum(notSigAll$Score.pval.HispanicLatino < 5e-9, na.rm=TRUE) ## 7
sum(notSigAll$Score.pval.Asian < 5e-9, na.rm=TRUE) ## 1
sum(notSigAll$Score.pval.Brazilian < 5e-9, na.rm=TRUE) ## 2
sum(notSigAll$Score.pval.Samoan < 5e-9, na.rm=TRUE) ## 0

idsByRound<-condRes %>% select(chr, grep("snpID_", colnames(condRes))) %>% pivot_longer(-chr, names_to="round", values_to ="snpID") %>% mutate(round = sub("snpID_", "", round))
pvalsByRound <- condRes %>% select(chr, grep("Score.pval_", colnames(condRes))) %>% pivot_longer(-chr, names_to="round", values_to ="Score.pval") %>% mutate(round = sub("Score.pval_", "", round))
posByRound<- condRes %>% select(chr, grep("pos_", colnames(condRes))) %>% pivot_longer(-chr, names_to="round", values_to ="pos") %>% mutate(round = sub("pos_", "", round))

allByRound <- idsByRound %>% left_join(pvalsByRound) %>% left_join(posByRound)
## remove NAs and positions with p>5e-9
allByRound <- allByRound %>% filter(!is.na(Score.pval), Score.pval < 5e-9) 

### OLD CODE ####
newPeaks %>% filter(!snpID %in% condRes$snpID)

allWithLocusGroup<-read.csv("~/Research/OneDrive/telomere/gwas/results/ALL_GWASResults/allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop_p_lt_5e-8_withLocusGroups.csv", header = TRUE, stringsAsFactors = FALSE)

## get start and stop of locus regions
allWithLocusGroup %>% filter(c(locusGroup[1:(nrow(allWithLocusGroup)-1)] != locusGroup[2:nrow(allWithLocusGroup)], TRUE) | c(TRUE, locusGroup[2:nrow(allWithLocusGroup)] != locusGroup[1:(nrow(allWithLocusGroup) - 1)])) %>% select(chr, pos)

condRes <- allWithLocusGroup %>% select(chr, pos, newLOCUS) %>% right_join(condRes, by=c("chr", "pos"))
condRes %>% select(newLOCUS, snpID, roundNum_Cond)
condRes[condRes$snpID == "4:163126692:T:G", "newLOCUS"]<-"NAF1"
condRes[condRes$snpID == "5:1292843:C:T", "newLOCUS"]<-"TERT"
condRes[condRes$snpID == "10:103915847:C:T", "newLOCUS"]<-"SH3PXD2A_OBFC1_SLK"
condRes[condRes$snpID == "18:676473:C:T", "newLOCUS"]<-"TYMS"
condRes[condRes$snpID == "18:650764:C:T", "newLOCUS"]<-"TYMS"
condRes[condRes$snpID == "22:40023952:C:T", "newLOCUS"]<-"TNRC6B"

## want to incorporate locus names
locusInfo <- read_excel("~/Research/OneDrive/telomere/gwas/results/forRasikaOASISLookup-withnovelty_v2.xlsx")
locusInfo <- locusInfo %>% mutate(novelty = ifelse(is.na(novelty), "Known", "Novel")) %>% select(snpID, newLOCUS = LocusName.Final, novelty)

all(locusInfo$snpID %in% condRes$snpID)
condRes <- condRes 

# select min p-value per locus
condRes %>% arrange(as.numeric(chr), newLOCUS, roundNum_Cond)
locusOrder <- cbind(condRes %>% group_by(newLOCUS) %>% filter(roundNum_Cond == min(roundNum_Cond)) %>% select(chr, snpID, newLOCUS, roundNum_Cond) %>% arrange(as.numeric(chr), roundNum_Cond) %>% select(newLOCUS), locusOrder = 1:length(unique(condRes$newLOCUS))) %>% print(n=Inf)

condRes <- condRes %>% left_join(locusOrder)
condRes <- condRes %>% arrange(as.numeric(chr), locusOrder, roundNum_Cond) 

forTable1 <- condRes %>% mutate(Est_Primary_bp = 1000 * Est_Primary, Est_Cond_bp = 1000 * Est_Cond, position=paste0("Chr", chr, ":", pos)) %>% select(LocusName = newLOCUS, position, ref, alt, freq, Score.pval_Primary, Est_Primary_bp, PVE_Primary, Score.pval_Cond, Est_Cond_bp, PVE_Cond, roundOfConditioning = roundNum_Cond)

write.csv(forTable1, file="~/Research/OneDrive/telomere/gwas/results/Table1.csv", row.names=FALSE, quote=TRUE, na="-")

write.csv(condRes %>% group_by(newLOCUS) %>% summarise(chr = unique(chr), minSNPPos = min(pos), maxSNPPos=max(pos), locusOrder = unique(locusOrder)) %>% arrange(locusOrder),  file="~/Research/OneDrive/telomere/gwas/results/lociStartStop.csv", row.names=FALSE, quote=TRUE, na="-")

annoInfo<-read.csv(gzfile("~/Research/OneDrive/telomere/gwas/results/MERGED_allChrs_telomere_adjagesexstudyseqctrbatchPCs_minDP0_BRAVODepthDrop_p_lt_5e-8_OASISAnno.csv.gz"), header = TRUE, stringsAsFactors = FALSE)
condRes <- condRes %>% left_join(annoInfo)

forRasika <- condRes %>% select(LocusName = newLOCUS, snpID, rsNum, Chr38, Pos38, Pos37, freq, Score.pval_Primary, Est_Primary, PVE_Primary, Score.pval_Cond, Est_Cond, PVE_Cond, roundOfConditioning = roundNum_Cond)
write.csv(forRasika,  file="~/Research/OneDrive/telomere/gwas/results/forRasikaOASISLookup.csv", row.names=FALSE, quote=TRUE, na="-")


lociRes<- read_excel("~/Research/OneDrive/telomere/gwas/results/forRasikaOASISLookup.xlsx")
locusDef<-lociRes %>% filter(!is.na(`"-Kb"`)) %>% mutate(Pos38 = as.numeric(Pos38), start = Pos38 - `"-Kb"`*1000, end = Pos38 + `"+Kb"`*1000) %>% select(locusFINAL = LocusName.Final, chr = Chr38, start, end)
write.csv(locusDef, file="~/Research/OneDrive/telomere/gwas/results/locusDefinitions.csv", quote=FALSE, row.names=FALSE)

