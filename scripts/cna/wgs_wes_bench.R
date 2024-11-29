#!/usr/bin/env Rscript

####################################################################################
# Script to benchmark Dragen CNAs call on WES and IchorCNA CNAs call on low-pass WGS
####################################################################################

args <- commandArgs(trailingOnly=TRUE)

set.seed(1225)

wes_dir = args[1] # Directory containing Dragen output (.tsv file)
wgs_res = args[2] # Output of low_pass_wgs_analysis.R analysis
output_dir = args[3] # Directory for outputs tables and graphs


suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(gridExtra)) # For arranging plots in a grid
suppressPackageStartupMessages(library(plyr))


if (!dir.exists(output_dir)) {
   dir.create(output_dir)
  }

# Output files
output_bench=paste0(output_dir, '/CNAs_wes_wgs.txt')

res<-fread(wgs_res)[,1:14]

con <- list.files(path = wes_dir,pattern = '.tsv',full.names=T)

#WES CNAs calls with Dragen
owes<-data.frame()
for (i in con) {
wes=as.data.frame(fread(i))
wes$Sample<-gsub('.tsv','',i)
owes<-rbind(owes,wes)
}

#Matching with low pass WGS calls with IchorCNA
data_wes_E<-data.frame()
for (i in con) {
#Low WGS Tumor
data_wes<-res[res$Sample==gsub('.*/|.tsv','',i),]
tmp <- GRanges(data_wes$Chr, ranges =IRanges(as.numeric(data_wes$Start),as.numeric(data_wes$End)))
#WES Tumor
wes=as.data.frame(fread(i))
g1<-GRanges(wes$Chromosome, ranges =IRanges(wes$Start,wes$End))

type3 = findOverlaps(query = tmp, subject = g1)
wes_m = data.frame(data_wes[queryHits(type3),], wes[subjectHits(type3),2:7])
wes_m$Sample<-gsub('T1-T.cnv.annotated.tsv','',i)
data_wes_E<-rbind(data_wes_E,wes_m)
}
data_wes_E$Tumor_WES_Tumor_WGS_match<-ifelse(data_wes_E$Type=='DUP' & data_wes_E$Tumor.Call=='High-level Amplification',TRUE,ifelse(data_wes_E$Type=='DUP' & data_wes_E$Tumor.Call=='Amplification',TRUE,
ifelse(data_wes_E$Type=='DUP' & data_wes_E$Tumor.Call=='Copy Gain',TRUE,ifelse(data_wes_E$Type=='DEL' & data_wes_E$Tumor.Call=='Homozygous Deletions State',TRUE,
ifelse(data_wes_E$Type=='DEL' & data_wes_E$Tumor.Call=='Hemizygous Deletions',TRUE,FALSE)))))

data_wes_E$Tumor_WES_Plasma_WGS_match<-ifelse(data_wes_E$Type=='DUP' & data_wes_E$Plasma.Call=='High-level Amplification',TRUE,ifelse(data_wes_E$Type=='DUP' & data_wes_E$Plasma.Call=='Amplification',TRUE,
ifelse(data_wes_E$Type=='DUP' & data_wes_E$Plasma.Call=='Copy Gain',TRUE,ifelse(data_wes_E$Type=='DEL' & data_wes_E$Plasma.Call=='Homozygous Deletions State',TRUE,
ifelse(data_wes_E$Type=='DEL' & data_wes_E$Plasma.Call=='Hemizygous Deletions',TRUE,FALSE)))))
colnames(data_wes_E)[c(2:4,15:17)]<-c('chr_wgs','start_wgs','end_wgs','chr_wes','start_wes','end_wes')
data_wes_E$Sample<-gsub('.tsv','',data_wes_E$Sample)
data_wes_E$id_wgs<-paste0(data_wes_E$Sample,data_wes_E$chr_wgs,data_wes_E$start_wgs,data_wes_E$end_wgs)
data_wes_E$id_wes<-paste0(data_wes_E$Sample,data_wes_E$chr_wes,data_wes_E$start_wes,data_wes_E$end_wes)
print(length(unique(data_wes_E$id_wgs)))
print(length(unique(data_wes_E$id_wes)))

#Filter by keeping only regions that overlap by 50%

gr_wgs <- GRanges(seqnames = data_wes_E$chr_wgs,
                  ranges = IRanges(start = as.numeric(data_wes_E$start_wgs), end = as.numeric(data_wes_E$end_wgs)))
gr_wes <- GRanges(seqnames = data_wes_E$chr_wes,
                  ranges = IRanges(start = as.numeric(data_wes_E$start_wes), end = as.numeric(data_wes_E$end_wes)))
# Find overlaps
overlaps <- findOverlaps(gr_wgs, gr_wes)
# Calculate percentage overlap in a vectorized way
query_ranges <- gr_wgs[queryHits(overlaps)]
subject_ranges <- gr_wes[subjectHits(overlaps)]
intersection_ranges <- pintersect(query_ranges, subject_ranges)

percent_query <- width(intersection_ranges) / width(query_ranges)
percent_subject <- width(intersection_ranges) / width(subject_ranges)

# Combine results into a data frame
overlap_results <- data.frame(
  query_index = queryHits(overlaps),
  subject_index = subjectHits(overlaps),
  percent_query = percent_query,
  percent_subject = percent_subject
)

# Filter for overlaps with at least 50% overlap in either range
filtered_results <- overlap_results[
  overlap_results$percent_query >= 0.5, ]
filtered_results<-filtered_results[!duplicated(filtered_results$query_index),]
data_wes_E50<-data_wes_E[filtered_results$query_index,]

# No filter concordant CNAs
data_wes_E50_ALL<-data_wes_E50[!duplicated(data_wes_E50$id_wes),]
# Filter only concordant CNAs
data_wes_E50<-data_wes_E50[!duplicated(data_wes_E50$id_wes) & data_wes_E50$Tumor_WES_Tumor_WGS_match== TRUE,]
#print(nrow(data_wes_E50)/nrow(owes))

results<-data.frame(Sample=unique(owes$Sample),
               Match=apply(as.array(unique(owes$Sample)),1,function(x) sum(data_wes_E50[['Sample']]==x)),
               Match_perc=apply(as.array(unique(owes$Sample)),1,function(x) round((sum(data_wes_E50[['Sample']]==x)/ sum(owes[['Sample']]==x))*100,digits=2)))

print(results)
fwrite(results,file=output_bench, row.names=F, col.names=T,quote=F, sep='\t')

