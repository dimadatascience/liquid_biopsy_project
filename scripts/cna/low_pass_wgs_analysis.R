#!/usr/bin/env Rscript

#################################################################################
# Script to compare plasma and tumor CNAs call as output from DeSeq-Free pipeline
# It also generates stats and graphs for results evaluation
#################################################################################

args <- commandArgs(trailingOnly=TRUE)

set.seed(1225)

input_dir = args[1] # Directory containing IchorCNA output (.seg file) from https://github.com/mdelcorvo/DeSeq-Free pipeline
output_dir = args[2] # Directory for outputs tables and graphs
annot = args[3] # Gene annotation (.bed)

suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(gridExtra)) # For arranging plots in a grid
suppressPackageStartupMessages(library(plyr))

##########
#Functions
##########

# Function to compute correlation and p-value
compute_correlation_and_pvalue <- function(df) {
  if (nrow(df) > 1) { # Ensure there are at least two points to calculate correlation
    test <- cor.test(df$Tumor.LogR, df$Plasma.LogR, method = "pearson")
    return(data.frame(correlation = test$estimate, p_value = test$p.value))
  } else {
    return(data.frame(correlation = NA, p_value = NA)) # Return NA if not computable
  }
}

calculate_counts <- function(data, column, neutral_value = "Copy Neutral") {
  if (nrow(data) == 0) return(0)
  return(sum(data[[column]] != neutral_value))
}

##########

#setwd(input_dir)
if (!dir.exists(output_dir)) {
   dir.create(output_dir)
  }

# Output files
output_CNAs_comp=paste0(output_dir, '/CNAs_calls.txt')
output_plasma_tumor_corr=paste0(output_dir, '/plasma_tumor_corr.tiff')
output_plasma_tumor_sample_corr=paste0(output_dir, '/plasma_tumor_sample_corr.tiff')
output_CNAs_perc=paste0(output_dir, '/CNAs_perc.txt')
##############

ann<-fread(annot) #'canonical_exon_transcripts_hg38.bed'
tmp <- GRanges(ann$V1, ranges =IRanges(ann$V2,ann$V3))

###################################
# Low WGS plasma + tumor comparison
###################################

con <- list.files(path = input_dir,pattern = '.cna.seg',full.names=T)
data_all<-data.frame()
for (i in unique(gsub('-plasma.*|-tumor.*','',con))) {

#Plasma seg
plasma=paste0(i,'-plasma.cna.seg')
cna_p<-as.data.frame(fread(plasma))
cna_p <- cna_p[,c(1:3,6,8:10)]
cna_p$chr<-paste0('chr',cna_p$chr)
g1<-GRanges(cna_p$chr, ranges =IRanges(cna_p$start,cna_p$end))
type3 = findOverlaps(query = tmp, subject = g1)

cap = data.frame(cna_p[subjectHits(type3),], ann[queryHits(type3),])
cap<-cap[,c(1:7,11)]
colnames(cap)<-c('Chr','Start','End','Plasma.LogR','Plasma.Copy_Number','Plasma.Call','Plasma.LogR_Copy_Number','Gene_name')
cap1<-cap[!is.na(cap$Plasma.LogR_Copy_Number),]
cap1$pos<- paste(cap1$Chr,cap1$Start,cap1$End,sep='-')
cap2<- as.data.frame(group_by(cap1, pos) %>%
  summarise_all(funs(paste(unique(.), collapse = "|"))))
cap2$Plasma.Call <- ifelse(cap2$Plasma.Call=='NEUT','Copy Neutral', ifelse(cap2$Plasma.Call=='HETD','Hemizygous Deletions',
					ifelse(cap2$Plasma.Call=='GAIN','Copy Gain',ifelse(cap2$Plasma.Call=='HLAMP','High-level Amplification',
					ifelse(cap2$Plasma.Call=='AMP','Amplification',ifelse(cap2$Plasma.Call=='HOMD','Homozygous Deletions State',cap2$Plasma.Call))))))

#Tumor seg
tumor=paste0(i,'-tumor.cna.seg')
cat<-as.data.frame(fread(tumor))
cat <- cat[,c(1:3,6,8:10)]
cat$chr<-paste0('chr',cat$chr)
colnames(cat)<-c('Chr','Start','End','Tumor.LogR','Tumor.Copy_Number','Tumor.Call','Tumor.LogR_Copy_Number')
cat1<-cat[!is.na(cat$Tumor.LogR_Copy_Number),]
cat1$pos<- paste(cat1$Chr,cat1$Start,cat1$End,sep='-')
cat2<- as.data.frame(group_by(cat1, pos) %>%
  summarise_all(funs(paste(unique(.), collapse = "|"))))
cat2$Tumor.Call <- ifelse(cat2$Tumor.Call=='NEUT','Copy Neutral', ifelse(cat2$Tumor.Call=='HETD','Hemizygous Deletions',
					ifelse(cat2$Tumor.Call=='GAIN','Copy Gain',ifelse(cat2$Tumor.Call=='HLAMP','High-level Amplification',
					ifelse(cat2$Tumor.Call=='AMP','Amplification',ifelse(cat2$Tumor.Call=='HOMD','Homozygous Deletions State',cat2$Tumor.Call))))))

#Comparison
ca<-merge(cat2,cap2,by.x='pos',by.y='pos')
ca$pos<-NULL
ca<-ca[,c(1:3,15,4:7,11:14)]
colnames(ca)[1:3]<-c('Chr','Start','End')
ca$Sample<-i
data_all<-rbind(data_all,ca)
}

res <- data_all[with(data_all, order(data_all$Sample, data_all$Tumor.LogR,decreasing=T)), ]
res<-res[,c(ncol(res),1:(ncol(res)-1))]
res$Tumor_Plamsa_match<-ifelse(res$Tumor.Call==res$Plasma.Call,TRUE,FALSE)

res$Tumor.LogR <- as.numeric(res$Tumor.LogR)
res$Plasma.LogR <- as.numeric(res$Plasma.LogR)

fwrite(res,file=output_CNAs_comp, row.names=F, col.names=T,quote=F, sep='\t')

##################
#Correlation graph
##################

filtered_res <- subset(res, Tumor_Plamsa_match == TRUE)

plot_original <- ggscatter(res, x = "Tumor.LogR", y = "Plasma.LogR",
   add = "reg.line",  # Add regressin line
   add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
   conf.int = TRUE # Add confidence interval
   ) + stat_cor(method = "pearson") +
  labs(
    title = paste(
      "Tumor LogR vs Plasma LogR (Original Data)\n"
    ),
    x = "Tumor LogR",
    y = "Plasma LogR"
  ) +
  theme_minimal()

plot_filtered <- ggscatter(filtered_res, x = "Tumor.LogR", y = "Plasma.LogR",
   add = "reg.line",  # Add regressin line
   add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
   conf.int = TRUE # Add confidence interval
   ) + stat_cor(method = "pearson") +
  labs(
    title = paste(
      "Tumor LogR vs Plasma LogR (matched CNA Data)\n"
    ),
    x = "Tumor LogR",
    y = "Plasma LogR"
  ) +
  theme_minimal()

tiff(output_plasma_tumor_corr, width = 9, height = 7, units = 'in', res = 300, compression = 'lzw')
grid.arrange(plot_original, plot_filtered, ncol = 1)
dev.off()


# Calculate correlation and p-value for original data
cor_original <- res %>%
  group_by(Sample) %>%
  group_modify(~ compute_correlation_and_pvalue(.x)) %>%
  mutate(type = "Original")

# Filter data where Tumor.Plamsa.CN.Match is TRUE and calculate correlation and p-value
cor_filtered <- res %>%
  filter(Tumor_Plamsa_match == TRUE) %>%
  group_by(Sample) %>%
  group_modify(~ compute_correlation_and_pvalue(.x)) %>%
  mutate(type = "Only matched CNA")

# Combine the two datasets
cor_combined <- bind_rows(cor_original, cor_filtered)
# Remove rows with NA correlations to avoid plotting issues
cor_combined <- cor_combined %>% drop_na(correlation)

# Create the grouped bar plot
p<- ggplot(cor_combined, aes(x = Sample, y = correlation, fill = type)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = paste0("p = ", signif(p_value, 2))),
            position = position_dodge(width = 0.9), vjust = -0.3, size = 3) +
  labs(
    title = "Correlation per Sample: Original vs Only matched CNA Data",
    x = "Sample",
    y = "Correlation",
    fill = "Data Type"
  ) +
  scale_y_continuous(
    breaks = seq(-1, 1, by = 0.1) # Set y-axis breaks to show all values from -1 to 1
  ) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(face = "bold", size = 12), # Bold x-axis label
    axis.title.y = element_text(face = "bold", size = 12), # Bold y-axis label
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5), # Bold plot title
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 10), # Bold x-axis text (sample names)
    legend.title = element_text(face = "bold"), # Bold legend title
    legend.text = element_text(size = 10) # Adjust legend text size
  )

tiff(output_plasma_tumor_sample_corr, width = 9, height = 7, units = 'in', res = 300, compression = 'lzw')
p
dev.off()

###############################################
# Get % of CNAs that differ from 'Copy Neutral'
###############################################

# Calculate results for each sample
results <- res %>%
  group_by(Sample) %>%
  group_modify(~ {
    total_rows <- nrow(.x)
    total_tumor_diff <- calculate_counts(.x, "Tumor.Call")
    total_plasma_diff <- calculate_counts(.x, "Plasma.Call")
    matched_data <- .x[.x$Tumor_Plamsa_match == TRUE, ]

    matched_tumor_diff <- calculate_counts(matched_data, "Tumor.Call")
    matched_plasma_diff <- calculate_counts(matched_data, "Plasma.Call")

    tibble(
      Tumor_Diff_Count = total_tumor_diff,
      Tumor_Diff_Percent = ifelse(total_rows > 0, (total_tumor_diff / total_rows) * 100, 0),
      Plasma_Diff_Count = total_plasma_diff,
      Plasma_Diff_Percent = ifelse(total_rows > 0, (total_plasma_diff / total_rows) * 100, 0),
      Tumor_PLasma_Diff_Count_Matched = matched_tumor_diff,
      Tumor_Diff_Percent_Matched = ifelse(total_tumor_diff > 0,
                                          (matched_tumor_diff / total_tumor_diff) * 100, 0),
      Plasma_Diff_Percent_Matched = ifelse(total_plasma_diff > 0,
                                           (matched_plasma_diff / total_plasma_diff) * 100, 0)
    )
  }) %>%
  ungroup()

# Display the results
print(results)
fwrite(results,file=output_CNAs_perc, row.names=F, col.names=T,quote=F, sep='\t')

