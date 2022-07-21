
# P. Murat, MRC Laboratory of Molecular Biology, July 2022
# Code for generating Figures 5 and 6 with associated supporting Figures
# in Murat et al. DNA replication initiation shapes the mutational landscape and expression of the human genome

library(dplyr)
library(biomaRt)
library(ggpubr)
library(TFBSTools)
library(ggseqlogo)
library(scales)
library(zoo)
library(Matrix)
library(readr)
library(readxl)
library(stringr)

setwd("/Volumes/pmurat/OriMut")

############################################################################
############################################################################
###               Association between functional SNPs                    ###
###                   and constitutive origins                           ###
############################################################################
############################################################################

############################################################################
# Load and format datasets reporting functional SNPs 

# SNP functional annotation is from RegulomeDB
# RegulonmeDB is a database that annotates SNPs with known and predicted regulatory elements in the H. sapiens genome
# Known and predicted regulatory DNA elements include regions of DNase hypersensitivity, binding sites of transcription factors,
# and promoter regions that have been biochemically characterized to regulation transcription.
# Sources of these data include public datasets from GEO, the ENCODE project, and published literature.
# File: TSTFF344324.tsv
# Downloaded from https://regulomedb.org/regulome-help/

regulome <- read.csv("./Dataset/TSTFF344324.tsv", sep = '\t')
# refer to https://regulomedb.org/regulome-help/ for the definition of rank
# are considered functional only SNPs with a rank of 2 or higher (eQTL or found in ENCODE ChIP-seq signals)
# assign a simplified rank system
regulome.score <- regulome %>% mutate(SNP.function = case_when(ranking %in% c("1a", "1b", "1c", "1d", "1e", "1f", "2a", "2b", "2c") ~ "functional",
                                                               TRUE ~ "non functional"))

# Prepare bed files of sequence variants for closest with SNP ids

population.SNP <- read.table("./Dataset/REDUCED_common_all_20180418.SNP.recode.vcf", header = F)
colnames(population.SNP) <- c("chr", "pos", "SNP.id", "REF", "ALT")
nrow(population.SNP) # 33,629,539 SNPs

population.SNP.bed <- population.SNP %>% mutate(start = pos, end = as.integer(pos+1)) %>% dplyr::select(chr, start, end, SNP.id, REF, ALT)
write.table(population.SNP.bed, "./Figure_5/population.SNP.bed", sep="\t", col.names = F, row.names = F, quote = F)
# using zsh interactive shell on Mac
# sort -k1,1 -k2,2n ./Figure_5/population.SNP.bed > ./Figure_4/population.SNP.sorted.bed

############################################################################
# Assess the enrichment of functional SNPs around origins

# Compute distance of all SNPs to closest origins, gene starts and splice sites

# sort -k1,1 -k2,2n ./Figure_4/IniSeq.ori.center.bed > ./Figure_4/IniSeq.ori.center.sorted.bed
# sort -k1,1 -k2,2n ./Figure_4/SNS.core.ori.center.bed > ./Figure_4/SNS.core.ori.center.sorted.bed
# sort -k1,1 -k2,2n ./Figure_4/SNS.stochastic.ori.center.bed > ./Figure_4/SNS.stochastic.ori.center.sorted.bed

# bedtools closest -a ./Figure_5/population.SNP.sorted.bed -b ./Figure_4/hsapiens.genes.start.sorted.bed > ./Figure_5/population.SNP.closest.gene.start.bed
# bedtools closest -a ./Figure_5/population.SNP.sorted.bed -b ./Figure_4/hsapiens.exon.intron.junctions.sorted.bed > ./Figure_5/population.SNP.closest.exon.intron.junctions.bed
# bedtools closest -a ./Figure_5/population.SNP.sorted.bed -b ./Figure_4/hsapiens.intron.exon.junctions.sorted.bed > ./Figure_5/population.SNP.closest.intron.exon.junctions.bed
# bedtools closest -a ./Figure_5/population.SNP.sorted.bed -b ./Figure_4/IniSeq.ori.center.sorted.bed > ./Figure_5/population.SNP.closest.IniSeq.ori.bed
# bedtools closest -a ./Figure_5/population.SNP.sorted.bed -b ./Figure_4/SNS.core.ori.center.sorted.bed > ./Figure_5/population.SNP.closest.SNS.core.ori.bed
# bedtools closest -a ./Figure_5/population.SNP.sorted.bed -b ./Figure_4/SNS.stochastic.ori.center.sorted.bed > ./Figure_5/population.SNP.closest.SNS.stochastic.ori.bed

population.SNP.closest.IniSeq.ori <- read.table("./Figure_5/population.SNP.closest.IniSeq.ori.bed", header = F) %>% dplyr::select(V1, V2, V4, V5, V6, V8, V9, V10)
colnames(population.SNP.closest.IniSeq.ori) <- c("chr", "SNP.pos", "SNP.id", "REF", "ALT", "ori.start", "ori.end", "ori.eff")
population.SNP.closest.IniSeq.ori.10kb <- population.SNP.closest.IniSeq.ori %>% mutate(SNP.ori.dist = SNP.pos-round((ori.start+ori.end)/2)) %>% filter(SNP.ori.dist >= -10000 & SNP.ori.dist <= 10000)
nrow(population.SNP.closest.IniSeq.ori.10kb) # 4,119,879 SNPs
# Save
write.table(population.SNP.closest.IniSeq.ori.10kb, "./Figure_5/population.SNP.closest.IniSeq.ori.10kb.bed", quote=FALSE, sep='\t', col.names = T)

population.SNP.closest.SNS.core.ori <- read.table("./Figure_5/population.SNP.closest.SNS.core.ori.bed", header = F) %>% dplyr::select(V1, V2, V4, V5, V6, V8, V9, V10)
colnames(population.SNP.closest.SNS.core.ori) <- c("chr", "SNP.pos", "SNP.id", "REF", "ALT", "ori.start", "ori.end", "ori.eff")
population.SNP.closest.SNS.core.ori.10kb <- population.SNP.closest.SNS.core.ori %>% mutate(SNP.ori.dist = SNP.pos-round((ori.start+ori.end)/2)) %>% filter(SNP.ori.dist >= -10000 & SNP.ori.dist <= 10000)
nrow(population.SNP.closest.SNS.core.ori.10kb) # 4,119,443 SNPs
# Save
write.table(population.SNP.closest.SNS.core.ori.10kb, "./Figure_5/population.SNP.closest.SNS.core.ori.10kb.bed", quote=FALSE, sep='\t', col.names = T)

population.SNP.closest.SNS.stochastic.ori <- read.table("./Figure_5/population.SNP.closest.SNS.stochastic.ori.bed", header = F) %>% dplyr::select(V1, V2, V4, V5, V6, V8, V9, V10)
colnames(population.SNP.closest.SNS.stochastic.ori) <- c("chr", "SNP.pos", "SNP.id", "REF", "ALT", "ori.start", "ori.end", "ori.eff")
population.SNP.closest.SNS.stochastic.ori.10kb <- population.SNP.closest.SNS.stochastic.ori %>% mutate(SNP.ori.dist = SNP.pos-round((ori.start+ori.end)/2)) %>% filter(SNP.ori.dist >= -10000 & SNP.ori.dist <= 10000)
nrow(population.SNP.closest.SNS.stochastic.ori.10kb) # 4,119,443 SNPs
# Save
write.table(population.SNP.closest.SNS.stochastic.ori.10kb, "./Figure_5/population.SNP.closest.SNS.stochastic.ori.10kb.bed", quote=FALSE, sep='\t', col.names = T)

# Load information for population SNP distribution around origins

population.SNP.closest.IniSeq.ori.10kb <- read.csv("./Figure_5/population.SNP.closest.IniSeq.ori.10kb.bed", sep = '\t')
population.SNP.closest.SNS.core.ori.10kb <- read.csv("./Figure_5/population.SNP.closest.SNS.core.ori.10kb.bed", sep = '\t')
population.SNP.closest.SNS.stochastic.ori.10kb <- read.csv("./Figure_5/population.SNP.closest.SNS.stochastic.ori.10kb.bed", sep = '\t')

# Count the percentage of functional SNPs

# All SNPs
table(regulome.score$SNP.function)
nrow(regulome.score) # 13275023 SNPs
# functional   non functional 
# 524563       12750460 
# 4.114071     95.88593%

# SNPs at IniSeq origins (+- 2.5kb)
population.SNP.closest.IniSeq.ori.2.5kb <- population.SNP.closest.IniSeq.ori.10kb %>% filter(SNP.ori.dist >= -2500 & SNP.ori.dist <= 2500)
regulome.score.IniSeq.ori <- regulome.score %>% filter(rsid %in% population.SNP.closest.IniSeq.ori.2.5kb$SNP.id) # 480,757 SNPs
table(regulome.score.IniSeq.ori$SNP.function)
# functional non functional 
#  92553         388204 
# 19.25151     80.74849%
# 2x2 chi-squared significance test (comparing IniSeq SNPs and all SNPs)
df1 <- regulome.score %>% dplyr::select(SNP.function) %>% mutate(condition = "all")
df2 <- regulome.score.IniSeq.ori %>% dplyr::select(SNP.function) %>% mutate(condition = "ori")
df.test <- rbind(df1, df2)
pval <- chisq.test(df.test$SNP.function, df.test$condition) # p-value < 2.2e-16, enr = 4.679431

# SNPs at SNS core origins (+- 2.5kb)
population.SNP.closest.SNS.core.ori.2.5kb <- population.SNP.closest.SNS.core.ori.10kb %>% filter(SNP.ori.dist >= -2500 & SNP.ori.dist <= 2500)
regulome.score.SNS.core.ori <- regulome.score %>% filter(rsid %in% population.SNP.closest.SNS.core.ori.2.5kb$SNP.id) # 609,716 SNPs
table(regulome.score.SNS.core.ori$SNP.function)
# functional non functional 
#  38675         571204 
# 6.343117      93.68362%
# 2x2 chi-squared significance test (comparing SNS core SNPs and all SNPs)
df1 <- regulome.score %>% dplyr::select(SNP.function) %>% mutate(condition = "all")
df2 <- regulome.score.SNS.core.ori %>% dplyr::select(SNP.function) %>% mutate(condition = "ori")
df.test <- rbind(df1, df2)
pval <- chisq.test(df.test$SNP.function, df.test$condition) # p-value < 2.2e-16, enr = 1.54181

# SNPs at SNS stochastic origins (+- 2.5kb)
population.SNP.closest.SNS.stochastic.ori.2.5kb <- population.SNP.closest.SNS.stochastic.ori.10kb %>% filter(SNP.ori.dist >= -2500 & SNP.ori.dist <= 2500)
regulome.score.SNS.stochastic.ori <- regulome.score %>% filter(rsid %in% population.SNP.closest.SNS.stochastic.ori.2.5kb$SNP.id) # 1,236,126 SNPs
table(regulome.score.SNS.stochastic.ori$SNP.function)
# functional non functional 
#  66908         1169218 
# 5.412717      94.58728%
# 2x2 chi-squared significance test (comparing SNS core SNPs and all SNPs)
df1 <- regulome.score %>% dplyr::select(SNP.function) %>% mutate(condition = "all")
df2 <- regulome.score.SNS.stochastic.ori %>% dplyr::select(SNP.function) %>% mutate(condition = "ori")
df.test <- rbind(df1, df2)
pval <- chisq.test(df.test$SNP.function, df.test$condition) # p-value < 2.2e-16, enr = 1.31566

# SNPs at low efficiency origins (+- 2.5kb)
population.SNP.closest.IniSeq.ori.low.2.5kb <- population.SNP.closest.IniSeq.ori.2.5kb %>% filter(ori.eff < 0.8418) %>% filter(SNP.ori.dist >= -2500 & SNP.ori.dist <= 2500)
regulome.score.IniSeq.ori.low <- regulome.score %>% filter(rsid %in% population.SNP.closest.IniSeq.ori.low.2.5kb$SNP.id) # 172,094 SNPs
table(regulome.score.IniSeq.ori.low$SNP.function)
#  functional non functional 
#    29908         142186 
#    17.37887     82.62113%

# SNPs at medium efficiency origins (+- 2.5kb)
population.SNP.closest.IniSeq.ori.medium.2.5kb <- population.SNP.closest.IniSeq.ori.2.5kb %>% filter(ori.eff < 0.9065 & ori.eff >= 0.8418) %>% filter(SNP.ori.dist >= -2500 & SNP.ori.dist <= 2500)
regulome.score.IniSeq.ori.medium <- regulome.score %>% filter(rsid %in% population.SNP.closest.IniSeq.ori.medium.2.5kb$SNP.id) # 159,850 SNPs
table(regulome.score.IniSeq.ori.medium$SNP.function)
#  functional non functional 
#    30425         129425 
#    19.03347     80.966531%

# SNPs at high efficiency origins (+- 2.5kb)
population.SNP.closest.IniSeq.ori.high.2.5kb <- population.SNP.closest.IniSeq.ori.2.5kb %>% filter(ori.eff >= 0.9065) %>% filter(SNP.ori.dist >= -2500 & SNP.ori.dist <= 2500)
regulome.score.IniSeq.ori.high <- regulome.score %>% filter(rsid %in% population.SNP.closest.IniSeq.ori.high.2.5kb$SNP.id) # 148,829 SNPs
table(regulome.score.IniSeq.ori.high$SNP.function)
#  functional non functional 
#    32224         116605 
#    21.65169     78.34831%

# chi-squared significance test to assess the contribution of origin efficiency to SNP function
df1 <- regulome.score.IniSeq.ori.low %>% dplyr::select(SNP.function) %>% mutate(condition = "low")
df2 <- regulome.score.IniSeq.ori.medium %>% dplyr::select(SNP.function) %>% mutate(condition = "medium")
df3 <- regulome.score.IniSeq.ori.high %>% dplyr::select(SNP.function) %>% mutate(condition = "high")
df.test <- rbind(df1, df2, df3)
pval <- chisq.test(df.test$SNP.function, df.test$condition) # p-value = 7.490841e-206

# Prepare plots

ori.df1 <- tibble(Percentage = 4.114071, class = "All SNPs")
ori.df2 <- tibble(Percentage = 19.25151, class = "IniSeq SNPs")
ori.df3 <- tibble(Percentage = 6.343117, class = "SNS core SNPs")
ori.df4 <- tibble(Percentage = 5.412717, class = "SNS stochastic SNPs")
ori.df <- rbind(ori.df1, ori.df2, ori.df3, ori.df4)

ori.plot <- ori.df %>% mutate(Class = fct_relevel(class, "All SNPs", "SNS stochastic SNPs", "SNS core SNPs","IniSeq SNPs")) %>%
  ggplot(aes(x=Class, y=Percentage, fill=Class)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=c("#D4D2D2", "#56B4E9", "#0072B2", "#DE1C14")) +
  ggtitle("all pval < 2.2e-16\nenr = 1.31, 1.54, 4.68") + ylab("Percentage of functional SNPs") +
  theme_bw() + theme(aspect.ratio=2, plot.title = element_text(size = 9)) + theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 0.5)) + theme(axis.title.x = element_blank()) +
  annotate("text", x = 2, y = 8, label = "***", size = 3) + annotate("text", x = 3, y = 8, label = "***", size = 3) + annotate("text", x = 4, y = 21, label = "***", size = 3)

eff.df1 <- tibble(Percentage = 17.37887, class = "low")
eff.df2 <- tibble(Percentage = 19.03347, class = "medium")
eff.df3 <- tibble(Percentage = 21.65169, class = "high")
eff.df <- rbind(eff.df1, eff.df2, eff.df3)

eff.plot <- eff.df %>% mutate(Class = fct_relevel(class, "low", "medium", "high")) %>%
  ggplot(aes(x=Class, y=Percentage, fill=Class)) +
  geom_bar(stat="identity") + ylim(0,30) +
  scale_fill_manual(values=c("#E69F00", "#D55E00", "#DE1C14")) +
  ggtitle("pval(assoc) = 7.490841e-206") + ylab("Percentage of functional SNPs") +
  theme_bw() + theme(aspect.ratio=2, plot.title = element_text(size = 9)) + theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 0.5)) + theme(axis.title.x = element_blank())

pdf("./Rplot/Fig_5A.pdf", width=4, height=4, useDingbats=FALSE)
ori.plot
dev.off()

pdf("./Rplot/Figure_S5A.pdf", width=4, height=4, useDingbats=FALSE)
eff.plot
dev.off()

############################################################################
############################################################################
###               Enrichment of functional SNPs at the                   ###
###                start of genes marked with origins                    ###
############################################################################
############################################################################

############################################################################
# Load data
population.SNP.closest.gene <- read.table("./Figure_5/population.SNP.closest.gene.start.bed", header = F) %>% dplyr::select(V1, V2, V4, V5, V6, V8, V10, V11, V12)
colnames(population.SNP.closest.gene) <- c("chr", "SNP.pos", "SNP.id", "REF", "ALT", "gene.start", "gene.strand", "gene.id", "gene.name")
# Select only SNP in 10 kb around gene start
population.SNP.closest.gene.10kb <- population.SNP.closest.gene %>% mutate(SNP.gene.dist = SNP.pos-gene.start) %>% filter(SNP.gene.dist >= -10000 & SNP.gene.dist <= 10000)
nrow(population.SNP.closest.gene.10kb) # 36,466,299 SNPs
# Save
write.table(population.SNP.closest.gene.10kb, "./Figure_4/population.SNP.closest.gene.start.10kb.bed", quote=FALSE, sep='\t', col.names = T)
population.SNP.closest.gene.10kb  <- read.csv("./Figure_4/population.SNP.closest.gene.start.10kb.bed", sep = '\t')

############################################################################
# Recover functional SNP ID
SNP.functional <- regulome.score %>% filter(SNP.function == "functional")

# Select functional SNPs
functional.SNP.closest.gene.10kb <- population.SNP.closest.gene.10kb %>% filter(SNP.id %in% SNP.functional$rsid) %>% dplyr::select(-SNP.gene.dist) %>% 
  mutate(dist = case_when(gene.strand == 1 ~ (SNP.pos-gene.start), TRUE ~ (gene.start-SNP.pos)))

# Load results from Fig. 4
hsapiens.gene.start.closest.IniSeq.ori <- read.csv("./Figure_4/hsapiens.genes.start.closest.IniSeq.center.bed", header = F, sep = "\t") %>% dplyr::select(V1, V2, V4, V5, V6, V8, V10, V11)
colnames(hsapiens.gene.start.closest.IniSeq.ori) <- c("chr", "gene.start", "gene.strand", "gene.id", "gene.name", "ori.pos", "ori.eff", "dist")
hsapiens.gene.start.closest.IniSeq.ori <- hsapiens.gene.start.closest.IniSeq.ori %>% mutate(dist.stranded = case_when(gene.strand == 1 ~ dist, TRUE ~ -dist))
hsapiens.gene.start.closest.IniSeq.ori.10kb <- hsapiens.gene.start.closest.IniSeq.ori %>% filter(dist.stranded >= -10000 & dist.stranded <= 10000) %>%
  mutate(dist.stranded = ((as.integer(cut(dist.stranded, breaks = 4000)))-2000)*5) # per 5bp bin

# Genes marked by origins at TSS
genes.TSS.with.ori <- hsapiens.gene.start.closest.IniSeq.ori.10kb %>% filter(dist <= 1250)  # 11,410 genes
genes.TSS.without.ori <- hsapiens.gene.start.closest.IniSeq.ori %>% filter(!gene.id %in% genes.TSS.with.ori$gene.id)  # 7,925 genes

# Save gene lists
write.csv(genes.TSS.with.ori, "./Figure_5/genes.TSS.with.ori.bed", quote=FALSE)
write.csv(genes.TSS.without.ori, "./Figure_5/genes.TSS.without.ori.bed", quote=FALSE)

# Select functional SNPs associated with genes marked by an origin
functional.SNP.closest.gene.10kb.with.ori <- functional.SNP.closest.gene.10kb %>% filter(gene.id %in% unique(genes.TSS.with.ori$gene.id)) # 11,032 genes
functional.SNP.closest.gene.10kb.without.ori <- functional.SNP.closest.gene.10kb %>% filter(!gene.id %in% unique(genes.TSS.with.ori$gene.id))  # 6,384 genes

############################################################################
# Compute enrichment profile from background in windows of 5nt and averaged over 20 windows

functional.SNP.closest.gene.10kb.with.ori.dens <- functional.SNP.closest.gene.10kb.with.ori %>% mutate(dist = ((as.integer(cut(dist, breaks = 4000)))-2000)*5) %>% 
  group_by(dist) %>% summarise(count = n()) %>% mutate(enr = log2(count/mean(count[c(1:100,900:1000)])), enr.average=rollapply(enr,20,mean,align='right',fill=0), class = "with origin")

functional.SNP.closest.gene.10kb.without.ori.dens <- functional.SNP.closest.gene.10kb.without.ori %>% mutate(dist = ((as.integer(cut(dist, breaks = 4000)))-2000)*5) %>% 
  group_by(dist) %>% summarise(count = n()) %>% mutate(enr = log2(count/mean(count[c(1:100,900:1000)])), enr.average=rollapply(enr,20,mean,align='right',fill=0), class = "without origin")

TSS.plot <- rbind(functional.SNP.closest.gene.10kb.with.ori.dens, functional.SNP.closest.gene.10kb.without.ori.dens)

pdf("./Rplot/Fige_5B.pdf", width=6, height=4, useDingbats=FALSE)
TSS.plot %>% ggplot(aes(x=dist, y=enr.average, color=class)) +
  geom_line(size = 0.5) + xlim(-5000,5000) + ylim(-0.2,3) +
  scale_color_manual(values=c("#DE1C14", "#BDBDBD")) +
  xlab("Distance from gene start (bp)") + ylab("Functional SNP enrichment (log2 FC)") + ggtitle("Functional SNPs around gene starts") +
  theme_bw() + theme(aspect.ratio=1.1, plot.title = element_text(size=12))
dev.off()

############################################################################
# Quantification and statistics

# Select functional SNPs associated with genes marked by an origin binned by efficiency
genes.TSS.with.ori.low <- genes.TSS.with.ori %>% filter(ori.eff < 0.8418)  # 1,424 genes
genes.TSS.with.ori.medium <- genes.TSS.with.ori %>% filter(ori.eff >= 0.8418 & ori.eff < 0.9065)  # 3,849 genes
genes.TSS.with.ori.high <- genes.TSS.with.ori %>% filter(ori.eff >= 0.9065)  # 6,137 genes
functional.SNP.closest.gene.10kb.with.ori.low <- functional.SNP.closest.gene.10kb %>% filter(gene.id %in% unique(genes.TSS.with.ori.low$gene.id)) # 1,360 genes
functional.SNP.closest.gene.10kb.with.ori.medium <- functional.SNP.closest.gene.10kb %>% filter(gene.id %in% unique(genes.TSS.with.ori.medium$gene.id)) # 3,711 genes
functional.SNP.closest.gene.10kb.with.ori.high <- functional.SNP.closest.gene.10kb %>% filter(gene.id %in% unique(genes.TSS.with.ori.high$gene.id)) # 5,961 genes

# Compute the number of functional SNPs in the vicinity of genes (+- 2.5kb around gene start) marked by an origins

functional.SNP.closest.gene.start.with.ori <- functional.SNP.closest.gene.10kb.with.ori %>% filter(dist >= -2500 & dist <= 2500)  # 41,670 mutations for 10,327 genes
functional.SNP.closest.gene.start.without.ori <- functional.SNP.closest.gene.10kb.without.ori %>% filter(dist >= -2500 & dist <= 2500)  # 15,533 mutations for 5,086 genes
functional.SNP.closest.gene.start.with.ori.low <- functional.SNP.closest.gene.10kb.with.ori.low %>% filter(dist >= -2500 & dist <= 2500)  # 4,876 mutations for 1,243 genes
functional.SNP.closest.gene.start.with.ori.medium <- functional.SNP.closest.gene.10kb.with.ori.medium %>% filter(dist >= -2500 & dist <= 2500)  # 13,780 mutations for 3,460 genes
functional.SNP.closest.gene.start.with.ori.high <- functional.SNP.closest.gene.10kb.with.ori.high %>% filter(dist >= -2500 & dist <= 2500)  # 23,014 mutations for 5,624 genes

# Prepare count tables

count.df1 <- tibble(Count = 41670/(10327*5000), class = "with ori")
count.df2 <- tibble(Count = 15533/(5086*5000), class = "without ori")
count.df <- rbind(count.df1, count.df2)
# 2x2 chi-squared significance test (comparing genes with or without ori)
pval.table <- rbind(c(41670,(10327*5000)-41670),c(15533,(5086*5000)-15533))
pval <- chisq.test(pval.table) # p-value = 2.964766e-194, enr = 1.321205

count.plot <- count.df %>% mutate(Class = fct_relevel(class, "without ori", "with ori")) %>%
  ggplot(aes(x=Class, y=Count, fill=Class)) +
  geom_bar(stat="identity") + ylim(0,0.001) +
  scale_fill_manual(values=c("#D4D2D2", "#DE1C14")) +
  ggtitle("pval = 2.964766e-194\nenr = 1.321205") + ylab("Number of functional SNPs per bp\n(TSS +- 2.5kb)") +
  theme_bw() + theme(aspect.ratio=2, plot.title = element_text(size = 11)) + theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 0.5)) + theme(axis.title.x = element_blank())

count.eff.df1 <- tibble(Count = 4876/(1243*5000), class = "low")
count.eff.df2 <- tibble(Count = 13780/(3460*5000), class = "medium")
count.eff.df3 <- tibble(Count = 23014/(5624*5000), class = "high")
count.eff.df <- rbind(count.eff.df1, count.eff.df2, count.eff.df3)
# 2x2 chi-squared significance test (comparing genes low with high)
pval.table <- rbind(c(4876,(1243*5000)-4876),c(23014,(5624*5000)-23014))
pval <- chisq.test(pval.table) # p-value = 0.007489555, enr = 1.043168

count.eff.plot <- count.eff.df %>% mutate(Class = fct_relevel(class, "low", "medium", "high")) %>%
  ggplot(aes(x=Class, y=Count, fill=Class)) +
  geom_bar(stat="identity") + coord_cartesian(ylim=c(0.0007,0.00085)) +
  scale_fill_manual(values=c("#E69F00", "#D55E00", "#DE1C14")) +
  ggtitle("pval = 7.489555e-03\nenr = 1.043168") + ylab("Number of functional SNPs per bp\n(TSS +- 2.5kb)") +
  theme_bw() + theme(aspect.ratio=2, plot.title = element_text(size = 11)) + theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 0.5)) + theme(axis.title.x = element_blank())

# Save plots

pdf("./Rplot/Fig_S5B_S5C.pdf", width=6, height=4, useDingbats=FALSE)
ggarrange(count.plot, count.eff.plot, ncol=2, nrow=1)
dev.off()

############################################################################
############################################################################
###                 Enrichment of functional SNPs at                     ###
###                 splice sites marked with origins                     ###
############################################################################
############################################################################

############################################################################
# Load data
population.SNP.closest.exon.intron.junction <- read.table("./Figure_4/population.SNP.closest.exon.intron.junctions.bed", header = F) %>% dplyr::select(V1, V2, V4, V5, V6, V8, V10, V11, V12)
colnames(population.SNP.closest.exon.intron.junction) <- c("chr", "SNP.pos", "SNP.id", "REF", "ALT", "splice.site.pos", "gene.strand", "gene.id", "gene.name")
population.SNP.closest.intron.exon.junction <- read.table("./Figure_4/population.SNP.closest.intron.exon.junctions.bed", header = F) %>% dplyr::select(V1, V2, V4, V5, V6, V8, V10, V11, V12)
colnames(population.SNP.closest.intron.exon.junction) <- c("chr", "SNP.pos", "SNP.id", "REF", "ALT", "splice.site.pos", "gene.strand", "gene.id", "gene.name")
# Combine exon/intron and intron/exon SNPs
population.SNP.closest.splice.site <- rbind(population.SNP.closest.exon.intron.junction, population.SNP.closest.intron.exon.junction)
# Select only SNP in 10 kb around exon/intron junction
population.SNP.closest.splice.site.10kb <- population.SNP.closest.splice.site %>% mutate(SNP.gene.dist = SNP.pos-splice.site.pos) %>% filter(SNP.gene.dist >= -10000 & SNP.gene.dist <= 10000)
nrow(population.SNP.closest.splice.site.10kb) # 22,600,903 SNPs
# Save
write.table(population.SNP.closest.splice.site.10kb, "./Figure_4/population.SNP.closest.splice.site.10kb.bed", quote=FALSE, sep='\t', col.names = T)
population.SNP.closest.splice.site.10kb  <- read.csv("./Figure_4/population.SNP.closest.splice.site.10kb.bed", sep = '\t')

############################################################################
# Select functional SNPs at splice sites
functional.SNP.closest.splice.site.10kb <- population.SNP.closest.splice.site.10kb %>% filter(SNP.id %in% SNP.functional$rsid) # 560,099 functional SNPs

# Load results from Fig. 4
hsapiens.exon.intron.junctions.closest.IniSeq.ori <- read.csv("./Figure_4/hsapiens.exon.intron.junctions.closest.IniSeq.center.bed", header = F, sep = "\t") %>% dplyr::select(V1, V2, V4, V5, V6, V8, V10, V11)
colnames(hsapiens.exon.intron.junctions.closest.IniSeq.ori) <- c("chr", "exon.intron.pos", "gene.strand", "gene.id", "gene.name", "ori.pos", "ori.eff", "dist")
hsapiens.intron.exon.junctions.closest.IniSeq.ori <- read.csv("./Figure_4/hsapiens.intron.exon.junctions.closest.IniSeq.center.bed", header = F, sep = "\t") %>% dplyr::select(V1, V2, V4, V5, V6, V8, V10, V11)
colnames(hsapiens.intron.exon.junctions.closest.IniSeq.ori) <- c("chr", "intron.exon.pos", "gene.strand", "gene.id", "gene.name", "ori.pos", "ori.eff", "dist")
hsapiens.exon.intron.junctions.closest.IniSeq.ori.10kb <- hsapiens.exon.intron.junctions.closest.IniSeq.ori %>% filter(dist >= -10000 & dist <= 10000) %>%
  mutate(dist = ((as.integer(cut(dist, breaks = 4000)))-2000)*5) # per 5bp bin
hsapiens.intron.exon.junctions.closest.IniSeq.ori.10kb <- hsapiens.intron.exon.junctions.closest.IniSeq.ori %>% filter(dist >= -10000 & dist <= 10000) %>%
  mutate(dist = ((as.integer(cut(dist, breaks = 4000)))-2000)*5) # per 5bp bin

# Genes marked by origins at exon/intron junctions

genes.exon.intron.junction.with.ori <- hsapiens.exon.intron.junctions.closest.IniSeq.ori.10kb %>% filter(dist <= 1250) %>%
  mutate(junction.id = paste(chr, exon.intron.pos, sep ="_")) # 50,276 junctions
colnames(genes.exon.intron.junction.with.ori) <- c("chr", "splice.site.pos", "gene.strand", "gene.id", "gene.name", "ori.pos", "ori.eff", "dist", "junction.id")
genes.exon.intron.junction.without.ori <- hsapiens.exon.intron.junctions.closest.IniSeq.ori %>%
  mutate(junction.id = paste(chr, exon.intron.pos, sep ="_")) %>% filter(!junction.id %in% genes.exon.intron.junction.with.ori$junction.id) # 163,867 junctions
colnames(genes.exon.intron.junction.without.ori) <- c("chr", "splice.site.pos", "gene.strand", "gene.id", "gene.name", "ori.pos", "ori.eff", "dist", "junction.id")
# 23.48% of exon/intron junctions have an origin

genes.intron.exon.junction.with.ori <- hsapiens.intron.exon.junctions.closest.IniSeq.ori.10kb %>% filter(dist <= 1250) %>%
  mutate(junction.id = paste(chr, intron.exon.pos, sep ="_")) # 50,094 junctions
colnames(genes.intron.exon.junction.with.ori) <- c("chr", "splice.site.pos", "gene.strand", "gene.id", "gene.name", "ori.pos", "ori.eff", "dist", "junction.id")
genes.intron.exon.junction.without.ori <- hsapiens.intron.exon.junctions.closest.IniSeq.ori %>%
  mutate(junction.id = paste(chr, intron.exon.pos, sep ="_")) %>% filter(!junction.id %in% genes.intron.exon.junction.with.ori$junction.id) # 164,049 junctions
colnames(genes.intron.exon.junction.without.ori) <- c("chr", "splice.site.pos", "gene.strand", "gene.id", "gene.name", "ori.pos", "ori.eff", "dist", "junction.id")
# 23.39% of intron/exon junctions have an origin

# Define splice sites marked by origins

splice.sites.with.ori <- rbind(genes.exon.intron.junction.with.ori, genes.intron.exon.junction.with.ori)
splice.sites.without.ori <- rbind(genes.exon.intron.junction.without.ori, genes.intron.exon.junction.without.ori)

# Save lists
write.csv(splice.sites.with.ori, "./Figure_5/splice.sites.with.ori.bed", quote=FALSE)
write.csv(splice.sites.without.ori, "./Figure_5/splice.sites.without.ori.bed", quote=FALSE)

# Select functional SNPs associated with splice sites marked by an origin
functional.SNP.closest.splice.sites.10kb.with.ori <- functional.SNP.closest.splice.site.10kb %>% mutate(junction.id = paste(chr, splice.site.pos, sep = "_")) %>% 
  filter(junction.id %in% unique(splice.sites.with.ori$junction.id)) # 55,903 junctions
functional.SNP.closest.splice.sites.10kb.without.ori <- functional.SNP.closest.splice.site.10kb %>% mutate(junction.id = paste(chr, splice.site.pos, sep = "_")) %>% 
  filter(junction.id %in% unique(splice.sites.without.ori$junction.id)) # 121,511 junctions

############################################################################
# Compute enrichment profile from background in windows of 5nt and averaged over 20 windows

functional.SNP.closest.splice.sites.10kb.with.ori.dens <- functional.SNP.closest.splice.sites.10kb.with.ori %>% mutate(dist = ((as.integer(cut(SNP.gene.dist, breaks = 4000)))-2000)*5) %>% 
  group_by(dist) %>% summarise(count = n()) %>% mutate(enr = log2(count/mean(count[c(1:100,900:1000)])), enr.average=rollapply(enr,20,mean,align='right',fill=0), class = "with origin")

functional.SNP.closest.splice.sites.10kb.without.ori.dens <- functional.SNP.closest.splice.sites.10kb.without.ori %>% mutate(dist = ((as.integer(cut(SNP.gene.dist, breaks = 4000)))-2000)*5) %>% 
  group_by(dist) %>% summarise(count = n()) %>% mutate(enr = log2(count/mean(count[c(1:100,900:1000)])), enr.average=rollapply(enr,20,mean,align='right',fill=0), class = "without origin")

Splice.sites.plot <- rbind(functional.SNP.closest.splice.sites.10kb.with.ori.dens, functional.SNP.closest.splice.sites.10kb.without.ori.dens)

pdf("./Rplot/Fig_5B_2.pdf", width=6, height=4, useDingbats=FALSE)
Splice.sites.plot %>% ggplot(aes(x=dist, y=enr.average, color=class)) +
  geom_line(size = 0.5) + xlim(-5000,5000) + ylim(0,4.6) +
  scale_color_manual(values=c("#009E73", "#BDBDBD")) +
  xlab("Distance from splice junction (bp)") + ylab("Functional SNP enrichment (log2 FC)") + ggtitle("Functional SNPs around gene starts") +
  theme_bw() + theme(aspect.ratio=1.1, plot.title = element_text(size=12))
dev.off()

############################################################################
# Quantification and statistics

# Select functional SNPs associated with splice sites marked by an origin binned by efficiency
splice.sites.with.ori.low <- splice.sites.with.ori %>% filter(ori.eff < 0.8418)  # 30,746 sites
splice.sites.with.ori.medium <- splice.sites.with.ori %>% filter(ori.eff >= 0.8418 & ori.eff < 0.9065)  # 33,010 sites
splice.sites.with.ori.high <- splice.sites.with.ori %>% filter(ori.eff >= 0.9065)  # 31,065 sites
functional.SNP.closest.splice.site.10kb.with.ori <- functional.SNP.closest.splice.site.10kb %>% mutate(junction.id = paste(chr, splice.site.pos, sep = "_")) %>% filter(junction.id %in% unique(splice.sites.with.ori$junction.id)) # 55,903 sites
functional.SNP.closest.splice.site.10kb.without.ori <- functional.SNP.closest.splice.site.10kb %>% mutate(junction.id = paste(chr, splice.site.pos, sep = "_")) %>% filter(junction.id %in% unique(splice.sites.without.ori$junction.id)) # 121,511 sites
functional.SNP.closest.splice.site.10kb.with.ori.low <- functional.SNP.closest.splice.site.10kb %>% mutate(junction.id = paste(chr, splice.site.pos, sep = "_")) %>% filter(junction.id %in% unique(splice.sites.with.ori.low$junction.id)) # 16,215 sites
functional.SNP.closest.splice.site.10kb.with.ori.medium <- functional.SNP.closest.splice.site.10kb %>% mutate(junction.id = paste(chr, splice.site.pos, sep = "_")) %>% filter(junction.id %in% unique(splice.sites.with.ori.medium$junction.id)) # 19,233 sites
functional.SNP.closest.splice.site.10kb.with.ori.high <- functional.SNP.closest.splice.site.10kb %>% mutate(junction.id = paste(chr, splice.site.pos, sep = "_")) %>% filter(junction.id %in% unique(splice.sites.with.ori.high$junction.id)) # 20,455 sites

# Compute the number of functional SNPs in the vicinity of splice sites (+- 2.5kb) marked by an origins

functional.SNP.closest.splice.site.with.ori <- functional.SNP.closest.splice.site.10kb.with.ori %>% filter(SNP.gene.dist >= -2500 & SNP.gene.dist <= 2500)  # 144,577 mutations for 53,105 sites
functional.SNP.closest.splice.site.without.ori <- functional.SNP.closest.splice.site.10kb.without.ori %>% filter(SNP.gene.dist >= -2500 & SNP.gene.dist <= 2500)  # 198,426 mutations for 97,466 sites
functional.SNP.closest.splice.site.with.ori.low <- functional.SNP.closest.splice.site.10kb.with.ori.low %>% filter(SNP.gene.dist >= -2500 & SNP.gene.dist <= 2500)  # 39,593 mutations for 15,574 sites
functional.SNP.closest.splice.site.with.ori.medium <- functional.SNP.closest.splice.site.10kb.with.ori.medium %>% filter(SNP.gene.dist >= -2500 & SNP.gene.dist <= 2500)  # 48,769 mutations for 18,302 sites
functional.SNP.closest.splice.site.with.ori.high <- functional.SNP.closest.splice.site.10kb.with.ori.high %>% filter(SNP.gene.dist >= -2500 & SNP.gene.dist <= 2500)  # 56,215 mutations for 19,229 sites

# Prepare count tables

count.df1 <- tibble(Count = 144577/(53105*5000), class = "with ori")
count.df2 <- tibble(Count = 198426/(97466*5000), class = "without ori")
count.df <- rbind(count.df1, count.df2)
# 2x2 chi-squared significance test (comparing genes with or without ori)
pval.table <- rbind(c(144577,(53105*5000)-144577),c(198426,(97466*5000)-198426))
pval <- chisq.test(pval.table) # p-value < 2.2e-16, enr = 1.337268

count.splice.plot <- count.df %>% mutate(Class = fct_relevel(class, "without ori", "with ori")) %>%
  ggplot(aes(x=Class, y=Count, fill=Class)) +
  geom_bar(stat="identity") + ylim(0,0.0008) +
  scale_fill_manual(values=c("#D4D2D2", "#009E73")) +
  ggtitle("p-value < 2.2e-16\nenr = 1.337268") + ylab("Number of functional SNPs per bp\n(Splice sites +- 2.5kb)") +
  theme_bw() + theme(aspect.ratio=2, plot.title = element_text(size = 11)) + theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 0.5)) + theme(axis.title.x = element_blank())

count.eff.df1 <- tibble(Count = 39593/(15574*5000), class = "low")
count.eff.df2 <- tibble(Count = 48769/(18302*5000), class = "medium")
count.eff.df3 <- tibble(Count = 56215/(19229*5000), class = "high")
count.eff.df <- rbind(count.eff.df1, count.eff.df2, count.eff.df3)
# 2x2 chi-squared significance test (comparing genes low with high)
pval.table <- rbind(c(39593,(15574*5000)-39593),c(56215,(19229*5000)-56215))
pval <- chisq.test(pval.table) # p-value = 8.263174e-101, enr = 1.149946

count.splice.eff.plot <- count.eff.df %>% mutate(Class = fct_relevel(class, "low", "medium", "high")) %>%
  ggplot(aes(x=Class, y=Count, fill=Class)) +
  geom_bar(stat="identity") + coord_cartesian(ylim=c(0.00035,0.00065)) +
  scale_fill_manual(values=c("#C3D7A4", "#ADCC54", "#009E73")) +
  ggtitle("pval = 8.263174e-101\nenr = 1.149946") + ylab("Number of functional SNPs per bp\n(Splice sites +- 2.5kb)") +
  theme_bw() + theme(aspect.ratio=2, plot.title = element_text(size = 11)) + theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 0.5)) + theme(axis.title.x = element_blank())

# Save plots

pdf("./Rplot/Fig_S5B_S5C_2.pdf", width=6, height=4, useDingbats=FALSE)
ggarrange(count.splice.plot, count.splice.eff.plot, ncol=2, nrow=1)
dev.off()

############################################################################
############################################################################
###                  Functional SNPs enrichment at                       ###
###                      TSS and splice sites                            ###
###    - influence of origin efficiency and germline gene expression -  ###
############################################################################
############################################################################

############################################################################
# Germline gene expression data is from
# Sohni et al. Cell reports 2019, 26, 1501–1517
# We considered single cell rna-seq data from two unfractioned adult testes
# sample (A1 and A2) from the Spermatogonia (SPG), Spermatocytes (SPC) and
# the Spermatid (St) lineages
# from the GSE124263_RAW.tar file available at
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE124263

# Load A1 data
counts.A1 <- readMM("./Misc/GSE124263/GSE124263_RAW/GSM3526588_A1Total_matrix.mtx.gz")
genes.A1 <- read_tsv("./Misc/GSE124263/GSE124263_RAW/GSM3526588_A1Total_genes.tsv.gz", col_names = FALSE)
gene.ids.A1 <- genes.A1$X1
cell.ids.A1 <- read_tsv("./Misc/GSE124263/GSE124263_RAW/GSM3526588_A1Total_barcodes.tsv.gz", col_names = FALSE)$X1
# Make the column names as the cell IDs and the row names as the gene IDs
rownames(counts.A1) <- gene.ids.A1
colnames(counts.A1) <- cell.ids.A1
# Transform into a dataframe
counts.A1.df <- as.data.frame(as.matrix(counts.A1))

# Load A2 data
counts.A2 <- readMM("./Misc/GSE124263/GSE124263_RAW/GSM3526590_A2_total_matrix.mtx.gz")
genes.A2 <- read_tsv("./Misc/GSE124263/GSE124263_RAW/GSM3526590_A2_total_genes.tsv.gz", col_names = FALSE)
gene.ids.A2 <- genes.A2$X1
cell.ids.A2 <- read_tsv("./Misc/GSE124263/GSE124263_RAW/GSM3526590_A2_total_barcodes.tsv.gz", col_names = FALSE)$X1
# Make the column names as the cell IDs and the row names as the gene IDs
rownames(counts.A2) <- gene.ids.A2
colnames(counts.A2) <- cell.ids.A2
# Transform into a dataframe
counts.A2.df <- as.data.frame(as.matrix(counts.A2))

# Load cell identity information

# From the supplementary table, save a tsv for Adult cells
# cell.identity.csv
cell.identity.df <- read.csv("./Misc/GSE124263/cell.identity.csv")
unique(cell.identity.df$AdultCellAnnotation)
# Select germ cells (SPG, SPC and St)
cell.identity.GM.df <- cell.identity.df %>% filter(AdultCellAnnotation == "Spermatogonia (SPG) (37%)" | AdultCellAnnotation == "Spermatocytes (SPC) (6%)" | AdultCellAnnotation == "Spermatid (St) (2%)")
# Rename cell type
cell.identity.GM.df <- cell.identity.GM.df %>%
  mutate(CellType = case_when(AdultCellAnnotation == "Spermatogonia (SPG) (37%)" ~ "SPG",
                              AdultCellAnnotation == "Spermatocytes (SPC) (6%)" ~ "SPC",
                              TRUE ~ "St")) %>% dplyr::select(-AdultCellAnnotation)
# Select cell identity for A1 and A2 (considering only unfractioned samples)
cell.identity.GM.A1.df <- cell.identity.GM.df %>% filter(orig.ident == "A1_T") %>% dplyr::select(-orig.ident) %>% separate(Cell.Barcode, c("Sample", "UMI")) %>% dplyr::select(-Sample)
cell.identity.GM.A2.df <- cell.identity.GM.df %>% filter(orig.ident == "A2_T") %>% dplyr::select(-orig.ident) %>% separate(Cell.Barcode, c("Sample", "UMI")) %>% dplyr::select(-Sample)

# Parse count matrices

# Transform column name to recover UMI

UMI.A1 <- vector()
for (i in colnames(counts.A1.df)) {
  p <- str_split(colnames(counts.A1.df[i]), "-")[[1]][1]
  UMI.A1 <- append(UMI.A1, p)
}
colnames(counts.A1.df) <- UMI.A1

UMI.A2 <- vector()
for (i in colnames(counts.A2.df)) {
  p <- str_split(colnames(counts.A2.df[i]), "-")[[1]][1]
  UMI.A2 <- append(UMI.A2, p)
}
colnames(counts.A2.df) <- UMI.A2

# Prepare individual matrices

A1.SPG <- counts.A1.df %>% dplyr::select(cell.identity.GM.A1.df[which(cell.identity.GM.A1.df$CellType == "SPG"),]$UMI) # 911 cells
A1.SPC <- counts.A1.df %>% dplyr::select(cell.identity.GM.A1.df[which(cell.identity.GM.A1.df$CellType == "SPC"),]$UMI) # 324 cells
A1.St <- counts.A1.df %>% dplyr::select(cell.identity.GM.A1.df[which(cell.identity.GM.A1.df$CellType == "St"),]$UMI) # 105 cells
A2.SPG <- counts.A2.df %>% dplyr::select(cell.identity.GM.A2.df[which(cell.identity.GM.A2.df$CellType == "SPG"),]$UMI) # 2867 cells
A2.SPC <- counts.A2.df %>% dplyr::select(cell.identity.GM.A2.df[which(cell.identity.GM.A2.df$CellType == "SPC"),]$UMI) # 397 cells
A2.St <- counts.A2.df %>% dplyr::select(cell.identity.GM.A2.df[which(cell.identity.GM.A2.df$CellType == "St"),]$UMI) # 139 cells
# Compute mean counts
A1.SPG.mean <- cbind.data.frame(gene.id = row.names(A1.SPG), SPG.A1 = rowMeans(A1.SPG))
A1.SPC.mean <- cbind.data.frame(gene.id = row.names(A1.SPC), SPC.A1 = rowMeans(A1.SPC))
A1.St.mean <- cbind.data.frame(gene.id = row.names(A1.St), St.A1 = rowMeans(A1.St))
A2.SPG.mean <- cbind.data.frame(gene.id = row.names(A2.SPG), SPG.A2 = rowMeans(A2.SPG))
A2.SPC.mean <- cbind.data.frame(gene.id = row.names(A2.SPC), SPC.A2 = rowMeans(A2.SPC))
A2.St.mean <- cbind.data.frame(gene.id = row.names(A2.St), St.A2 = rowMeans(A2.St))
# Combine all data
germcell.expression.df <- A1.SPG.mean %>% full_join(A1.SPC.mean, by = "gene.id") %>% 
  full_join(A1.St.mean, by = "gene.id") %>% 
  full_join(A2.SPG.mean, by = "gene.id") %>% 
  full_join(A2.SPC.mean, by = "gene.id") %>% 
  full_join(A2.St.mean, by = "gene.id")
# Combine replicates
germcell.exp.df <- germcell.expression.df %>% mutate(SPG = (SPG.A1+SPG.A2)/2, SPC = (SPC.A1+SPC.A2)/2, St = (St.A1+St.A2)/2) %>% dplyr::select(gene.id, SPG, SPC, St)
# Compute averages between the three lineages
germcell.exp.df <- germcell.exp.df %>% mutate(EXP = (SPG+SPC+St)/3)
# Compute quantiles (excluding 0 values)
germcell.exp.df.2 <- germcell.exp.df %>% filter(EXP != 0)
quantile(germcell.exp.df.2$EXP, c(1/3,2/3))
# 33.33333% 66.66667% 
# 0.005922367 0.089219808 
# Add gene expression information
germcell.exp.df <- germcell.exp.df %>% mutate(exp = case_when(EXP == 0 ~ "non",
                                                              EXP > 0 & EXP <= 0.005922367 ~ "low",
                                                              EXP > 0.005922367 & EXP <= 0.089219808 ~ "medium",
                                                              T ~ "high"))
table(germcell.exp.df$exp)
#   high    low   medium  non 
#   9098    9098   9097   3604
germcell.exp.df.2 <- germcell.exp.df %>% dplyr::select(gene.id, EXP, exp)

############################################################################
# Prepare functional SNP coverage bigwig file
# Original coordinates are in hg19

SNP.functional.hg19.bed <- SNP.functional %>% dplyr::select(chrom, start, end)
write.table(SNP.functional.hg19.bed, "./Dataset/SNP.functional.hg19.bed", sep="\t", col.names = F, row.names = F, quote = F)

# using zsh interactive shell on Mac with bedtools2 and UCSC liftOver
# liftOver SNP.functional.hg19.bed hg19ToHg38.over.chain SNP.functional.hg38.bed unMapped.bed
# sort -k1,1 -k2,2n SNP.functional.hg38.bed > SNP.functional.hg38.sorted.bed
# windowBed -c -w 0 -a hg38.tile.50nt.bed -b SNP.functional.hg38.sorted.bed > SNP.functional.hg38.count.50nt.bedgraph
# LC_COLLATE=C sort -k1,1 -k2,2n SNP.functional.hg38.count.50nt.bedgraph > SNP.functional.hg38.count.50nt.sorted.bedgraph
# bedGraphToBigWig SNP.functional.hg38.count.50nt.sorted.bedgraph hg38.chrom.sizes SNP.functional.hg38.count.50nt.sorted.bw

############################################################################
# Count the number of functional SNP at each TSS/junctions

# Load grange object for functional SNP
functional.SNP.gr <- import("./Dataset/SNP.functional.hg38.count.50nt.sorted.bw")

# Load results from Fig. 4
genes.TSS.with.ori <- read.csv("./Figure_4/genes.TSS.with.ori.bed", row.names = NULL) %>% dplyr::select(-X) # 11,410
genes.TSS.without.ori <- read.csv("./Figure_4/genes.TSS.without.ori.bed", row.names = NULL) %>% dplyr::select(-X) # 7,925
splice.sites.with.ori <- read.csv("./Figure_4/splice.sites.with.ori.bed", row.names = NULL) %>% dplyr::select(-X) # 100,370
splice.sites.without.ori <- read.csv("./Figure_4/splice.sites.without.ori.bed", row.names = NULL) %>% dplyr::select(-X) # 327,916
# Transform
genes.TSS.with.ori.bed <- genes.TSS.with.ori %>% mutate(ori.id = paste(chr, ori.pos, sep = "_"), start = as.integer(gene.start), end = as.integer(gene.start+1)) %>% 
  dplyr::select(chr, start, end, gene.id, gene.name, strand = gene.strand, ori.id, EFF = ori.eff, dist) %>% mutate(strand = case_when(strand == 1 ~ "+", T ~ "-"))  %>% mutate(eff = case_when(EFF < 0.8418 ~ "low",
                                                                                                                                                                                               EFF >= 0.8418 & EFF < 0.9065 ~ "medium",
                                                                                                                                                                                               T ~ "high"))
genes.TSS.without.ori.bed <- genes.TSS.without.ori %>% mutate(start = as.integer(gene.start), end = as.integer(gene.start+1)) %>% 
  dplyr::select(chr, start, end, gene.id, gene.name, strand = gene.strand) %>% mutate(strand = case_when(strand == 1 ~ "+", T ~ "-"))

splice.sites.with.ori.bed <- splice.sites.with.ori %>% mutate(ori.id = paste(chr, ori.pos, sep = "_"), start = as.integer(splice.site.pos), end = as.integer(splice.site.pos+1)) %>% 
  dplyr::select(chr, start, end, gene.id, gene.name, strand = gene.strand, ori.id, EFF = ori.eff, dist) %>% mutate(strand = case_when(strand == 1 ~ "+", T ~ "-")) %>% mutate(strand = case_when(strand == 1 ~ "+", T ~ "-"))  %>% mutate(eff = case_when(EFF < 0.8418 ~ "low",
                                                                                                                                                                                                                                                          EFF >= 0.8418 & EFF < 0.9065 ~ "medium",
                                                                                                                                                                                                                                                          T ~ "high"))
splice.sites.without.ori.bed <- splice.sites.without.ori %>% mutate(start = as.integer(splice.site.pos), end = as.integer(splice.site.pos+1)) %>% 
  dplyr::select(chr, start, end, gene.id, gene.name, strand = gene.strand) %>% mutate(strand = case_when(strand == 1 ~ "+", T ~ "-"))

# Prepare TSS/junction gr files
# setting origin distances at 2.5 kb
d <- 2500

genes.TSS.with.ori.bed.2.5.kb <- genes.TSS.with.ori.bed %>% mutate(start = as.integer(start-d), end = as.integer(end+d))
genes.TSS.without.ori.bed.2.5.kb <- genes.TSS.without.ori.bed %>% mutate(start = as.integer(start-d), end = as.integer(end+d)) %>% 
  mutate(ori.id = NA, EFF = NA, dist = NA, eff = "non")
genes.TSS.bed.2.5.kb <- rbind(genes.TSS.with.ori.bed.2.5.kb, genes.TSS.without.ori.bed.2.5.kb) %>% mutate(TSS.id = paste(chr, start, end, sep = "_")) %>% left_join(germcell.exp.df.2, by = "gene.id")
TSS.2.5.kb.gr <- makeGRangesFromDataFrame(genes.TSS.bed.2.5.kb, keep.extra.columns=T)
seqlevelsStyle(TSS.2.5.kb.gr) <- "UCSC"

splice.sites.with.ori.bed.2.5.kb <- splice.sites.with.ori.bed %>% mutate(start = as.integer(start-d), end = as.integer(end+d))
splice.sites.without.ori.bed.2.5.kb <- splice.sites.without.ori.bed %>% mutate(start = as.integer(start-d), end = as.integer(end+d)) %>% 
  mutate(ori.id = NA, EFF = NA, dist = NA, eff = "non")
junctions.bed.2.5.kb <- rbind(splice.sites.with.ori.bed.2.5.kb, splice.sites.without.ori.bed.2.5.kb) %>% mutate(junction.id = paste(chr, start, end, sep = "_")) %>% left_join(germcell.exp.df.2, by = "gene.id")
junctions.2.5.kb.gr <- makeGRangesFromDataFrame(junctions.bed.2.5.kb, keep.extra.columns=T)
seqlevelsStyle(junctions.2.5.kb.gr) <- "UCSC"

# Compute number of functional SNP at each site

TSS.func.SNP <- as.data.frame(mergeByOverlaps(TSS.2.5.kb.gr, functional.SNP.gr)) %>% group_by(TSS.id) %>% summarise(SNP.count = sum(score))
genes.TSS.SNP.count <- genes.TSS.bed.2.5.kb %>% left_join(TSS.func.SNP, by = "TSS.id")
# Drop NA, add a pseudo count of one and compute mutational burden (expressed in x 10e-3 per bp)
genes.TSS.SNP.count.2 <- genes.TSS.SNP.count %>% drop_na(exp) %>% mutate(SNP.count = (SNP.count + 1)/5)

junctions.func.SNP <- as.data.frame(mergeByOverlaps(junctions.2.5.kb.gr, functional.SNP.gr)) %>% group_by(junction.id) %>% summarise(SNP.count = sum(score))
junctions.SNP.count <- junctions.bed.2.5.kb %>% left_join(junctions.func.SNP, by = "junction.id")
# Drop NA and add a pseudo count of one
junctions.SNP.count.2 <- junctions.SNP.count %>% drop_na(exp) %>% mutate(SNP.count = (SNP.count + 1)/5)

SNP.TSS.exp.plot <- genes.TSS.SNP.count.2 %>% mutate(exp2 = fct_relevel(exp, "non", "low", "medium", "high"), eff2 = fct_relevel(eff, "non","low", "medium", "high")) %>%
  ggplot(aes(exp2, SNP.count, fill = eff2)) +
  geom_boxplot() + ggtitle("TSS\nP = 2.2e-16, 2.2e-16, 2.2e-16, 2.2e-16") + ylab("Functional SNPs load") +
  scale_fill_manual(values = c("#D4D2D2", "#CC79A7", "#0072B2","#FC4E07")) +
  scale_y_continuous(trans = 'log10') + labs(fill = "Origin efficiency") +
  theme_bw() + theme(aspect.ratio=1, plot.title = element_text(size=14), axis.title=element_text(size=14), axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1), axis.text.y = element_text(size = 12), axis.title.x=element_blank(), axis.title.y=element_text(vjust=4))

SNP.junctions.exp.plot <- junctions.SNP.count.2 %>% mutate(exp2 = fct_relevel(exp, "non", "low", "medium", "high"), eff2 = fct_relevel(eff, "non", "low", "medium", "high")) %>%
  ggplot(aes(exp2, SNP.count, fill = eff2)) +
  geom_boxplot() + ggtitle("Junctions\nP = 2.2e-16, 2.2e-16, 2.2e-16, 2.2e-16") + ylab("Functional SNPs load") +
  scale_fill_manual(values = c("#D4D2D2", "#CC79A7", "#0072B2","#FC4E07")) +
  scale_y_continuous(trans = 'log10') + labs(fill = "Origin efficiency") +
  theme_bw() + theme(aspect.ratio=1, plot.title = element_text(size=14), axis.title=element_text(size=14), axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1), axis.text.y = element_text(size = 12), axis.title.x=element_blank(), axis.title.y=element_text(vjust=4))

pdf("./Rplot/Fig_5C.pdf", width=12, height=4, useDingbats=FALSE)
ggarrange(SNP.TSS.exp.plot, SNP.junctions.exp.plot, ncol = 2)
dev.off()

############################################################################
############################################################################
###                 Enrichment of functional SNPs at                     ###
###          TFBS within the promoter of genes marked by origins         ###
############################################################################
############################################################################

############################################################################
# TF binding sites were recovered from the JASPAR database at
# http://jaspar.genereg.net/
# All human non-redundant TF core binding sites PWM were downloaded

# The name and ID of 639 human core TF were manually curated from the JASPAR website
# to generate JASPAR_human_core_TF.csv

# Recover human TF binding sites (PWM) ID from the Jaspar database
human.core.TF <- read.csv("./Dataset/JASPAR_human_core_TF.csv")
# Split TF ID and remove duplicates
human.core.TF <- human.core.TF %>% tidyr::separate(ID, c("ID", "version", sep = ".")) %>% dplyr::select(-., -Logo, -version) %>% distinct(ID, .keep_all = TRUE)
# contains 639 individual transcription factors with associated PWM

# Prepare PWM matrix list
pwmList.human.TF <- PWMatrixList(use.names=TRUE)
TF.PFM.jaspar <- list.files("./Dataset/JASPAR2020_CORE_non-redundant_pfms_jaspar/")
for(name in human.core.TF$ID){
  num <- grep(name,TF.PFM.jaspar)
  path <- paste("./Dataset/JASPAR2020_CORE_non-redundant_pfms_jaspar/", TF.PFM.jaspar[num], sep = "")
  p <- readJASPARMatrix(path, matrixClass= "PFM")[[1]]
  q <- toPWM(p)
  r <- PWMatrixList(q, use.names=TRUE)
  pwmList.human.TF <- c(pwmList.human.TF, r)
}

# Analyze length of binding motifs in the PWM list to determine the length of the sequences to screen
TFBS.length <- vector()
for (i in 1:length(pwmList.human.TF)) {
  p <- length(pwmList.human.TF[[i]])
  TFBS.length <- append(TFBS.length, p)
}
# summary(TFBS.length)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 6.00   10.00   12.00   12.05   14.00   24.00 
# Then screen 24nt sequences (+- 12 nt sequences around functional SNPs) and select only TFs whose binding sites overlap with the functional SNPs (start <= 12 & end >= 12)

############################################################################
# Recover genomic location of functional SNPs (+- 2.5kb around gene starts)

functional.SNP.closest.gene.2.5.kb <- functional.SNP.closest.gene.10kb %>% filter(dist >= -2500 & dist <= 2500) %>% mutate(start = SNP.pos-12, end = SNP.pos+12) %>% 
  filter(start > 0 & end > 0) %>% dplyr::select(chr, start, end, SNP.id) %>% mutate(chr = paste("chr",chr,sep=""))
# 57,203 functional SNPs
# Save
write.table(functional.SNP.closest.gene.2.5.kb, "./Figure_4/functional.SNP.for.TF.analysis.bed", sep="\t", col.names = F, row.names = F, quote = F)
# Create grange object
functional.SNP.gr <- import("./Figure_4/functional.SNP.for.TF.analysis.bed")
# Recover sequences
functional.SNP.seq <- getSeq(Hsapiens, functional.SNP.gr)
# Recover names
functional.SNP.names <- functional.SNP.gr$name

############################################################################
# Identify TFBSs within sequences
# Multiple transcription factors can be associated with a given SNPs
# TF hits are selected based on a minimum score of 95%

TF.df <- tibble(SNP.id = character(), TF = character())
for (i in 1:length(functional.SNP.seq)) {
  print(i/length(functional.SNP.seq))
  print(Sys.time())
  SNP.seq <- functional.SNP.seq[i]
  SNP.name <- functional.SNP.gr$name[i]
  TF <- searchSeq(pwmList.human.TF, SNP.seq, seqname= SNP.name, min.score="95%", strand="*")
  TF.GFF3 <- writeGFF3(TF)
  if (nrow(TF.GFF3) == 0) {
    p <- tibble(SNP.id = SNP.name, TF = NA)
  } else {
    TF.result <- TF.GFF3 %>% filter(score >= 10 & start <= 12 & end >= 12) %>% tidyr::separate(attributes, c("TF", "class", "sequence"), sep = ";") %>% dplyr::select(TF)
    if (nrow(TF.result) == 0) {
      p <- tibble(SNP.id = SNP.name, TF = NA)
    } else {
      TF.result.2 <- unique(str_replace(TF.result$TF, "TF=", ""))
      p <- tibble(SNP.id = SNP.name, TF = TF.result.2)
    }
  }
  TF.df <- rbind(TF.df, p)
}

# Save results
write.csv(TF.df, "./Figure_4/functional.SNP.TFBS.csv")
SNP.functional.TFBS <- read.csv("./Figure_4/functional.SNP.TFBS.csv") %>% dplyr::select(-X)

# Prepare table for TFBS enrichment
SNP.functional.TFBS.gene.with.ori <- SNP.functional.TFBS %>% filter(SNP.id %in% functional.SNP.closest.gene.start.with.ori$SNP.id) # 10,327 genes
SNP.functional.TFBS.gene.without.ori <- SNP.functional.TFBS %>% filter(SNP.id %in% functional.SNP.closest.gene.start.without.ori$SNP.id) # 5,086 genes

SNP.functional.TFBS.gene.with.ori.table <- as.data.frame(table(SNP.functional.TFBS.gene.with.ori$TF)) %>% dplyr::select(TF = Var1, Count.with = Freq)
SNP.functional.TFBS.gene.without.ori.table <- as.data.frame(table(SNP.functional.TFBS.gene.without.ori$TF)) %>% dplyr::select(TF = Var1, Count.without = Freq)
SNP.functional.TFBS.gene.summary <- SNP.functional.TFBS.gene.with.ori.table %>% full_join(SNP.functional.TFBS.gene.without.ori.table, by = "TF")
SNP.functional.TFBS.gene.summary[is.na(SNP.functional.TFBS.gene.summary)] <- 0

############################################################################
# Compute enrichment and associated pvalues pvalues (2x2 chi-squared test)
pval <- vector()
enr <- vector()
for (i in 1:nrow(SNP.functional.TFBS.gene.summary)) {
  TF <- SNP.functional.TFBS.gene.summary[i,]
  chi.table <- rbind(c(TF$Count.with,10154), c(TF$Count.without,5259))
  p <- chisq.test(chi.table)$p.value
  pval <- append(pval,p)
  q <- (TF$Count.with/10154)/(TF$Count.without/5259)
  enr <- append(enr,q)
}
SNP.functional.TFBS.gene.summary.stat <- cbind.data.frame(SNP.functional.TFBS.gene.summary, pval,enr)

# Compute expected pvalues for quantile-quantile plots of TF associations in gene starts

# Define function that input a dataframe with a list of pvalues associated to TFs and return the observed and expected pvalues in -log10
QQfunct = function (stat.df) {  # Input is a dataframe
  observedPValues <- stat.df$pval
  pval.obs.log10 <- -log10(sort(observedPValues))
  pval.exp.log10 <- -log10(1:length(observedPValues)/length(observedPValues))
  # Reorder original table
  stat.sort.df <- stat.df[order(stat.df$pval),]
  # Bind results
  QQ.df <- cbind.data.frame(stat.sort.df, pval.obs.log10, pval.exp.log10)
  return(QQ.df)
}
SNP.functional.TFBS.gene.summary.stat.2 <- QQfunct(SNP.functional.TFBS.gene.summary.stat)

# Prepare plot

SNP.functional.TFBS.plot <- SNP.functional.TFBS.gene.summary.stat.2 %>% mutate(enrichment = case_when(enr >= 1 ~ 1, TRUE ~ 0)) %>%
  mutate(enr.2 = case_when(enrichment == 1 ~ enr, TRUE ~ 1/enr))

pdf("./Rplot/Fig_5D.pdf", width=5, height=4, useDingbats=FALSE)
SNP.functional.TFBS.plot %>% ggplot(aes(x=pval.exp.log10, y=pval.obs.log10, size = enr.2, color = enrichment)) +
  geom_point(alpha=0.7) + scale_size(range = c(1, 20), name="enrichment", breaks = c(2,4,8)) + geom_abline(intercept = 0, slope = 1, size = 0.5, linetype = "dashed") +
  scale_colour_gradient(low = "#0072B2", high = "#DE1C14") + xlab("Expected P value (-log10)") + ylab("Observed P value (-log10)") +
  ggtitle("TFBS enrichment/depletion in genes marked with origins") + xlim(0,3) +
  geom_text(aes(label=ifelse(pval.obs.log10 >= 3,as.character(TF),'')), size=3.5, hjust=2, vjust=0) +
  theme_bw() + theme(aspect.ratio=1, plot.title = element_text(size = 11), axis.title=element_text(size=11))
dev.off()

############################################################################
# Plot seq logo of representative hits

KLF15.PFM <- readJASPARMatrix("./Dataset/JASPAR2020_CORE_non-redundant_pfms_jaspar/MA1513.1.jaspar", matrixClass= "PFM")[[1]]
TFAP2E.PFM <- readJASPARMatrix("./Dataset/JASPAR2020_CORE_non-redundant_pfms_jaspar/MA1569.1.jaspar", matrixClass= "PFM")[[1]]
MAZ.PFM <- readJASPARMatrix("./Dataset/JASPAR2020_CORE_non-redundant_pfms_jaspar/MA1522.1.jaspar", matrixClass= "PFM")[[1]]
NRF1.PFM <- readJASPARMatrix("./Dataset/JASPAR2020_CORE_non-redundant_pfms_jaspar/MA0506.1.jaspar", matrixClass= "PFM")[[1]]
ETV4.PFM <- readJASPARMatrix("./Dataset/JASPAR2020_CORE_non-redundant_pfms_jaspar/MA0764.2.jaspar", matrixClass= "PFM")[[1]]

cs1 <- make_col_scheme(chars=c('A', 'T', 'C', 'G'),cols=c('#009E73', '#DE1C14', '#0072B2', '#F0E442'))

a <- ggplot() + geom_logo(KLF15.PFM@profileMatrix, col_scheme=cs1) + theme_logo() + theme(aspect.ratio=0.4) + ggtitle("KLF15 - Sp/KLF factors")
b <- ggplot() + geom_logo(TFAP2E.PFM@profileMatrix, col_scheme=cs1) + theme_logo() + theme(aspect.ratio=0.4) + ggtitle("TFAP2E - AP2 factors")
c <- ggplot() + geom_logo(MAZ.PFM@profileMatrix, col_scheme=cs1) + theme_logo() + theme(aspect.ratio=0.4) + ggtitle("MAZ - C2H2 zinc finger factors")
d <- ggplot() + geom_logo(NRF1.PFM@profileMatrix, col_scheme=cs1) + theme_logo() + theme(aspect.ratio=0.4) + ggtitle("NRF1 - Basic leucine zipper factors")
e <- ggplot() + geom_logo(ETV4.PFM@profileMatrix, col_scheme=cs1) + theme_logo() + theme(aspect.ratio=0.4) + ggtitle("ETV4 - Tryptophan cluster factor")

pdf("./Rplot/Fig_5E.pdf", width=5, height=10, useDingbats=FALSE)
ggarrange(a,b,c,d,e, ncol = 1, nrow = 5)
dev.off()

############################################################################
# Assess the contribution of replication efficiency to the enrichment of selected TFBS

SNP.functional.TFBS.gene.with.ori.low <- SNP.functional.TFBS %>% filter(SNP.id %in% functional.SNP.closest.gene.start.with.ori.low$SNP.id) # 1,243 genes
SNP.functional.TFBS.gene.with.ori.medium <- SNP.functional.TFBS %>% filter(SNP.id %in% functional.SNP.closest.gene.start.with.ori.medium$SNP.id) # 3,460 genes
SNP.functional.TFBS.gene.with.ori.high <- SNP.functional.TFBS %>% filter(SNP.id %in% functional.SNP.closest.gene.start.with.ori.high$SNP.id) # 5,624 genes

# KLF15
nrow(SNP.functional.TFBS.gene.with.ori.low[which(SNP.functional.TFBS.gene.with.ori.low$TF == "KLF15"),])  # 111/(1243*5000) = 1.786002e-05 per bp
nrow(SNP.functional.TFBS.gene.with.ori.medium[which(SNP.functional.TFBS.gene.with.ori.medium$TF == "KLF15"),])  # 400/(3460*5000) = 2.312139e-05 per bp
nrow(SNP.functional.TFBS.gene.with.ori.high[which(SNP.functional.TFBS.gene.with.ori.high$TF == "KLF15"),])  # 732/(5624*5000) = 2.603129e-05 per bp
chi.table <- rbind(c(111,1243), c(400,3460),c(732,5624))
chisq.test(chi.table)$p.value # 0.0010616110, enr = 1.457517

KLF15.eff.df1 <- tibble(Count = 111/(1243*5000), class = "low")
KLF15.eff.df2 <- tibble(Count = 400/(3460*5000), class = "medium")
KLF15.eff.df3 <- tibble(Count = 732/(5624*5000), class = "high")
KLF15.eff.df <- rbind(KLF15.eff.df1, KLF15.eff.df2, KLF15.eff.df3)
KLF15.plot <- KLF15.eff.df %>% mutate(Class = fct_relevel(class, "low", "medium", "high")) %>%
  ggplot(aes(x=Class, y=Count, fill=Class)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=c("#E69F00", "#D55E00", "#DE1C14")) +
  ggtitle("KLF15\npval = 0.0010616110\nenr = 1.457517") + ylab("Number of functional SNPs per bp") +
  theme_bw() + theme(aspect.ratio=2, plot.title = element_text(size = 11)) + theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 0.5)) + theme(axis.title.x = element_blank()) +
  scale_y_continuous(labels = comma, limits=c(0, 0.00003))

# TFAP2E
nrow(SNP.functional.TFBS.gene.with.ori.low[which(SNP.functional.TFBS.gene.with.ori.low$TF == "TFAP2E"),])  # 98/(1243*5000) = 1.57683e-05 per bp
nrow(SNP.functional.TFBS.gene.with.ori.medium[which(SNP.functional.TFBS.gene.with.ori.medium$TF == "TFAP2E"),])  # 392/(3460*5000) = 2.265896e-05 per bp
nrow(SNP.functional.TFBS.gene.with.ori.high[which(SNP.functional.TFBS.gene.with.ori.high$TF == "TFAP2E"),])  # 725/(5624*5000) = 2.578236e-05 per bp
chi.table <- rbind(c(98,1243), c(392,3460),c(725,5624))
chisq.test(chi.table)$p.value # 3.328952e-05, enr = 1.635075

TFAP2E.eff.df1 <- tibble(Count = 98/(1243*5000), class = "low")
TFAP2E.eff.df2 <- tibble(Count = 392/(3460*5000), class = "medium")
TFAP2E.eff.df3 <- tibble(Count = 725/(5624*5000), class = "high")
TFAP2E.eff.df <- rbind(TFAP2E.eff.df1, TFAP2E.eff.df2, TFAP2E.eff.df3)
TFAP2E.plot <- TFAP2E.eff.df %>% mutate(Class = fct_relevel(class, "low", "medium", "high")) %>%
  ggplot(aes(x=Class, y=Count, fill=Class)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=c("#E69F00", "#D55E00", "#DE1C14")) +
  ggtitle("TFAP2E\npval =3.328952e-05\nenr =1.635075") + ylab("Number of functional SNPs per bp") +
  theme_bw() + theme(aspect.ratio=2, plot.title = element_text(size = 11)) + theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 0.5)) + theme(axis.title.x = element_blank()) +
  scale_y_continuous(labels = comma, limits=c(0, 0.00003))

# MAZ
nrow(SNP.functional.TFBS.gene.with.ori.low[which(SNP.functional.TFBS.gene.with.ori.low$TF == "MAZ"),])  # 112/(1243*5000) = 1.802092e-05 per bp
nrow(SNP.functional.TFBS.gene.with.ori.medium[which(SNP.functional.TFBS.gene.with.ori.medium$TF == "MAZ"),])  # 328/(3460*5000) = 1.895954e-05 per bp
nrow(SNP.functional.TFBS.gene.with.ori.high[which(SNP.functional.TFBS.gene.with.ori.high$TF == "MAZ"),])  # 618/(5624*5000) = 2.197724e-05 per bp
chi.table <- rbind(c(112,1243), c(328,3460),c(618,5624))
chisq.test(chi.table)$p.value # 0.04420774, enr = 1.21954

MAZ.eff.df1 <- tibble(Count = 112/(1243*5000), class = "low")
MAZ.eff.df2 <- tibble(Count = 328/(3460*5000), class = "medium")
MAZ.eff.df3 <- tibble(Count = 618/(5624*5000), class = "high")
MAZ.eff.df <- rbind(MAZ.eff.df1, MAZ.eff.df2, MAZ.eff.df3)
MAZ.plot <- MAZ.eff.df %>% mutate(Class = fct_relevel(class, "low", "medium", "high")) %>%
  ggplot(aes(x=Class, y=Count, fill=Class)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=c("#E69F00", "#D55E00", "#DE1C14")) +
  ggtitle("MAZ\npval = 0.04420774\nenr = 1.21954") + ylab("Number of functional SNPs per bp") +
  theme_bw() + theme(aspect.ratio=2, plot.title = element_text(size = 11)) + theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 0.5)) + theme(axis.title.x = element_blank()) +
  scale_y_continuous(labels = comma, limits=c(0, 0.00003))

# NRF1
nrow(SNP.functional.TFBS.gene.with.ori.low[which(SNP.functional.TFBS.gene.with.ori.low$TF == "NRF1"),])  # 28/(1243*5000) = 4.505229e-06 per bp
nrow(SNP.functional.TFBS.gene.with.ori.medium[which(SNP.functional.TFBS.gene.with.ori.medium$TF == "NRF1"),])  # 72/(3460*5000) = 4.16185e-06 per bp
nrow(SNP.functional.TFBS.gene.with.ori.high[which(SNP.functional.TFBS.gene.with.ori.high$TF == "NRF1"),])  # 185/(5624*5000) = 6.578947e-06 per bp
chi.table <- rbind(c(28,1243), c(72,3460),c(185,5624))
chisq.test(chi.table)$p.value # 0.002103182, enr = 1.460291

NRF1.eff.df1 <- tibble(Count = 28/(1243*5000), class = "low")
NRF1.eff.df2 <- tibble(Count = 72/(3460*5000), class = "medium")
NRF1.eff.df3 <- tibble(Count = 185/(5624*5000), class = "high")
NRF1.eff.df <- rbind(NRF1.eff.df1, NRF1.eff.df2, NRF1.eff.df3)
NRF1.plot <- NRF1.eff.df %>% mutate(Class = fct_relevel(class, "low", "medium", "high")) %>%
  ggplot(aes(x=Class, y=Count, fill=Class)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=c("#E69F00", "#D55E00", "#DE1C14")) +
  ggtitle("NRF1\npval = 0.002103182\nenr = 1.460291") + ylab("Number of functional SNPs per bp") +
  theme_bw() + theme(aspect.ratio=2, plot.title = element_text(size = 11)) + theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 0.5)) + theme(axis.title.x = element_blank()) +
  scale_y_continuous(labels = comma, limits=c(0, 0.00001))

# ETV4
nrow(SNP.functional.TFBS.gene.with.ori.low[which(SNP.functional.TFBS.gene.with.ori.low$TF == "ETV4"),])  # 82/(1243*5000) = 1.319389e-05 per bp
nrow(SNP.functional.TFBS.gene.with.ori.medium[which(SNP.functional.TFBS.gene.with.ori.medium$TF == "ETV4"),])  # 253/(3460*5000) = 1.484512e-05 per bp
nrow(SNP.functional.TFBS.gene.with.ori.high[which(SNP.functional.TFBS.gene.with.ori.high$TF == "ETV4"),])  # 501/(5624*5000) = 1.78165e-05 per bp
chi.table <- rbind(c(82,1243), c(253,3460),c(501,5624))
chisq.test(chi.table)$p.value # 0.007109866, enr = 1.35036

ETV4.eff.df1 <- tibble(Count = 82/(1243*5000), class = "low")
ETV4.eff.df2 <- tibble(Count = 253/(3460*5000), class = "medium")
ETV4.eff.df3 <- tibble(Count = 501/(5624*5000), class = "high")
ETV4.eff.df <- rbind(ETV4.eff.df1, ETV4.eff.df2, ETV4.eff.df3)
ETV4.plot <- ETV4.eff.df %>% mutate(Class = fct_relevel(class, "low", "medium", "high")) %>%
  ggplot(aes(x=Class, y=Count, fill=Class)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=c("#E69F00", "#D55E00", "#DE1C14")) +
  ggtitle("ETV4\npval = 0.007109866\nenr = 1.35036") + ylab("Number of functional SNPs per bp") +
  theme_bw() + theme(aspect.ratio=2, plot.title = element_text(size = 11)) + theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 0.5)) + theme(axis.title.x = element_blank()) +
  scale_y_continuous(labels = comma, limits=c(0, 0.000025))

# Prepare plot reporting the influence of origin efficiency on the density of functional SNPs overlapping TFBS

pdf("./Rplot/Fig_5E_2.pdf", width=11, height=5, useDingbats=FALSE)
ggarrange(KLF15.plot,TFAP2E.plot,MAZ.plot,NRF1.plot,ETV4.plot, ncol = 5, nrow = 1)
dev.off()

############################################################################
############################################################################
###                  Quantitative trait loci analysis                    ###
############################################################################
############################################################################

############################################################################
# Single tissue cis-eQTL and cis-sQTL data
# downloaded from the GTEx portal (https://www.gtexportal.org/home/datasets)
# GTEx_Analysis_v8_eQTL.tar and GTEx_Analysis_v8_eQTL.tar

# eGene and significant variant-gene associations based on permutations.
# The archive contains a *.egenes.txt.gz and *.signif_variant_gene_pairs.txt.gz file for each tissue.
# Note that the *.egenes.txt.gz files contain data for all genes tested;
# to obtain the list of eGenes, select the rows with 'qval' ≤ 0.05.
# README file describing the contents of the eQTL and sQTL files can be found at the same URL
# README_eQTL_v8.txt

# sGene and significant variant-gene associations are based on LeafCutter intron excision phenotypes

############################################################################
# eQTLs analysis

# Assess the fraction of eGenes containing an origin within 2.5kb of their promoters

# Load gene lists
genes.TSS.with.ori <- read.csv("./Figure_4/genes.TSS.with.ori.bed") %>% dplyr::select(-X)
genes.TSS.without.ori <- read.csv("./Figure_4/genes.TSS.without.ori.bed") %>% dplyr::select(-X)

# Select eGenes analysis files 
GTEx.files <- list.files("Dataset/GTEx_Analysis_v8_eQTL/")
# Recover list of tissues
GTEx.files.names <- GTEx.files[grep("*.v8.egenes.txt.gz", GTEx.files)]
GTEx.tissues <- sub(".v8.egenes.txt.gz", "", GTEx.files.names) # 49 tissues
# Compute the number of eGenes with or without origins
eQTL.gene.analysis <- tibble(tissue = character(), eGene.with.ori = numeric(), eGene.without.ori = numeric())
for (i in GTEx.tissues) {
  print(i)
  tissue.eQTL <- read.table(gzfile(paste("./Dataset/GTEx_Analysis_v8_eQTL/", i, ".v8.egenes.txt.gz", sep = "")), sep="\t", header = T)
  tissue.eQTL.sig.genes <- tissue.eQTL %>% dplyr::filter(qval <= 0.05 & tss_distance >= -2500 & tss_distance <= 2500) %>% tidyr::separate(gene_id, c("gene.id", "version"))
  tissue.eQTL.sig.genes.with.ori <- tissue.eQTL.sig.genes %>% dplyr::filter(gene.id %in% genes.TSS.with.ori$gene.id)
  tissue.eQTL.sig.genes.without.ori <- tissue.eQTL.sig.genes %>% dplyr::filter(gene.id %in% genes.TSS.without.ori$gene.id)
  p <- tibble(tissue = i, eGene.with.ori = nrow(tissue.eQTL.sig.genes.with.ori), eGene.without.ori = nrow(tissue.eQTL.sig.genes.without.ori))
  eQTL.gene.analysis <- rbind(eQTL.gene.analysis,p)
}
# Compute enrichment and associated pvalues (2x2 chi-squared test) to test the association for a gene to be a egene and the presence of origins
pval.log10 <- vector()
enr <- vector()
for (i in 1:nrow(eQTL.gene.analysis)) {
  tissue <- eQTL.gene.analysis[i,]
  nbr.gene <- tissue$eGene.with.ori + tissue$eGene.without.ori
  chi.table.ori <- rbind(c(tissue$eGene.with.ori,nbr.gene), c(tissue$eGene.without.ori,nbr.gene))
  p <- -log10(chisq.test(chi.table.ori)$p.value)
  pval.log10 <- append(pval.log10,p)
  q <- (tissue$eGene.with.ori/nbr.gene)/(tissue$eGene.without.ori/nbr.gene)
  enr <- append(enr,q)
}
eQTL.gene.analysis.summary.stat <- cbind.data.frame(eQTL.gene.analysis, pval.log10, enr)
# Prepare bar plot comparing eQTL density at genes with ori or not
eQTL.density <- eQTL.gene.analysis.summary.stat %>% mutate(eQTL.density.with.ori = eGene.with.ori/(nrow(genes.TSS.with.ori)*5000),
                                                           eQTL.density.without.ori = eGene.without.ori/(nrow(genes.TSS.without.ori)*5000))
# Rewrite tissue name
tissue.2 <- str_replace_all(eQTL.density$tissue, "_", " ")
tissue.2 <- str_replace_all(tissue.2, " c-1", "")
tissue.2 <- str_replace_all(tissue.2, " BA24", "")
tissue.2 <- str_replace_all(tissue.2, " BA9", "")
eQTL.density$tissue <- tissue.2
#
eQTL.density.ori.plot.1 <- eQTL.density %>% mutate(enr = round(enr, digits = 2), pval.log10 = scientific(10^-pval.log10, digits = 2), info = paste(enr, " (", pval.log10, ")", sep = "")) %>% dplyr::select(tissue, density = eQTL.density.with.ori, info) %>% mutate(class = "with ori")
eQTL.density.ori.plot.2 <- eQTL.density %>% dplyr::select(tissue, density = eQTL.density.without.ori) %>% mutate(info = "", class = "without ori")
eQTL.density.ori.plot <- rbind(eQTL.density.ori.plot.1, eQTL.density.ori.plot.2) %>% mutate(exp = "eQTL")
# Save plot
pdf("./Rplot/Fig_S6_A.pdf", width=5, height=7, useDingbats=FALSE)
eQTL.density.ori.plot %>% mutate(tissue = fct_reorder(tissue, density)) %>% 
  ggplot(aes(x=tissue, y=density, fill = class)) +
  scale_fill_manual(values=c("#DE1C14", "#D4D2D2")) + 
  geom_bar(position="dodge", stat = "identity") + scale_y_continuous(labels = comma, limits = c(0,0.000035)) + ylab("Number of eQTLs per bp") +
  coord_flip() + theme_bw() + theme(aspect.ratio=5, axis.title.y=element_blank(), axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5)) + geom_text(aes(label=info), hjust=-0.1,  size=1.5)
dev.off()  

############################################################################
# sQTLs analysis

# Assess the fraction of sGenes containing an origin within 2.5kb of an intron/exon or exon/inron junction

# Load splice sites lists
splice.sites.with.ori <- read.csv("./Figure_4/splice.sites.with.ori.bed") %>% dplyr::select(-X)
splice.sites.without.ori <- read.csv("./Figure_4/splice.sites.without.ori.bed") %>% dplyr::select(-X)

# Recover ids of SNPs associated with junctions

population.SNP.closest.splice.site.10kb  <- read.csv("./Figure_4/population.SNP.closest.splice.site.10kb.bed", sep = '\t')
junctions.SNP <- population.SNP.closest.splice.site.10kb %>% filter(SNP.gene.dist >= -2500 & SNP.gene.dist <= 2500) %>% mutate(junction.id = paste(chr, splice.site.pos, sep = "_")) %>% dplyr::select(rs_id_dbSNP151_GRCh38p7 = SNP.id, dist = SNP.gene.dist, junction.id)
# 11,637,314 SNPs

# Recover junction ids
junction.with.ori <- unique(c(splice.sites.with.ori$junction.id)) # 94,821 junctions
junction.without.ori <- unique(c(splice.sites.without.ori$junction.id)) # 314,293 junctions

# Select eGenes analysis files 
GTEx.files <- list.files("Dataset/GTEx_Analysis_v8_sQTL/")
# Recover list of tissues
GTEx.files.names <- GTEx.files[grep("*.v8.sgenes.txt.gz", GTEx.files)]
GTEx.tissues <- sub(".v8.sgenes.txt.gz", "", GTEx.files.names) # 49 tissues
# Compute the number of eGenes with or without origins
sQTL.gene.analysis <- tibble(tissue = character(), sGene.with.ori = numeric(), sGene.without.ori = numeric())
for (i in GTEx.tissues) {
  print(i)
  tissue.sQTL <- read.table(gzfile(paste("./Dataset/GTEx_Analysis_v8_sQTL/", i, ".v8.sgenes.txt.gz", sep = "")), sep="\t", header = T)
  tissue.sQTL.sig.genes <- tissue.sQTL %>% dplyr::filter(qval <= 0.05) %>% dplyr::filter(rs_id_dbSNP151_GRCh38p7 %in% junctions.SNP$rs_id_dbSNP151_GRCh38p7) %>% tidyr::separate(gene_id, c("gene.id", "version")) %>% left_join(junctions.SNP, by = "rs_id_dbSNP151_GRCh38p7")
  tissue.sQTL.sig.genes.with.ori <- tissue.sQTL.sig.genes %>% dplyr::filter(junction.id %in% junction.with.ori)
  tissue.sQTL.sig.genes.without.ori <- tissue.sQTL.sig.genes %>% dplyr::filter(junction.id %in% junction.without.ori)
  p <- tibble(tissue = i, sGene.with.ori = nrow(tissue.sQTL.sig.genes.with.ori), sGene.without.ori = nrow(tissue.sQTL.sig.genes.without.ori))
  sQTL.gene.analysis <- rbind(sQTL.gene.analysis,p)
}
# Compute enrichment and associated pvalues (2x2 chi-squared test) to test the association for a gene to be a sgene and the presence of origins
pval.log10 <- vector()
enr <- vector()
for (i in 1:nrow(sQTL.gene.analysis)) {
  tissue <- sQTL.gene.analysis[i,]
  chi.table.ori <- rbind(c(tissue$sGene.with.ori,length(junction.with.ori)), c(tissue$sGene.without.ori,length(junction.without.ori)))
  p <- -log10(chisq.test(chi.table.ori)$p.value)
  pval.log10 <- append(pval.log10,p)
  q <- (tissue$sGene.with.ori/length(junction.with.ori))/(tissue$sGene.without.ori/length(junction.without.ori))
  enr <- append(enr,q)
}
sQTL.gene.analysis.summary.stat <- cbind.data.frame(sQTL.gene.analysis, pval.log10, enr)
# Prepare bar plot comparing eQTL density at genes with ori or not
sQTL.density <- sQTL.gene.analysis.summary.stat %>% mutate(sQTL.density.with.ori = sGene.with.ori/(length(junction.with.ori)*5000),
                                                           sQTL.density.without.ori = sGene.without.ori/(length(junction.without.ori)*5000))
# Rewrite tissue name
tissue.2 <- str_replace_all(sQTL.density$tissue, "_", " ")
tissue.2 <- str_replace_all(tissue.2, " c-1", "")
tissue.2 <- str_replace_all(tissue.2, " BA24", "")
tissue.2 <- str_replace_all(tissue.2, " BA9", "")
sQTL.density$tissue <- tissue.2
#
sQTL.density.ori.plot.1 <- sQTL.density %>% mutate(enr = round(enr, digits = 2), pval.log10 = scientific(10^-pval.log10, digits = 2), info = paste(enr, " (", pval.log10, ")", sep = "")) %>% dplyr::select(tissue, density = sQTL.density.with.ori, info) %>% mutate(class = "with ori")
sQTL.density.ori.plot.2 <- sQTL.density %>% dplyr::select(tissue, density = sQTL.density.without.ori) %>% mutate(info = "", class = "without ori")
sQTL.density.ori.plot <- rbind(sQTL.density.ori.plot.1, sQTL.density.ori.plot.2) %>% mutate(exp = "sQTL")
# Save plot
pdf("./Rplot/Fig_S6_A_2.pdf", width=5, height=7, useDingbats=FALSE)
sQTL.density.ori.plot %>% mutate(tissue = fct_reorder(tissue, density)) %>% 
  ggplot(aes(x=tissue, y=density, fill = class)) +
  scale_fill_manual(values=c("#009E73", "#D4D2D2")) + 
  geom_bar(position="dodge", stat = "identity") + scale_y_continuous(labels = comma, limits = c(0,0.000010)) + ylab("Number of sQTLs per bp") +
  coord_flip() + theme_bw() + theme(aspect.ratio=5, axis.title.y=element_blank(), axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5)) + geom_text(aes(label=info), hjust=-0.1,  size=1.5)
dev.off()  

############################################################################
# Combine results and prepare selection plot

QTL.ori <- rbind(eQTL.density.ori.plot, sQTL.density.ori.plot)

# Select 10 tissues with data for all conditions 
tissue.selected <- c("Testis", "Thyroid", "Muscle Skeletal", "Esophagus Mucosa", "Whole Blood", "Artery Aorta", "Lung", "Colon Sigmoid",  "Pancreas", "Heart Left Ventricle",
                     "Spleen", "Stomach", "Prostate", "Ovary", "Minor Salivary Gland", "Liver", "Uterus" , "Vagina", "Brain Amygdala", "Kidney Cortex")

# Prepare plots
eQTL.plot <- QTL.ori %>% filter(tissue %in% tissue.selected & exp == "eQTL")
sQTL.plot <- QTL.ori %>% filter(tissue %in% tissue.selected & exp == "sQTL")

# Plot
eQTL.plot.2 <- eQTL.plot %>% mutate(Tissue = factor(tissue, rev(tissue.selected))) %>% 
  ggplot(aes(x=Tissue, y=density, fill = class)) +
  scale_fill_manual(values=c("#DE1C14", "#D4D2D2")) + 
  geom_bar(position="dodge", stat = "identity") + scale_y_continuous(labels = comma, limits = c(0,0.00003)) + ylab("Number of eQTLs per bp") +
  coord_flip() + theme_bw() + theme(aspect.ratio=2, axis.title.y=element_blank(), axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5)) + geom_text(aes(label=info), hjust=-0.1,  size=3)
sQTL.plot.2 <- sQTL.plot %>% mutate(Tissue = factor(tissue, rev(tissue.selected))) %>% 
  ggplot(aes(x=Tissue, y=density, fill = class)) +
  scale_fill_manual(values=c("#009E73", "#D4D2D2")) + 
  geom_bar(position="dodge", stat = "identity") + scale_y_continuous(labels = comma, limits = c(0,0.000009)) + ylab("Number of sQTLs per bp") +
  coord_flip() + theme_bw() + theme(aspect.ratio=2, axis.title.y=element_blank(), axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5)) + geom_text(aes(label=info), hjust=-0.1,  size=3)

# Save plots
pdf("./Rplot/Fig_6.pdf", width=10, height=5, useDingbats=FALSE)
ggarrange(eQTL.plot.2, sQTL.plot.2, nrow = 1 , ncol = 2)
dev.off()






