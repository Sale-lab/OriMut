
# P. Murat, MRC Laboratory of Molecular Biology, July 2022
# Code for generating Figure 4 with associated supporting Figures
# in Murat et al. DNA replication initiation shapes the mutational landscape and expression of the human genome

library(dplyr)
library(biomaRt)
library(ggpubr)
library(tidyr)
library(GenomicRanges)
library(rtracklayer)
library(forcats)

setwd("/Volumes/pmurat/OriMut")

############################################################################
############################################################################
###               Overlap quantification between origins                 ###
###                         and gene features                            ###
############################################################################
############################################################################

############################################################################
# Load origins coordinates

SNS.core <- read.table("./Figure_1/Core_origins_unique.bed") # only SNS core origin non overlapping with IniSeq origins are considered
nrow(SNS.core)  # 31,523 origins
colnames(SNS.core) <- c("chr", "start", "end", "ID")
SNS.core <- SNS.core %>% mutate(chr = str_remove(chr, "chr"))
write.table(SNS.core, "./Figure_4/SNS.core.ori.bed", sep="\t", col.names = F, row.names = F, quote = F)

SNS.stochastic <- read.table("./Dataset/GSE128477_Stochastic_origins_hg38.bed")
nrow(SNS.stochastic)  # 256,600 origins
colnames(SNS.stochastic) <- c("chr", "start", "end", "ID")
SNS.stochastic <- SNS.stochastic %>% mutate(chr = str_remove(chr, "chr"))
write.table(SNS.stochastic, "./Figure_4/SNS.stochastic.ori.bed", sep="\t", col.names = F, row.names = F, quote = F)

IniSeq.ori <- read.table("./Dataset/GSM5658908_Ini-seq2.called.replication.origins.bed", header = F)
colnames(IniSeq.ori) <- c("chr", "start", "end", "EFF")
nrow(IniSeq.ori) # 23,905 origins
write.table(IniSeq.ori, "./Figure_4/IniSeq.ori.bed", sep="\t", col.names = F, row.names = F, quote = F)

# Prepare bed files reporting the center of the origins

SNS.core.center <- SNS.core %>% mutate(start = round((start+end)/2), end = start+1) %>%  filter(start > 0 & end > 0)
SNS.stochastic.center <- SNS.stochastic %>% mutate(start = round((start+end)/2), end = start+1) %>%  filter(start > 0 & end > 0)
IniSeq.center <- IniSeq.ori %>% mutate(start = round((start+end)/2), end = start+1) %>%  filter(start > 0 & end > 0)

write.table(SNS.core.center, "./Figure_4/SNS.core.ori.center.bed", sep="\t", col.names = F, row.names = F, quote = F)
write.table(SNS.stochastic.center, "./Figure_4/SNS.stochastic.ori.center.bed", sep="\t", col.names = F, row.names = F, quote = F)
write.table(IniSeq.center, "./Figure_4/IniSeq.ori.center.bed", sep="\t", col.names = F, row.names = F, quote = F)

############################################################################
# Recover gene coordinates

# Recover all genes coordinates from Biomart
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
att <- listAttributes(mart)
# Select only protein coding genes
hsapiens.genes <-  getBM(attributes=c("ensembl_gene_id", "hgnc_symbol", "chromosome_name", "start_position", "end_position", "strand", "description"), filters='biotype', values=c('protein_coding'), mart=mart)
# Remove novel protein/transcript and unusual chromosomes
my.chr <- c(1:22, "X", "Y")
hsapiens.genes.2 <- hsapiens.genes %>% filter(hgnc_symbol != "" & chromosome_name %in% my.chr )  # 19,335 unique genes
# Prepare bed file reporting gene coordinates
hsapiens.genes.coordinate <- hsapiens.genes.2 %>% mutate(start = as.integer(start_position), end = as.integer(end_position)) %>%
  dplyr::select(chr = chromosome_name, start, end, strand, ensembl_gene_id,hgnc_symbol) %>% filter(start > 0 & end > 0)
# Prepare bed file reporting the start of the genes
hsapiens.genes.start <- hsapiens.genes.2 %>% mutate(start = ifelse(strand == 1, as.integer(start_position), as.integer(end_position))) %>% mutate(end = as.integer(start +1)) %>% 
  dplyr::select(chr = chromosome_name, start, end, strand, ensembl_gene_id,hgnc_symbol) %>% filter(start > 0 & end > 0)
# Prepare bed file reporting gene promoters (-1.0 kb, downstream)
hsapiens.promoter <- hsapiens.genes.2 %>% mutate(start = ifelse(strand == 1, as.integer(start_position-1000), as.integer(end_position))) %>%
  mutate(end = ifelse(strand == 1, as.integer(start_position), as.integer(end_position+1000))) %>% 
  dplyr::select(chr = chromosome_name, start, end, strand, ensembl_gene_id,hgnc_symbol) %>% filter(start > 0 & end > 0)

# Save bed files
write.table(hsapiens.genes.coordinate, "./Figure_4/hsapiens.genes.coordinate.bed", sep="\t", col.names = F, row.names = F, quote = F)
write.table(hsapiens.genes.start, "./Figure_4/hsapiens.genes.start.bed", sep="\t", col.names = F, row.names = F, quote = F)
write.table(hsapiens.promoter, "./Figure_4/hsapiens.promoters.1kb.bed", sep="\t", col.names = F, row.names = F, quote = F)

# Recover all exon coordinates from Biomart
# Select only gencode basic transcripts
hsapiens.exons <-  getBM(attributes=c("ensembl_transcript_id", "ensembl_gene_id", "chromosome_name", "exon_chrom_start", "exon_chrom_end", "strand", "description"), filters=c("transcript_gencode_basic"), values=c(TRUE), mart=mart)
# Remove unusual chromosomes
my.chr <- c(1:22, "X", "Y")
hsapiens.exons.2 <- hsapiens.exons %>% filter(chromosome_name %in% my.chr )   # 60,570 unique genes, 110,877 transcripts
# Select only protein coding genes
hsapiens.exons.3 <- hsapiens.exons.2 %>% filter(ensembl_gene_id %in% hsapiens.genes.2$ensembl_gene_id) # 19,334 unique genes, 59,563 transcripts
# Prepare bed file reporting exon coordinates
hsapiens.exons.coordinate <- hsapiens.exons.3 %>% mutate(start = exon_chrom_start, end = exon_chrom_end) %>%
  dplyr::select(chr = chromosome_name, start, end, strand, ensembl_transcript_id, ensembl_gene_id) %>% filter(start > 0 & end > 0)

# Recover all UTR5 coordinates from Biomart
# Select only gencode basic transcripts
hsapiens.UTR5 <-  getBM(attributes=c("ensembl_transcript_id", "ensembl_gene_id", "chromosome_name", "5_utr_start", "5_utr_end"), filters=c("transcript_gencode_basic"), values=c(TRUE), mart=mart)
# Remove unusual chromosomes
my.chr <- c(1:22, "X", "Y")
hsapiens.UTR5.2 <- hsapiens.UTR5 %>% filter(chromosome_name %in% my.chr )   # 60,570 unique genes, 110,877 transcripts
# Select only protein coding genes
hsapiens.UTR5.3 <- hsapiens.UTR5.2 %>% filter(ensembl_gene_id %in% hsapiens.genes.2$ensembl_gene_id) # 19,334 unique genes, 59,563 transcripts
# Remove entries with no information (NA values)
hsapiens.UTR5.4 <- hsapiens.UTR5.3 %>% filter(`5_utr_start` >= 0 & `5_utr_end` >= 0) %>%  # 18,929 unique genes, 56,691 transcripts
  dplyr::select(chromosome_name, `5_utr_start`, `5_utr_end`, ensembl_transcript_id, ensembl_gene_id)

# Recover all UTR3 coordinates from Biomart
# Select only gencode basic transcripts
hsapiens.UTR3 <-  getBM(attributes=c("ensembl_transcript_id", "ensembl_gene_id", "chromosome_name", "3_utr_start", "3_utr_end"), filters=c("transcript_gencode_basic"), values=c(TRUE), mart=mart)
# Remove unusual chromosomes
my.chr <- c(1:22, "X", "Y")
hsapiens.UTR3.2 <- hsapiens.UTR3 %>% filter(chromosome_name %in% my.chr )   # 60,570 unique genes, 110,877 transcripts
# Select only protein coding genes
hsapiens.UTR3.3 <- hsapiens.UTR3.2 %>% filter(ensembl_gene_id %in% hsapiens.genes.2$ensembl_gene_id) # 19,334 unique genes, 59,563 transcripts
# Remove entries with no information (NA values)
hsapiens.UTR3.4 <- hsapiens.UTR3.3 %>% filter(`3_utr_start` >= 0 & `3_utr_end` >= 0) %>%   # 19,002 unique genes, 57,092 transcripts
  dplyr::select(chromosome_name, `3_utr_start`, `3_utr_end`, ensembl_transcript_id, ensembl_gene_id)

# Save bed files
write.table(hsapiens.exons.coordinate, "./Figure_4/hsapiens.exons.UTRs.coordinate.bed", sep="\t", col.names = F, row.names = F, quote = F)
write.table(hsapiens.UTR5.4, "./Figure_4/hsapiens.UTR5.coordinate.bed", sep="\t", col.names = F, row.names = F, quote = F)
write.table(hsapiens.UTR3.4, "./Figure_4/hsapiens.UTR3.coordinate.bed", sep="\t", col.names = F, row.names = F, quote = F)

# Use bedtools substract to recover intron coordinates
# bedtools subtract -a ./Figure_4/hsapiens.genes.coordinate.bed -b ./Figure_4/hsapiens.exons.UTRs.coordinate.bed > ./Figure_4/hsapiens.introns.coordinate.bed

# Use bedtools substract to recover exon coordinates (without UTRs)
# bedtools subtract -a ./Figure_4/hsapiens.exons.UTRs.coordinate.bed -b ./Figure_4/hsapiens.UTR5.coordinate.bed > ./Figure_4/hsapiens.exons.no.UTR5.coordinate.bed
# bedtools subtract -a ./Figure_4/hsapiens.exons.no.UTR5.coordinate.bed -b ./Figure_4/hsapiens.UTR3.coordinate.bed > ./Figure_4/hsapiens.exons.no.UTRs.coordinate.bed

############################################################################
# Compute overlaps with bedtools closest function

# sort -k1,1 -k2,2n ./Figure_4/SNS.core.ori.center.bed > ./Figure_4/SNS.core.ori.center.sorted.bed
# sort -k1,1 -k2,2n ./Figure_4/SNS.stochastic.ori.center.bed > ./Figure_4/SNS.stochastic.ori.center.sorted.bed
# sort -k1,1 -k2,2n ./Figure_4/IniSeq.ori.center.bed > ./Figure_4/IniSeq.ori.center.sorted.bed
# sort -k1,1 -k2,2n ./Figure_4/hsapiens.UTR5.coordinate.bed > ./Figure_4/hsapiens.UTR5.coordinate.sorted.bed
# sort -k1,1 -k2,2n ./Figure_4/hsapiens.introns.coordinate.bed > ./Figure_4/hsapiens.introns.coordinate.sorted.bed
# sort -k1,1 -k2,2n ./Figure_4/hsapiens.UTR3.coordinate.bed > ./Figure_4/hsapiens.UTR3.coordinate.sorted.bed
# sort -k1,1 -k2,2n ./Figure_4/hsapiens.genes.coordinate.bed > ./Figure_4/hsapiens.genes.coordinate.sorted.bed

# Analse results

# For IniSeq origins
# bedtools closest -a ./Figure_4/IniSeq.ori.center.sorted.bed -b ./Figure_4/hsapiens.promoters.1kb.sorted.bed -d > ./Figure_4/IniSeq.ori.center.promoters.overlap.bed
IniSeq.ori.center.promoters.overlap <- read.csv("./Figure_4/IniSeq.ori.center.promoters.overlap.bed", header = F, sep = "\t") %>% distinct(V1, V2, V3, V4, .keep_all = T)
nrow(IniSeq.ori.center.promoters.overlap[which(IniSeq.ori.center.promoters.overlap$V11 == 0),])   # 2,216 origins
# bedtools closest -a ./Figure_4/IniSeq.ori.center.sorted.bed -b ./Figure_4/hsapiens.UTR5.coordinate.sorted.bed -d > ./Figure_4/IniSeq.ori.center.UTR5.overlap.bed
IniSeq.ori.center.UTR5.overlap <- read.csv("./Figure_4/IniSeq.ori.center.UTR5.overlap.bed", header = F, sep = "\t") %>% distinct(V1, V2, V3, V4, .keep_all = T)
nrow(IniSeq.ori.center.UTR5.overlap[which(IniSeq.ori.center.UTR5.overlap$V10 == 0),])   # 3,897 origins
# bedtools closest -a ./Figure_4/IniSeq.ori.center.sorted.bed -b ./Figure_4/hsapiens.introns.coordinate.sorted.bed -d > ./Figure_4/IniSeq.ori.center.introns.overlap.bed
IniSeq.ori.center.introns.overlap <- read.csv("./Figure_4/IniSeq.ori.center.introns.overlap.bed", header = F, sep = "\t") %>% distinct(V1, V2, V3, V4, .keep_all = T)
nrow(IniSeq.ori.center.introns.overlap[which(IniSeq.ori.center.introns.overlap$V11 == 0),])   # 9,741 origins
# bedtools closest -a ./Figure_4/IniSeq.ori.center.sorted.bed -b ./Figure_4/hsapiens.UTR3.coordinate.sorted.bed -d > ./Figure_4/IniSeq.ori.center.UTR3.overlap.bed
IniSeq.ori.center.UTR3.overlap <- read.csv("./Figure_4/IniSeq.ori.center.UTR3.overlap.bed", header = F, sep = "\t") %>% distinct(V1, V2, V3, V4, .keep_all = T)
nrow(IniSeq.ori.center.UTR3.overlap[which(IniSeq.ori.center.UTR3.overlap$V10 == 0),])   # 526 origins
# bedtools closest -a ./Figure_4/IniSeq.ori.center.sorted.bed -b ./Figure_4/hsapiens.genes.coordinate.sorted.bed -d > ./Figure_4/IniSeq.ori.center.genes.overlap.bed
IniSeq.ori.center.genes.overlap <- read.csv("./Figure_4/IniSeq.ori.center.genes.overlap.bed", header = F, sep = "\t") %>% distinct(V1, V2, V3, V4, .keep_all = T)
nrow(IniSeq.ori.center.genes.overlap[which(IniSeq.ori.center.genes.overlap$V11 == 0),])   # 17,428 origins

# IniSeq origins summary (23,905 origins)
# promoter  UTR5   intron     UTR3      exon      intergenic
# 2216      3897   9741       526       3264      4261

# For SNS core origins
# bedtools closest -a ./Figure_4/SNS.core.ori.center.sorted.bed -b ./Figure_4/hsapiens.promoters.1kb.sorted.bed -d > ./Figure_4/SNS.core.ori.center.promoters.overlap.bed
SNS.core.ori.center.promoters.overlap <- read.csv("./Figure_4/SNS.core.ori.center.promoters.overlap.bed", header = F, sep = "\t") %>% distinct(V1, V2, V3, V4, .keep_all = T)
nrow(SNS.core.ori.center.promoters.overlap[which(SNS.core.ori.center.promoters.overlap$V11 == 0),])   # 564 origins
# bedtools closest -a ./Figure_4/SNS.core.ori.center.sorted.bed -b ./Figure_4/hsapiens.UTR5.coordinate.sorted.bed -d > ./Figure_4/SNS.core.ori.center.UTR5.overlap.bed
SNS.core.ori.center.UTR5.overlap <- read.csv("./Figure_4/SNS.core.ori.center.UTR5.overlap.bed", header = F, sep = "\t") %>% distinct(V1, V2, V3, V4, .keep_all = T)
nrow(SNS.core.ori.center.UTR5.overlap[which(SNS.core.ori.center.UTR5.overlap$V10 == 0),])   # 166 origins
# bedtools closest -a ./Figure_4/SNS.core.ori.center.sorted.bed -b ./Figure_4/hsapiens.introns.coordinate.sorted.bed -d > ./Figure_4/SNS.core.ori.center.introns.overlap.bed
SNS.core.ori.center.introns.overlap <- read.csv("./Figure_4/SNS.core.ori.center.introns.overlap.bed", header = F, sep = "\t") %>% distinct(V1, V2, V3, V4, .keep_all = T)
nrow(SNS.core.ori.center.introns.overlap[which(SNS.core.ori.center.introns.overlap$V11 == 0),])   # 13,206 origins
# bedtools closest -a ./Figure_4/SNS.core.ori.center.sorted.bed -b ./Figure_4/hsapiens.UTR3.coordinate.sorted.bed -d > ./Figure_4/SNS.core.ori.center.UTR3.overlap.bed
SNS.core.ori.center.UTR3.overlap <- read.csv("./Figure_4/SNS.core.ori.center.UTR3.overlap.bed", header = F, sep = "\t") %>% distinct(V1, V2, V3, V4, .keep_all = T)
nrow(SNS.core.ori.center.UTR3.overlap[which(SNS.core.ori.center.UTR3.overlap$V10 == 0),])   # 373 origins
# bedtools closest -a ./Figure_4/SNS.core.ori.center.sorted.bed -b ./Figure_4/hsapiens.genes.coordinate.sorted.bed -d > ./Figure_4/SNS.core.ori.center.genes.overlap.bed
SNS.core.ori.center.genes.overlap <- read.csv("./Figure_4/SNS.core.ori.center.genes.overlap.bed", header = F, sep = "\t") %>% distinct(V1, V2, V3, V4, .keep_all = T)
nrow(SNS.core.ori.center.genes.overlap[which(SNS.core.ori.center.genes.overlap$V11 == 0),])   # 14,404 origins

# SNS core origins summary (31,523 origins)
# promoter  UTR5   intron     UTR3      exon      intergenic
# 564       166    13206      373       659       12294

# For SNS stochastic origins
# bedtools closest -a ./Figure_4/SNS.stochastic.ori.center.sorted.bed -b ./Figure_4/hsapiens.promoters.1kb.sorted.bed -d > ./Figure_4/SNS.stochastic.ori.center.promoters.overlap.bed
SNS.stochastic.ori.center.promoters.overlap <- read.csv("./Figure_4/SNS.stochastic.ori.center.promoters.overlap.bed", header = F, sep = "\t") %>% distinct(V1, V2, V3, V4, .keep_all = T)
nrow(SNS.stochastic.ori.center.promoters.overlap[which(SNS.stochastic.ori.center.promoters.overlap$V11 == 0),])   # 1579 origins
# bedtools closest -a ./Figure_4/SNS.stochastic.ori.center.sorted.bed -b ./Figure_4/hsapiens.UTR5.coordinate.sorted.bed -d > ./Figure_4/SNS.stochastic.ori.center.UTR5.overlap.bed
SNS.stochastic.ori.center.UTR5.overlap <- read.csv("./Figure_4/SNS.stochastic.ori.center.UTR5.overlap.bed", header = F, sep = "\t") %>% distinct(V1, V2, V3, V4, .keep_all = T)
nrow(SNS.stochastic.ori.center.UTR5.overlap[which(SNS.stochastic.ori.center.UTR5.overlap$V10 == 0),])   # 554 origins
# bedtools closest -a ./Figure_4/SNS.stochastic.ori.center.sorted.bed -b ./Figure_4/hsapiens.introns.coordinate.sorted.bed -d > ./Figure_4/SNS.stochastic.ori.center.introns.overlap.bed
SNS.stochastic.ori.center.introns.overlap <- read.csv("./Figure_4/SNS.stochastic.ori.center.introns.overlap.bed", header = F, sep = "\t") %>% distinct(V1, V2, V3, V4, .keep_all = T)
nrow(SNS.stochastic.ori.center.introns.overlap[which(SNS.stochastic.ori.center.introns.overlap$V11 == 0),])   # 34,006 origins
# bedtools closest -a ./Figure_4/SNS.stochastic.ori.center.sorted.bed -b ./Figure_4/hsapiens.UTR3.coordinate.sorted.bed -d > ./Figure_4/SNS.stochastic.ori.center.UTR3.overlap.bed
SNS.stochastic.ori.center.UTR3.overlap <- read.csv("./Figure_4/SNS.stochastic.ori.center.UTR3.overlap.bed", header = F, sep = "\t") %>% distinct(V1, V2, V3, V4, .keep_all = T)
nrow(SNS.stochastic.ori.center.UTR3.overlap[which(SNS.stochastic.ori.center.UTR3.overlap$V10 == 0),])   # 1000 origins
# bedtools closest -a ./Figure_4/SNS.stochastic.ori.center.sorted.bed -b ./Figure_4/hsapiens.genes.coordinate.sorted.bed -d > ./Figure_4/SNS.stochastic.ori.center.genes.overlap.bed
SNS.stochastic.ori.center.genes.overlap <- read.csv("./Figure_4/SNS.stochastic.ori.center.genes.overlap.bed", header = F, sep = "\t") %>% distinct(V1, V2, V3, V4, .keep_all = T)
nrow(SNS.stochastic.ori.center.genes.overlap[which(SNS.stochastic.ori.center.genes.overlap$V11 == 0),])   # 37,198 origins

# SNS stochastic origins summary (256,600 origins)
# promoter  UTR5   intron     UTR3      exon      intergenic
# 1579      554    34006      1000      1638      217823

# Report distribution within the human genome

hsapiens.promoters <- read.csv("./Figure_4/hsapiens.promoters.1kb.sorted.bed", header = F, sep = "\t") %>% mutate(length = V3-V2)
sum(hsapiens.promoters$length) # 19,335,000 bp
hsapiens.UTR5 <- read.csv("./Figure_4/hsapiens.UTR5.coordinate.sorted.bed", header = F, sep = "\t") %>% distinct(V5, .keep_all = T) %>% mutate(length = V3-V2)
sum(hsapiens.UTR5$length) # 2,797,471 bp
hsapiens.introns <- read.csv("./Figure_4/hsapiens.introns.coordinate.sorted.bed", header = F, sep = "\t") %>% distinct(V5, .keep_all = T) %>% mutate(length = V3-V2)
sum(hsapiens.introns$length) # 147,890,852 bp
hsapiens.UTR3 <- read.csv("./Figure_4/hsapiens.UTR3.coordinate.sorted.bed", header = F, sep = "\t") %>% distinct(V5, .keep_all = T)  %>% mutate(length = V3-V2)
sum(hsapiens.UTR3$length) # 24,911,816 bp
hsapiens.exons <- read.csv("./Figure_4/hsapiens.exons.UTRs.coordinate.bed", header = F, sep = "\t") %>% distinct(V5, .keep_all = T)  %>% mutate(length = V3-V2)
sum(hsapiens.exons$length) # 15,495,580 bp

# All genomic sites summary (3,272,116,950 bp GRCh38.p13)
# promoter     UTR5        intron       UTR3        exon         intergenic   in bp
# 19335000      2797471    147890852    24911816    15495580     3061686231

# Prepare barplots

Iniseq.data <- tibble(
  category = c("promoter", "5'-UTR", "exon", "intron", "3'-UTR", "intergenic"),
  count = c(2216,3897,3264,9741,526,4261)) %>% mutate(class = "constitutive")
Iniseq.data$fraction <- Iniseq.data$count/sum(Iniseq.data$count)
SNS.core.data <- tibble(
  category = c("promoter", "5'-UTR", "exon", "intron", "3'-UTR", "intergenic"),
  count = c(564,166,659,13206,373,12294)) %>% mutate(class = "core")
SNS.core.data$fraction <- SNS.core.data$count/sum(SNS.core.data$count)
SNS.stochastic.data <- tibble(
  category = c("promoter", "5'-UTR", "exon", "intron", "3'-UTR", "intergenic"),
  count = c(1579,554,1638,34006,1000,217823)) %>% mutate(class = "stochastic")
SNS.stochastic.data$fraction <- SNS.stochastic.data$count/sum(SNS.stochastic.data$count)
All.sites.data <- tibble(
  category = c("promoter", "5'-UTR", "exon", "intron", "3'-UTR", "intergenic"),
  count = c(19335000,2797471,15495580,147890852,24911816,3061686231)) %>% mutate(class = "all")
All.sites.data$fraction <- All.sites.data$count/sum(All.sites.data$count)
Genomic.dist <- rbind(Iniseq.data, SNS.core.data, SNS.stochastic.data, All.sites.data)

pdf("./Rplot/Fig_4A.pdf", width=4, height=4, useDingbats=FALSE)
Genomic.dist %>% mutate(Class = fct_relevel(class, "all", "stochastic", "core", "constitutive")) %>%
  mutate(Category = fct_relevel(category, "promoter", "5'-UTR", "exon", "intron", "3'-UTR", "intergenic")) %>%
  ggplot(aes(x = Class, y = fraction, fill = Category)) + 
  geom_bar(position="stack", stat="identity") + ggtitle("All pval < 2.2e-16\nOR = 2.35, 8.54, 12.78") +
  scale_fill_manual(values = c("#D55E00", "#E69F00", "#ADCC54", "#2EBAED", "#F0D0CE", "#D4D2D2")) +
  theme_bw() + theme(aspect.ratio=2, plot.title = element_text(size=10), axis.title=element_text(size=14), axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1), axis.text.y = element_text(size = 12), axis.title.x=element_blank(), axis.title.y=element_text(vjust=4))
dev.off() 

table.IniSeq <- rbind(c(19644,23905),c(210430719, 3272116950))
table.SNS.core <- rbind(c(14968,27262),c(210430719, 3272116950))  # OR = 8.537402, p-value < 2.2e-16
table.SNS.stochastic <- rbind(c(38777,256600),c(210430719, 3272116950))
chisq.test(table.IniSeq)    # OR = 12.77794, p-value < 2.2e-16
chisq.test(table.SNS.core)   # OR = 8.537402, p-value < 2.2e-16
chisq.test(table.SNS.stochastic)  # OR = 2.349835, p-value < 2.2e-16

############################################################################
# Analyse distribution of origins around TSS

# Use bedtools closest
# For origins centered analysis
# sort -k1,1 -k2,2n ./Figure_4/hsapiens.genes.start.bed > ./Figure_4/hsapiens.genes.start.sorted.bed
# bedtools closest -a ./Figure_4/hsapiens.genes.start.sorted.bed -b ./Figure_4/IniSeq.ori.center.sorted.bed -D ref > ./Figure_4/hsapiens.genes.start.closest.IniSeq.center.bed
# Load results
hsapiens.gene.start.closest.IniSeq.ori <- read.csv("./Figure_4/hsapiens.genes.start.closest.IniSeq.center.bed", header = F, sep = "\t") %>% dplyr::select(V1, V2, V4, V5, V6, V8, V10, V11)
colnames(hsapiens.gene.start.closest.IniSeq.ori) <- c("chr", "gene.start", "gene.strand", "gene.id", "gene.name", "ori.pos", "ori.eff", "dist")
hsapiens.gene.start.closest.IniSeq.ori <- hsapiens.gene.start.closest.IniSeq.ori %>% mutate(dist.stranded = case_when(gene.strand == 1 ~ dist, TRUE ~ -dist))
# for 19,335 genes
# Transform data
hsapiens.gene.start.closest.IniSeq.ori.10kb <- hsapiens.gene.start.closest.IniSeq.ori %>% filter(dist.stranded >= -10000 & dist.stranded <= 10000) %>%
  mutate(dist.stranded = ((as.integer(cut(dist.stranded, breaks = 4000)))-2000)*5) # per 5bp bin
# Plot distribution of origins
hsapiens.gene.start.plot <- hsapiens.gene.start.closest.IniSeq.ori.10kb %>% ggplot(aes(x=dist.stranded)) + geom_histogram(binwidth=5, fill = "#DE1C1470") +
  theme_bw() + theme(aspect.ratio=0.5, plot.title = element_text(size=12)) + ylab("Gene count") + ggtitle("TSS, n = 19335") + xlab("Distance (bp)") +
  geom_vline(xintercept = 0, linetype="dashed", size = 0.5)

############################################################################
# Analyse distribution of origins at Exon/intron junctions

hsapiens.introns <- read.csv("./Figure_4/hsapiens.introns.coordinate.bed", header = F, sep = "\t")
colnames(hsapiens.introns) <- c("chr", "intron.start", "intron.end", "gene.strand", "gene.id", "gene.name")
# Prepare file reporting exon/intron and intron/exon junctions
hsapiens.exon.intron <- hsapiens.introns %>% mutate(start = intron.start, end = intron.start + 1) %>% dplyr::select(chr, start, end, gene.strand, gene.id, gene.name)
hsapiens.intron.exon <- hsapiens.introns %>% mutate(start = intron.end, end = intron.end + 1) %>% dplyr::select(chr, start, end, gene.strand, gene.id, gene.name)

write.table(hsapiens.exon.intron, "./Figure_4/hsapiens.exon.intron.junctions.bed", sep="\t", col.names = F, row.names = F, quote = F)
write.table(hsapiens.intron.exon, "./Figure_4/hsapiens.intron.exon.junctions.bed", sep="\t", col.names = F, row.names = F, quote = F)

# sort -k1,1 -k2,2n ./Figure_4/hsapiens.exon.intron.junctions.bed > ./Figure_4/hsapiens.exon.intron.junctions.sorted.bed
# sort -k1,1 -k2,2n ./Figure_4/hsapiens.intron.exon.junctions.bed > ./Figure_4/hsapiens.intron.exon.junctions.sorted.bed
# bedtools closest -a ./Figure_4/hsapiens.exon.intron.junctions.sorted.bed -b ./Figure_4/IniSeq.ori.center.sorted.bed -D ref > ./Figure_4/hsapiens.exon.intron.junctions.closest.IniSeq.center.bed
# bedtools closest -a ./Figure_4/hsapiens.intron.exon.junctions.sorted.bed -b ./Figure_4/IniSeq.ori.center.sorted.bed -D ref > ./Figure_4/hsapiens.intron.exon.junctions.closest.IniSeq.center.bed
# Load results
hsapiens.exon.intron.junctions.closest.IniSeq.ori <- read.csv("./Figure_4/hsapiens.exon.intron.junctions.closest.IniSeq.center.bed", header = F, sep = "\t") %>% dplyr::select(V1, V2, V4, V5, V6, V8, V10, V11)
colnames(hsapiens.exon.intron.junctions.closest.IniSeq.ori) <- c("chr", "exon.intron.pos", "gene.strand", "gene.id", "gene.name", "ori.pos", "ori.eff", "dist")
# 214,145 junctions
hsapiens.intron.exon.junctions.closest.IniSeq.ori <- read.csv("./Figure_4/hsapiens.intron.exon.junctions.closest.IniSeq.center.bed", header = F, sep = "\t") %>% dplyr::select(V1, V2, V4, V5, V6, V8, V10, V11)
colnames(hsapiens.intron.exon.junctions.closest.IniSeq.ori) <- c("chr", "intron.exon.pos", "gene.strand", "gene.id", "gene.name", "ori.pos", "ori.eff", "dist")
# 214,148 junctions
# Transform data
hsapiens.exon.intron.junctions.closest.IniSeq.ori.10kb <- hsapiens.exon.intron.junctions.closest.IniSeq.ori %>% filter(dist >= -10000 & dist <= 10000) %>%
  mutate(dist = ((as.integer(cut(dist, breaks = 4000)))-2000)*5) # per 5bp bin
hsapiens.intron.exon.junctions.closest.IniSeq.ori.10kb <- hsapiens.intron.exon.junctions.closest.IniSeq.ori %>% filter(dist >= -10000 & dist <= 10000) %>%
  mutate(dist = ((as.integer(cut(dist, breaks = 4000)))-2000)*5) # per 5bp bin
# Plot distribution of origins
hsapiens.exon.intron.plot <- hsapiens.exon.intron.junctions.closest.IniSeq.ori.10kb %>% ggplot(aes(x=dist)) + geom_histogram(binwidth=5, fill = "#009E73") +
  ylim(0,150) + theme_bw() + theme(aspect.ratio=0.5, plot.title = element_text(size=12)) + ylab("Gene count") + ggtitle("Exon/Intron junctions, n = 214,145") + xlab("Exon/Intron junctions") +
  geom_vline(xintercept = 0, linetype="dashed", size = 0.5)
hsapiens.intron.exon.plot <- hsapiens.intron.exon.junctions.closest.IniSeq.ori.10kb %>% ggplot(aes(x=dist)) + geom_histogram(binwidth=5, fill = "#009E73") +
  ylim(0,150) + theme_bw() + theme(aspect.ratio=0.5, plot.title = element_text(size=12))  + ylab("Gene count") + ggtitle("Intron/Exon junctions, n = 214,148") + xlab("Intron/Exon junctions") +
  geom_vline(xintercept = 0, linetype="dashed", size = 0.5)

# Save plots

pdf("./Rplot/Fig_4B_4C.pdf", width=15, height=4, useDingbats=FALSE)
ggarrange(hsapiens.gene.start.plot, hsapiens.exon.intron.plot, hsapiens.intron.exon.plot, ncol = 3, nrow = 1)
dev.off()

############################################################################
############################################################################
###                     DNA breaks distribution                          ###
###                         at origins                                   ###
############################################################################
############################################################################

# Raw and processed sequencing data used to contruct Fig 4D to 4I are available 
# at GEO (https://www.ncbi.nlm.nih.gov/geo/)
# under accession number GSE202799

############################################################################
# Format bed files for TSS and splice junctions

genes.TSS.with.ori <- read.csv("./Figure_4/genes.TSS.with.ori.bed", row.names = NULL) %>% dplyr::select(-X) # 11,410
genes.TSS.without.ori <- read.csv("./Figure_4/genes.TSS.without.ori.bed", row.names = NULL) %>% dplyr::select(-X) # 7,925
splice.sites.with.ori <- read.csv("./Figure_4/splice.sites.with.ori.bed", row.names = NULL) %>% dplyr::select(-X) # 100,370
splice.sites.without.ori <- read.csv("./Figure_4/splice.sites.without.ori.bed", row.names = NULL) %>% dplyr::select(-X) # 327,916

# Prepare bed files

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

# Add gene expression information

# Load gene expression data from H9 cells

H9.1 <- read.table("./Dataset/GSE202799_RAW/GSM6133301_H9_RNAseq_rep1.counts.txt.gz")
colnames(H9.1) <- c("gene.id", "H9.1")
H9.2 <- read.table("./Dataset/GSE202799_RAW/GSM6133302_H9_RNAseq_rep2.counts.txt.gz")
colnames(H9.2) <- c("gene.id", "H9.2")
H9.3 <- read.table("./Dataset/GSE202799_RAW/GSM6133303_H9_RNAseq_rep3.counts.txt.gz")
colnames(H9.3) <- c("gene.id", "H9.3")
H9.df <- H9.1 %>% full_join(H9.2, by = "gene.id") %>% full_join(H9.3, by = "gene.id")

# Compute TPM

# Recover gene length
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
att <- listAttributes(mart)
# Select only protein coding genes
hsapiens.genes <-  getBM(attributes=c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand"), mart=mart)
# Remove novel protein/transcript and unusual chromosomes
my.chr <- c(1:22, "X", "Y")
hsapiens.genes.2 <- hsapiens.genes %>% filter(chromosome_name %in% my.chr )  # 61,461 unique genes
# Compute gene length
hsapiens.genes.2 <- hsapiens.genes.2 %>% mutate(length = end_position - start_position) %>% 
  dplyr::select(gene.id = ensembl_gene_id, length)
# Compute TPM
H9.df.2 <- H9.df %>% left_join(hsapiens.genes.2, by = "gene.id") %>% 
  mutate(H9.1.RPK = H9.1/(length/1000), H9.2.RPK = H9.2/(length/1000), H9.3.RPK = H9.3/(length/1000)) %>% 
  mutate(H9.1.RPKM = (H9.1.RPK/sum(H9.1.RPK, na.rm = T))*1000000, H9.2.RPKM = (H9.2.RPK/sum(H9.2.RPK, na.rm = T))*1000000, H9.3.RPKM = (H9.3.RPK/sum(H9.3.RPK, na.rm = T))*1000000) %>% 
  mutate(H9.RPKM = (H9.1.RPKM+H9.2.RPKM+H9.3.RPKM)/3) %>% 
  dplyr::select(gene.id, H9.RPKM) %>% 
  rename(EXP = H9.RPKM) %>% drop_na()
# Set all TPM values below 0.01 to 0
H9.df.3 <- H9.df.2 %>% mutate(EXP = case_when(EXP < 0.01 ~ 0, T ~ EXP))
# Compute quantiles for EXP bins
H9.df.4 <- H9.df.3 %>% filter(EXP != 0)
quantile(H9.df.4$EXP, c(1/3,2/3))
# 33.33333% 66.66667% 
#   0.7756252 5.9283778 
# Add EXP quantile information
H9.df.5 <- H9.df.3 %>% mutate(exp = case_when(EXP == 0 ~ "non",
                                              EXP > 0 & EXP <= 0.7756252 ~ "low",
                                              EXP > 0.7756252 & EXP <= 5.9283778 ~ "medium",
                                              T ~ "high"))

############################################################################
# Define intergenic, genic and TSS origins

ori.center <- read.csv("./Figure_4/IniSeq.ori.center.bed", header = F, sep = "\t")
colnames(ori.center) <- c("chr", "start", "end", "EFF")
# Add ori id
ori.center <- ori.center %>% mutate(ori.id = paste(chr, start, sep = "_"))
# Assess ori/gene overlap
ori.gene <- read.csv("./Figure_4/IniSeq.ori.center.genes.overlap.bed", header = F, sep = "\t") %>% distinct(V1, V2, V3, V4, .keep_all = T) %>% dplyr::select(-V5, -V6, -V7) # 23905
colnames(ori.gene) <- c("chr", "start", "end", "EFF", "strand", "gene.id", "gene.name", "dist")
ori.gene <- ori.gene %>% mutate(ori.id = paste(chr, start, sep = "_")) %>% filter(dist == 0) # 17428
# Assess ori/TSS overlap
ori.promoter <- read.csv("./Figure_4/IniSeq.ori.center.promoters.overlap.bed", header = F, sep = "\t") %>% distinct(V1, V2, V3, V4, .keep_all = T)  %>% dplyr::select(-V5, -V6, -V7) # 23905
colnames(ori.promoter) <- c("chr", "start", "end", "EFF", "strand", "gene.id", "gene.name", "dist")
ori.promoter <- ori.promoter %>% mutate(ori.id = paste(chr, start, sep = "_")) %>% filter(dist == 0) # 2216
#
ori.class <- ori.center %>% mutate(class = case_when(ori.id %in% ori.gene$ori.id ~ "genic", T ~ "intergenic")) %>% mutate(class = case_when(ori.id %in% ori.promoter$ori.id ~ "TSS", T ~ class)) %>% 
  mutate(strand = "*")

############################################################################
# Add efficiency to intergenic/genic origins

ori.class <- ori.class %>% mutate(eff = case_when(EFF < 0.8418 ~ "low",
                                                  EFF >= 0.8418 & EFF < 0.9065 ~ "medium",
                                                  T ~ "high"))
ori.class.bed <- ori.class %>%  dplyr::select(chr, start, end, ori.id, EFF, class, eff)
write.table(ori.class.bed, "./Figure_4/ori.class.bed", sep="\t", col.names = F, row.names = F, quote = F)

############################################################################
# Format DSB and DNM datasets

# DSB and DNM counts in 50 nt window covering the human genomes were prepared as bigwig files

# Prepare tiled hg38 genome

seqi_38  <- SeqinfoForUCSCGenome("hg38")
hg38.50nt <- tileGenome(seqi_38, tilewidth = 50, cut.last.tile.in.chrom = TRUE)
my.chr <- c(paste("chr", seq(1,22,1), sep = ""), "chrX", "chrY", "chrM")
hg38.50nt.bed <- as.data.frame(hg38.50nt) %>% filter(seqnames %in% my.chr) %>% dplyr::select(seqnames, start, end)
write.table(hg38.50nt.bed, "./Dataset/hg38.tile.50nt.bed", sep="\t", col.names = F, row.names = F, quote = F)

# For DSB

# Combine replicates
# using zsh interactive shell on Mac with bedtools2 and UCSC bedGraphToBigWig
# cat GSM6133299_H9_INDUCEseq_rep1_breakends.bed GSM6133300_H9_INDUCEseq_rep2_breakends.bed > UD.breakends.bed
# sort -k1,1 -k2,2n UD.breakends.bed > UD.breakends.sorted.bed
# windowBed -c -w 0 -a hg38.tile.50nt.bed -b UD.breakends.sorted.bed > UD.breakends.count.50nt.bedgraph
# LC_COLLATE=C sort -k1,1 -k2,2n UD.breakends.count.50nt.bedgraph > UD.breakends.count.50nt.sorted.bedgraph
# bedGraphToBigWig UD.breakends.count.50nt.sorted.bedgraph hg38.chrom.sizes > UD.DSB.coverage.50nt.bw
# where hg38.chrom.sizes can be found at https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/

# For DNM

# Starting from DNM.SNV.hg38.bed as reported in Figure_1.R script
# sort -k1,1 -k2,2n DNM.SNV.hg38.bed > DNM.SNV.hg38.sorted.bed
# bedtools genomecov -bga -i DNM.SNV.hg38.sorted.bed -g hg38.chrom.sizes > DNM.SNV.hg38.bedgraph
# sort -k1,1 -k2,2n DNM.SNV.hg38.bedgraph > DNM.SNV.hg38.sorted.bedgraph
# bedGraphToBigWig DNM.SNV.hg38.sorted.bedgraph hg38.chrom.sizes DNM.SNV.hg38.bw

#######################################################################
# Assess correlation between DSB and DNM density at origins in H9 cells

# Prepare grange object
IniSeq.ori.eff <- read.table("./Dataset/GSM5658908_Ini-seq2.called.replication.origins.bed")
colnames(IniSeq.ori.eff) <- c("chr", "start", "end", "EFF")
IniSeq.ori.eff <- IniSeq.ori.eff %>% mutate(strand = "*", eff = case_when(EFF < 0.8418 ~ "low",
                                                                          EFF >= 0.8418 & EFF < 0.9065 ~ "medium",
                                                                          T ~ "high")) %>% mutate(ori.center = (as.integer(round(start+end)/2)), ori.id = paste(chr, ori.center, sep = "_")) %>% dplyr::select(-ori.center) %>% 
  mutate(chr = paste("chr", chr, sep = ""))
ori.class.2 <- ori.class %>% dplyr::select(ori.id, class)
IniSeq.ori.eff <- IniSeq.ori.eff %>% left_join(ori.class.2, by = "ori.id")
# Load bed file to grange
ori.gr <- makeGRangesFromDataFrame(IniSeq.ori.eff, keep.extra.columns=T)

# Open BW data
DSB.gr <- import("./Dataset/UD.DSB.coverage.50nt.bw")
DNM.gr <- import("./Dataset/DNM.SNV.hg38.bw")
TOP2B.gr <- import("./Dataset/GSE202796_H9_TOP2B_bs50.bw")
seqlevelsStyle(TOP2B.gr) <- "UCSC"
TOP2A.gr <- import("./Dataset/GSE202796_H9_TOP2A_bs50.bw")
seqlevelsStyle(TOP2A.gr) <- "UCSC"
ATAC.gr <- import("./Dataset/GSE202795_H9_ATACseq_bs50.bw")
seqlevelsStyle(ATAC.gr) <- "UCSC"

# Compute origin overlaps
ori.DSB.overlap <- as.data.frame(mergeByOverlaps(ori.gr, DSB.gr)) %>% group_by(ori.id) %>% summarise(DSB.count = sum(DSB.gr.score))
ori.DNM.overlap <- as.data.frame(mergeByOverlaps(ori.gr, DNM.gr)) %>% group_by(ori.id) %>% summarise(DNM.count = sum(DNM.gr.score))
ori.TOP2B.overlap <- as.data.frame(mergeByOverlaps(ori.gr, TOP2B.gr)) %>% group_by(ori.id) %>% summarise(TOP2B.enr = mean(TOP2B.gr.score, na.rm = T))
ori.TOP2A.overlap <- as.data.frame(mergeByOverlaps(ori.gr, TOP2A.gr)) %>% group_by(ori.id) %>% summarise(TOP2A.enr = mean(TOP2A.gr.score, na.rm = T))
ori.ATAC.overlap <- as.data.frame(mergeByOverlaps(ori.gr, ATAC.gr)) %>% group_by(ori.id) %>% summarise(ATAC.enr = mean(ATAC.gr.score, na.rm = T))
# join results
ori.DSB.DNM.df <- IniSeq.ori.eff %>% left_join(ori.DSB.overlap, by = "ori.id") %>% left_join(ori.DNM.overlap, by = "ori.id") %>% left_join(ori.TOP2A.overlap, by = "ori.id") %>% left_join(ori.TOP2B.overlap, by = "ori.id")  %>% left_join(ori.ATAC.overlap, by = "ori.id")

ks.test(ori.DSB.DNM.df[which(ori.DSB.DNM.df$eff == "low"),]$DSB.count,  ori.DSB.DNM.df[which(ori.DSB.DNM.df$eff == "high"),]$DSB.count) # p-value < 2.2e-16
DSB.eff.plot <- ori.DSB.DNM.df %>% mutate(eff2 = fct_relevel(eff, "low", "medium", "high")) %>%
  ggplot(aes(eff2, DSB.count, fill = eff2)) +
  geom_boxplot() + ggtitle("DSB count vs Ori efficiency\nP < 2.2e-16") + ylab("DSB count") +
  scale_fill_manual(values = c("#CC79A7", "#0072B2","#FC4E07")) +
  scale_y_continuous(trans = 'log2') +
  theme_bw() + theme(aspect.ratio=2, plot.title = element_text(size=14), axis.title=element_text(size=14), axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1), axis.text.y = element_text(size = 12), axis.title.x=element_blank(), axis.title.y=element_text(vjust=4), legend.position = "none")

# Stats
DSB.low <- ori.DSB.DNM.df %>% filter(eff == "low") %>% dplyr::select(eff, DSB.count)
DSB.medium <- ori.DSB.DNM.df %>% filter(eff == "medium") %>% dplyr::select(eff, DSB.count)
DSB.high <- ori.DSB.DNM.df %>% filter(eff == "high") %>% dplyr::select(eff, DSB.count)
test.df <- rbind(DSB.low, DSB.medium, DSB.high)
test <- chisq.test(t(table(test.df$eff, test.df$DSB.count)))
test  # p-value < 2.2e-16

DSB.eff.plot.inter <- ori.DSB.DNM.df %>% filter(class == "intergenic") %>% mutate(eff2 = fct_relevel(eff, "low", "medium", "high")) %>%
  ggplot(aes(eff2, DSB.count, fill = eff2)) +
  geom_boxplot() + ggtitle("DSB count vs Ori inter efficiency\nP < 2.2e-16") + ylab("DSB count") +
  scale_fill_manual(values = c("#CC79A7", "#0072B2","#FC4E07")) +
  scale_y_continuous(trans = 'log2') +
  theme_bw() + theme(aspect.ratio=2, plot.title = element_text(size=14), axis.title=element_text(size=14), axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1), axis.text.y = element_text(size = 12), axis.title.x=element_blank(), axis.title.y=element_text(vjust=4), legend.position = "none")

# Stats
DSB.low <- ori.DSB.DNM.df %>% filter(class == "intergenic") %>% filter(eff == "low") %>% dplyr::select(eff, DSB.count)
DSB.medium <- ori.DSB.DNM.df %>% filter(class == "intergenic") %>% filter(eff == "medium") %>% dplyr::select(eff, DSB.count)
DSB.high <- ori.DSB.DNM.df %>% filter(class == "intergenic") %>% filter(eff == "high") %>% dplyr::select(eff, DSB.count)
test.df <- rbind(DSB.low, DSB.medium, DSB.high)
test <- chisq.test(t(table(test.df$eff, test.df$DSB.count)))
test  # p-value < 2.2e-16

ks.test(ori.DSB.DNM.df[which(ori.DSB.DNM.df$eff == "low"),]$ATAC.enr,  ori.DSB.DNM.df[which(ori.DSB.DNM.df$eff == "high"),]$ATAC.enr) # p-value < 2.2e-16
ATAC.eff.plot <- ori.DSB.DNM.df %>% mutate(eff2 = fct_relevel(eff, "low", "medium", "high")) %>%
  ggplot(aes(eff2, ATAC.enr, fill = eff2)) +
  geom_boxplot() + ggtitle("ATAC signal vs Ori efficiency\nP < 2.2e-16") + ylab("ATAC signal") +
  scale_fill_manual(values = c("#CC79A7", "#0072B2","#FC4E07")) +
  scale_y_continuous(trans = 'log10') +
  theme_bw() + theme(aspect.ratio=2, plot.title = element_text(size=14), axis.title=element_text(size=14), axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1), axis.text.y = element_text(size = 12), axis.title.x=element_blank(), axis.title.y=element_text(vjust=4), legend.position = "none")

ATAC.eff.plot.inter <- ori.DSB.DNM.df %>% filter(class == "intergenic") %>% mutate(eff2 = fct_relevel(eff, "low", "medium", "high")) %>%
  ggplot(aes(eff2, ATAC.enr, fill = eff2)) +
  geom_boxplot() + ggtitle("ATAC signal vs Ori inter efficiency\nP < 2.2e-16") + ylab("ATAC signal") +
  scale_fill_manual(values = c("#CC79A7", "#0072B2","#FC4E07")) +
  scale_y_continuous(trans = 'log10') +
  theme_bw() + theme(aspect.ratio=2, plot.title = element_text(size=14), axis.title=element_text(size=14), axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1), axis.text.y = element_text(size = 12), axis.title.x=element_blank(), axis.title.y=element_text(vjust=4), legend.position = "none")

# Cut DSB count in quantiles
# Add a pseudo count to DNM count

ori.DSB.DNM.df.2 <- ori.DSB.DNM.df %>% mutate(DSB.bin = paste("Q", as.character(ntile(DSB.count, 10)), sep = "")) %>% mutate(DNM.pseudo.count = DNM.count + 1)

colfunc <- colorRampPalette(c("white", "#F21A00"))
colfunc(10)

ks.test(ori.DSB.DNM.df.2[which(ori.DSB.DNM.df.2$DSB.bin == "Q1"),]$ATAC.enr,  ori.DSB.DNM.df.2[which(ori.DSB.DNM.df.2$DSB.bin == "Q10"),]$ATAC.enr) # p-value < 2.2e-16
ATAC.DSB.plot <- ori.DSB.DNM.df.2 %>% mutate(bin = fct_relevel(DSB.bin, "Q1", "Q2", "Q3", "Q4", "Q5", "Q6", "Q7", "Q8", "Q9", "Q10")) %>%
  ggplot(aes(bin, ATAC.enr, fill = bin)) +
  geom_boxplot() +
  ggtitle("ATAC vs DSB count (Q1-10)\nP < 2.2e-16") + ylab("ATAC signal") +
  scale_fill_manual(values = colfunc(10)) +
  scale_y_continuous(trans = 'log10') +
  theme_bw() + theme(aspect.ratio=2, plot.title = element_text(size=14), axis.title=element_text(size=14), axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1), axis.text.y = element_text(size = 12), axis.title.x=element_blank(), axis.title.y=element_text(vjust=4), legend.position = "none")

ks.test(ori.DSB.DNM.df.2[which(ori.DSB.DNM.df.2$DSB.bin == "Q1"),]$DNM.pseudo.count,  ori.DSB.DNM.df.2[which(ori.DSB.DNM.df.2$DSB.bin == "Q10"),]$DNM.pseudo.count) # p-value < 2.2e-16
DNM.DSB.plot <- ori.DSB.DNM.df.2 %>% mutate(bin = fct_relevel(DSB.bin, "Q1", "Q2", "Q3", "Q4", "Q5", "Q6", "Q7", "Q8", "Q9", "Q10")) %>%
  ggplot(aes(bin, DNM.pseudo.count, fill = bin)) +
  geom_boxplot() +
  ggtitle("DNM vs DSB count (Q1-10)\nP < 2.2e-16") + ylab("DNM count") +
  scale_fill_manual(values = colfunc(10)) +
  scale_y_continuous(trans = 'log2') +
  theme_bw() + theme(aspect.ratio=2, plot.title = element_text(size=14), axis.title=element_text(size=14), axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1), axis.text.y = element_text(size = 12), axis.title.x=element_blank(), axis.title.y=element_text(vjust=4), legend.position = "none")

# Stats
test <- chisq.test(t(table(ori.DSB.DNM.df.2$DSB.bin, ori.DSB.DNM.df.2$DNM.pseudo.count)))
test  # p-value < 2.2e-16

ks.test(ori.DSB.DNM.df.2[which(ori.DSB.DNM.df.2$DSB.bin == "Q1"),]$TOP2B.enr,  ori.DSB.DNM.df.2[which(ori.DSB.DNM.df.2$DSB.bin == "Q10"),]$TOP2B.enr) # p-value < 2.2e-16
TOP2B.DSB.plot <- ori.DSB.DNM.df.2 %>% mutate(bin = fct_relevel(DSB.bin, "Q1", "Q2", "Q3", "Q4", "Q5", "Q6", "Q7", "Q8", "Q9", "Q10")) %>%
  ggplot(aes(bin, TOP2B.enr, fill = bin)) +
  geom_boxplot() +
  ggtitle("TOP2B vs DSB count (Q1-10)\nP < 2.2e-16") + ylab("TOP2B signal") +
  scale_fill_manual(values = colfunc(10)) +
  scale_y_continuous(trans = 'log10') +
  theme_bw() + theme(aspect.ratio=2, plot.title = element_text(size=14), axis.title=element_text(size=14), axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1), axis.text.y = element_text(size = 12), axis.title.x=element_blank(), axis.title.y=element_text(vjust=4), legend.position = "none")

# Stats
test <- chisq.test(t(table(ori.DSB.DNM.df.2$DSB.bin, ori.DSB.DNM.df.2$TOP2B.enr)))
test  # p-value = 0.0264

ks.test(ori.DSB.DNM.df.2[which(ori.DSB.DNM.df.2$DSB.bin == "Q1"),]$TOP2A.enr,  ori.DSB.DNM.df.2[which(ori.DSB.DNM.df.2$DSB.bin == "Q10"),]$TOP2A.enr) # p-value < 2.2e-16
TOP2A.DSB.plot <- ori.DSB.DNM.df.2 %>% mutate(bin = fct_relevel(DSB.bin, "Q1", "Q2", "Q3", "Q4", "Q5", "Q6", "Q7", "Q8", "Q9", "Q10")) %>%
  ggplot(aes(bin, TOP2A.enr, fill = bin)) +
  geom_boxplot() +
  ggtitle("TOP2A vs DSB count (Q1-10)\nP < 2.2e-16") + ylab("TOP2A signal") +
  scale_fill_manual(values = colfunc(10)) +
  scale_y_continuous(trans = 'log10') +
  theme_bw() + theme(aspect.ratio=2, plot.title = element_text(size=14), axis.title=element_text(size=14), axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1), axis.text.y = element_text(size = 12), axis.title.x=element_blank(), axis.title.y=element_text(vjust=4), legend.position = "none")

# Stats
test <- chisq.test(t(table(ori.DSB.DNM.df.2$DSB.bin, ori.DSB.DNM.df.2$TOP2A.enr)))
test  # p-value = 0.01193

TOP2B.DSB.plot.inter <- ori.DSB.DNM.df.2 %>% filter(class == "intergenic") %>% mutate(bin = fct_relevel(DSB.bin, "Q1", "Q2", "Q3", "Q4", "Q5", "Q6", "Q7", "Q8", "Q9", "Q10")) %>%
  ggplot(aes(bin, TOP2B.enr, fill = bin)) +
  geom_boxplot() +
  ggtitle("TOP2B vs DSB count (Q1-10)\nP < 2.2e-16, intergenic") + ylab("TOP2B signal") +
  scale_fill_manual(values = colfunc(10)) +
  scale_y_continuous(trans = 'log10') +
  theme_bw() + theme(aspect.ratio=2, plot.title = element_text(size=14), axis.title=element_text(size=14), axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1), axis.text.y = element_text(size = 12), axis.title.x=element_blank(), axis.title.y=element_text(vjust=4), legend.position = "none")

# Stats
ori.DSB.DNM.df.inter <- ori.DSB.DNM.df.2 %>% filter(class == "intergenic")
test <- chisq.test(t(table(ori.DSB.DNM.df.inter$DSB.bin, ori.DSB.DNM.df.inter$TOP2B.enr)))
test  # p-value = 0.008124

TOP2A.DSB.plot.inter <- ori.DSB.DNM.df.2 %>% filter(class == "intergenic") %>% mutate(bin = fct_relevel(DSB.bin, "Q1", "Q2", "Q3", "Q4", "Q5", "Q6", "Q7", "Q8", "Q9", "Q10")) %>%
  ggplot(aes(bin, TOP2A.enr, fill = bin)) +
  geom_boxplot() +
  ggtitle("TOP2A vs DSB count (Q1-10)\nP < 2.2e-16, intergenic") + ylab("TOP2A signal") +
  scale_fill_manual(values = colfunc(10)) +
  scale_y_continuous(trans = 'log10') +
  theme_bw() + theme(aspect.ratio=2, plot.title = element_text(size=14), axis.title=element_text(size=14), axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1), axis.text.y = element_text(size = 12), axis.title.x=element_blank(), axis.title.y=element_text(vjust=4), legend.position = "none")

# Stats
test <- chisq.test(t(table(ori.DSB.DNM.df.inter$DSB.bin, ori.DSB.DNM.df.inter$TOP2A.enr)))
test  # p-value = 0.007462

pdf("./Rplot/Fig_4F_4G_4H_S4B_S4D_S4E_S4F_S4G.pdf", width=18, height=8, useDingbats=FALSE)
ggarrange(DSB.eff.plot, DSB.eff.plot.inter, ATAC.eff.plot, ATAC.eff.plot.inter, ATAC.DSB.plot, DNM.DSB.plot, TOP2B.DSB.plot, TOP2B.DSB.plot.inter, TOP2A.DSB.plot, TOP2A.DSB.plot.inter, nrow = 2, ncol = 5)
dev.off()

############################################################################
############################################################################
###                 Assess the contribution of origin                    ###
###                      usage and efficiency                            ###
###                      on DSB accumulation                             ###
###             at origins found at TSS and splice sites                 ###
############################################################################
############################################################################

###################################################################
# Prepare gr object reporting the position of all TSS and splice sites

gene.starts.bed <- read.table("./Figure_4/hsapiens.genes.start.bed", sep="\t")
colnames(gene.starts.bed) <- c("chr", "start", "end", "strand", "gene.id", "gene.name")
gene.starts.2.5.kb.bed <- gene.starts.bed %>% mutate(start = as.integer(start-2500), end = as.integer(end+2500)) %>%
  mutate(strand = case_when(strand == 1 ~ "+", T ~ "-"))

exon.intron.bed <- read.table("./Figure_4/hsapiens.exon.intron.junctions.bed", sep="\t") 
intron.exon.bed <- read.table("./Figure_4/hsapiens.intron.exon.junctions.bed", sep="\t") 
splice.sites.bed <- rbind(exon.intron.bed, intron.exon.bed)
colnames(splice.sites.bed) <- c("chr", "start", "end", "strand", "gene.id", "gene.name")
splice.sites.2.5.kb.bed <- splice.sites.bed %>% mutate(start = as.integer(start-2500), end = as.integer(end+2500)) %>%
  mutate(strand = case_when(strand == 1 ~ "+", T ~ "-"))

all.TSS.gr <- makeGRangesFromDataFrame(gene.starts.2.5.kb.bed, keep.extra.columns = T)
seqlevelsStyle(all.TSS.gr) <- "UCSC"
all.splice.sites.gr <- makeGRangesFromDataFrame(splice.sites.2.5.kb.bed, keep.extra.columns = T)
seqlevelsStyle(all.splice.sites.gr) <- "UCSC"

# Load origins coordinates

# Constitutive origins
IniSeq.ori.eff <- read.table("./Dataset/GSM5658908_Ini-seq2.called.replication.origins.bed")
colnames(IniSeq.ori.eff) <- c("chr", "start", "end", "EFF")
IniSeq.ori.eff <- IniSeq.ori.eff %>% mutate(strand = "*", eff = case_when(EFF < 0.8418 ~ "low",
                                                                          EFF >= 0.8418 & EFF < 0.9065 ~ "medium",
                                                                          T ~ "high")) %>% mutate(ori.center = (as.integer(round(start+end)/2)), ori.id = paste(chr, ori.center, sep = "_")) %>% dplyr::select(-ori.center) %>% 
  mutate(chr = paste("chr", chr, sep = ""))
# Load bed file to grange
IniSeq.ori.gr <- makeGRangesFromDataFrame(IniSeq.ori.eff, keep.extra.columns=T)

# Core origins 
Core.ori.bed <- read.table("./Figure_1/Core_origins_unique.bed", sep="\t")
colnames(Core.ori.bed) <- c("chr", "start", "end", "ori.id")
Core.ori.bed.2 <- Core.ori.bed %>% mutate(EFF = NA, eff = "core", strand = "*") %>% dplyr::select(chr, start, end, strand, EFF, eff, ori.id)
Core.ori.gr <- makeGRangesFromDataFrame(Core.ori.bed.2, keep.extra.columns=T)

# Stochastic origins
Stochastic.ori.bed <- read.table("./Dataset/GSE128477_Stochastic_origins_hg38.bed")
colnames(Stochastic.ori.bed) <- c("chr", "start", "end", "ori.id")
Stochastic.ori.bed.2 <- Stochastic.ori.bed %>% mutate(EFF = NA, eff = "stochastic", strand = "*") %>% dplyr::select(chr, start, end, strand, EFF, eff, ori.id)
Stochastic.ori.gr <- makeGRangesFromDataFrame(Stochastic.ori.bed.2, keep.extra.columns=T)

###################################################################
# Identify origins associated with TSS and splice sites

IniSeq.TSS.ori.gr <- unique(IniSeq.ori.gr[queryHits(findOverlaps(IniSeq.ori.gr, all.TSS.gr))])
IniSeq.splice.sites.ori.gr <- unique(IniSeq.ori.gr[queryHits(findOverlaps(IniSeq.ori.gr, all.splice.sites.gr))])

Core.TSS.ori.gr <- unique(Core.ori.gr[queryHits(findOverlaps(Core.ori.gr, all.TSS.gr))])
Core.splice.sites.ori.gr <- unique(Core.ori.gr[queryHits(findOverlaps(Core.ori.gr, all.splice.sites.gr))])

Stochastic.TSS.ori.gr <- unique(Stochastic.ori.gr[queryHits(findOverlaps(Stochastic.ori.gr, all.TSS.gr))])
Stochastic.splice.sites.ori.gr <- unique(Stochastic.ori.gr[queryHits(findOverlaps(Stochastic.ori.gr, all.splice.sites.gr))])

###################################################################
# Compute number of breaks associated with each origin

IniSeq.TSS.ori.DSB.overlap <- as.data.frame(mergeByOverlaps(IniSeq.TSS.ori.gr, DSB.gr)) %>% group_by(ori.id) %>% summarise(DSB.count = sum(DSB.gr.score))
Core.TSS.ori.DSB.overlap <- as.data.frame(mergeByOverlaps(Core.TSS.ori.gr, DSB.gr)) %>% group_by(ori.id) %>% summarise(DSB.count = sum(DSB.gr.score))
Stochastic.TSS.ori.DSB.overlap <- as.data.frame(mergeByOverlaps(Stochastic.TSS.ori.gr, DSB.gr)) %>% group_by(ori.id) %>% summarise(DSB.count = sum(DSB.gr.score))

IniSeq.splice.sites.ori.DSB.overlap <- as.data.frame(mergeByOverlaps(IniSeq.splice.sites.ori.gr, DSB.gr)) %>% group_by(ori.id) %>% summarise(DSB.count = sum(DSB.gr.score))
Core.splice.sites.ori.DSB.overlap <- as.data.frame(mergeByOverlaps(Core.splice.sites.ori.gr, DSB.gr)) %>% group_by(ori.id) %>% summarise(DSB.count = sum(DSB.gr.score))
Stochastic.splice.sites.ori.DSB.overlap <- as.data.frame(mergeByOverlaps(Stochastic.splice.sites.ori.gr, DSB.gr)) %>% group_by(ori.id) %>% summarise(DSB.count = sum(DSB.gr.score))

###################################################################
# Find associated gene

hsapiens.genes.coordinate <- read.table("./Figure_4/hsapiens.genes.coordinate.bed", sep="\t")
colnames(hsapiens.genes.coordinate) <- c("chr", "start", "end", "strand", "gene.id", "gene.name")
hsapiens.genes.coordinate <- hsapiens.genes.coordinate %>% mutate(strand = case_when(strand == 1 ~ "+", T ~ "-"))
hsapiens.genes.coordinate.gr <- makeGRangesFromDataFrame(hsapiens.genes.coordinate, keep.extra.columns = T)
seqlevelsStyle(hsapiens.genes.coordinate.gr) <- "UCSC"

IniSeq.ori.nearest.gene <- hsapiens.genes.coordinate.gr[nearest(IniSeq.ori.gr, hsapiens.genes.coordinate.gr)]$gene.id
IniSeq.ori.gr$gene.id <- IniSeq.ori.nearest.gene

Core.ori.nearest.gene <- hsapiens.genes.coordinate.gr[nearest(Core.ori.gr, hsapiens.genes.coordinate.gr)]$gene.id
Core.ori.gr$gene.id <- Core.ori.nearest.gene

Stochastic.ori.nearest.gene <- hsapiens.genes.coordinate.gr[nearest(Stochastic.ori.gr, hsapiens.genes.coordinate.gr)]$gene.id
Stochastic.ori.gr$gene.id <- Stochastic.ori.nearest.gene

###################################################################
# Combine all information

IniSeq.TSS.ori.df <- as.data.frame(IniSeq.ori.gr) %>% filter(ori.id %in% IniSeq.TSS.ori.gr$ori.id) %>% left_join(IniSeq.TSS.ori.DSB.overlap, by = "ori.id") %>% left_join(H9.df.5, by = "gene.id")
Core.TSS.ori.df <- as.data.frame(Core.ori.gr) %>% filter(ori.id %in% Core.TSS.ori.gr$ori.id) %>% left_join(Core.TSS.ori.DSB.overlap, by = "ori.id") %>% left_join(H9.df.5, by = "gene.id")
Stochastic.TSS.ori.df <- as.data.frame(Stochastic.ori.gr) %>% filter(ori.id %in% Stochastic.TSS.ori.gr$ori.id) %>% left_join(Stochastic.TSS.ori.DSB.overlap, by = "ori.id") %>% left_join(H9.df.5, by = "gene.id")
TSS.ori.df <- rbind(IniSeq.TSS.ori.df, Core.TSS.ori.df, Stochastic.TSS.ori.df) %>% filter(exp != "NA")

IniSeq.splice.sites.ori.df <- as.data.frame(IniSeq.ori.gr) %>% filter(ori.id %in% IniSeq.splice.sites.ori.gr$ori.id) %>% left_join(IniSeq.splice.sites.ori.DSB.overlap, by = "ori.id") %>% left_join(H9.df.5, by = "gene.id")
Core.splice.sites.ori.df <- as.data.frame(Core.ori.gr) %>% filter(ori.id %in% Core.splice.sites.ori.gr$ori.id) %>% left_join(Core.splice.sites.ori.DSB.overlap, by = "ori.id") %>% left_join(H9.df.5, by = "gene.id")
Stochastic.splice.sites.ori.df <- as.data.frame(Stochastic.ori.gr) %>% filter(ori.id %in% Stochastic.splice.sites.ori.gr$ori.id) %>% left_join(Stochastic.splice.sites.ori.DSB.overlap, by = "ori.id") %>% left_join(H9.df.5, by = "gene.id")
splice.sites.ori.df <- rbind(IniSeq.splice.sites.ori.df, Core.splice.sites.ori.df, Stochastic.splice.sites.ori.df) %>% filter(exp != "NA")

####################################################################
# Plot

DSB.TSS.exp.plot.2 <- TSS.ori.df %>% mutate(exp2 = fct_relevel(exp, "non", "low", "medium", "high"), eff2 = fct_relevel(eff, "stochastic", "core", "low", "medium", "high")) %>%
  ggplot(aes(exp2, DSB.count, fill = eff2)) +
  geom_boxplot() + ggtitle("TSS\nP = 0.002159, 2.015e-06, 4.33e-15, 2.2e-16") + ylab("DSB count") +
  scale_fill_manual(values = c("#D4D2D2", "#828388", "#CC79A7", "#0072B2","#FC4E07")) +
  scale_y_continuous(trans = 'log2', limits = c(1,2500)) + labs(fill = "Origin efficiency") +
  theme_bw() + theme(aspect.ratio=1, plot.title = element_text(size=14), axis.title=element_text(size=14), axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1), axis.text.y = element_text(size = 12), axis.title.x=element_blank(), axis.title.y=element_text(vjust=4))

DSB.junctions.exp.plot.2 <- splice.sites.ori.df %>% mutate(exp2 = fct_relevel(exp, "non", "low", "medium", "high"), eff2 = fct_relevel(eff, "stochastic", "core", "low", "medium", "high")) %>%
  ggplot(aes(exp2, DSB.count, fill = eff2)) +
  geom_boxplot() + ggtitle("Junctions\nP = 2.566e-07, 2.2e-16, 2.2e-16, 2.2e-16") + ylab("DSB count") +
  scale_fill_manual(values = c("#D4D2D2", "#828388", "#CC79A7", "#0072B2","#FC4E07")) +
  scale_y_continuous(trans = 'log2', limits = c(1,4000)) + labs(fill = "Origin efficiency") +
  theme_bw() + theme(aspect.ratio=1, plot.title = element_text(size=14), axis.title=element_text(size=14), axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1), axis.text.y = element_text(size = 12), axis.title.x=element_blank(), axis.title.y=element_text(vjust=4))

pdf("./Rplot/Fig_4I.pdf", width=12, height=4, useDingbats=FALSE)
ggarrange(DSB.TSS.exp.plot.2, DSB.junctions.exp.plot.2, ncol = 2)
dev.off()

# Stats reported on Fig 4I are computed below

# Compare constitutive origins of different efficiencies
ks.test(TSS.ori.df[which(TSS.ori.df$exp == "non" & TSS.ori.df$eff == "low"),]$DSB.count,
        TSS.ori.df[which(TSS.ori.df$exp == "non" & TSS.ori.df$eff == "high"),]$DSB.count) # p-value = 0.002159
ks.test(TSS.ori.df[which(TSS.ori.df$exp == "low" & TSS.ori.df$eff == "low"),]$DSB.count,
        TSS.ori.df[which(TSS.ori.df$exp == "low" & TSS.ori.df$eff == "high"),]$DSB.count) # p-value = 2.015e-06
ks.test(TSS.ori.df[which(TSS.ori.df$exp == "medium" & TSS.ori.df$eff == "low"),]$DSB.count,
        TSS.ori.df[which(TSS.ori.df$exp == "medium" & TSS.ori.df$eff == "high"),]$DSB.count) # p-value = 4.33e-15
ks.test(TSS.ori.df[which(TSS.ori.df$exp == "high" & TSS.ori.df$eff == "low"),]$DSB.count,
        TSS.ori.df[which(TSS.ori.df$exp == "high" & TSS.ori.df$eff == "high"),]$DSB.count) # p-value < 2.2e-16
ks.test(splice.sites.ori.df[which(splice.sites.ori.df$exp == "non" & splice.sites.ori.df$eff == "low"),]$DSB.count,
        splice.sites.ori.df[which(splice.sites.ori.df$exp == "non" & splice.sites.ori.df$eff == "high"),]$DSB.count) # p-value = 2.566e-07
ks.test(splice.sites.ori.df[which(splice.sites.ori.df$exp == "low" & splice.sites.ori.df$eff == "low"),]$DSB.count,
        splice.sites.ori.df[which(splice.sites.ori.df$exp == "low" & splice.sites.ori.df$eff == "high"),]$DSB.count) # p-value < 2.2e-16
ks.test(splice.sites.ori.df[which(splice.sites.ori.df$exp == "medium" & splice.sites.ori.df$eff == "low"),]$DSB.count,
        splice.sites.ori.df[which(splice.sites.ori.df$exp == "medium" & splice.sites.ori.df$eff == "high"),]$DSB.count) # p-value < 2.2e-16
ks.test(splice.sites.ori.df[which(splice.sites.ori.df$exp == "high" & splice.sites.ori.df$eff == "low"),]$DSB.count,
        splice.sites.ori.df[which(splice.sites.ori.df$exp == "high" & splice.sites.ori.df$eff == "high"),]$DSB.count) # p-value < 2.2e-16

# Compare constitutive origins to other origins
ks.test(TSS.ori.df[which(TSS.ori.df$exp == "non" & (TSS.ori.df$eff == "low" | TSS.ori.df$eff == "medium" | TSS.ori.df$eff == "high")),]$DSB.count,
        TSS.ori.df[which(TSS.ori.df$exp == "non" & (TSS.ori.df$eff == "stochastic" | TSS.ori.df$eff == "core")),]$DSB.count) # p-value < 2.2e-16
ks.test(TSS.ori.df[which(TSS.ori.df$exp == "low" & (TSS.ori.df$eff == "low" | TSS.ori.df$eff == "medium" | TSS.ori.df$eff == "high")),]$DSB.count,
        TSS.ori.df[which(TSS.ori.df$exp == "low" & (TSS.ori.df$eff == "stochastic" | TSS.ori.df$eff == "core")),]$DSB.count) # p-value < 2.2e-16
ks.test(TSS.ori.df[which(TSS.ori.df$exp == "medium" & (TSS.ori.df$eff == "low" | TSS.ori.df$eff == "medium" | TSS.ori.df$eff == "high")),]$DSB.count,
        TSS.ori.df[which(TSS.ori.df$exp == "medium" & (TSS.ori.df$eff == "stochastic" | TSS.ori.df$eff == "core")),]$DSB.count) # p-value < 2.2e-16
ks.test(TSS.ori.df[which(TSS.ori.df$exp == "high" & (TSS.ori.df$eff == "low" | TSS.ori.df$eff == "medium" | TSS.ori.df$eff == "high")),]$DSB.count,
        TSS.ori.df[which(TSS.ori.df$exp == "high" & (TSS.ori.df$eff == "stochastic" | TSS.ori.df$eff == "core")),]$DSB.count) # p-value < 2.2e-16
ks.test(splice.sites.ori.df[which(splice.sites.ori.df$exp == "non" & (splice.sites.ori.df$eff == "low" | splice.sites.ori.df$eff == "medium" | splice.sites.ori.df$eff == "high")),]$DSB.count,
        splice.sites.ori.df[which(splice.sites.ori.df$exp == "non" & (splice.sites.ori.df$eff == "stochastic" | splice.sites.ori.df$eff == "core")),]$DSB.count) # p-value < 2.2e-16
ks.test(splice.sites.ori.df[which(splice.sites.ori.df$exp == "low" & (splice.sites.ori.df$eff == "low" | splice.sites.ori.df$eff == "medium" | splice.sites.ori.df$eff == "high")),]$DSB.count,
        splice.sites.ori.df[which(splice.sites.ori.df$exp == "low" & (splice.sites.ori.df$eff == "stochastic" | splice.sites.ori.df$eff == "core")),]$DSB.count) # p-value < 2.2e-16
ks.test(splice.sites.ori.df[which(splice.sites.ori.df$exp == "medium" & (splice.sites.ori.df$eff == "low" | splice.sites.ori.df$eff == "medium" | splice.sites.ori.df$eff == "high")),]$DSB.count,
        splice.sites.ori.df[which(splice.sites.ori.df$exp == "medium" & (splice.sites.ori.df$eff == "stochastic" | splice.sites.ori.df$eff == "core")),]$DSB.count) # p-value < 2.2e-16
ks.test(splice.sites.ori.df[which(splice.sites.ori.df$exp == "high" & (splice.sites.ori.df$eff == "low" | splice.sites.ori.df$eff == "medium" | splice.sites.ori.df$eff == "high")),]$DSB.count,
        splice.sites.ori.df[which(splice.sites.ori.df$exp == "high" & (splice.sites.ori.df$eff == "stochastic" | splice.sites.ori.df$eff == "core")),]$DSB.count) # p-value < 2.2e-16

# Compare core and stochastic origins
ks.test(TSS.ori.df[which(TSS.ori.df$exp == "non" & TSS.ori.df$eff == "core"),]$DSB.count,
        TSS.ori.df[which(TSS.ori.df$exp == "non" & TSS.ori.df$eff == "stochastic"),]$DSB.count) # p-value = 1.302e-10
ks.test(TSS.ori.df[which(TSS.ori.df$exp == "low" & TSS.ori.df$eff == "core"),]$DSB.count,
        TSS.ori.df[which(TSS.ori.df$exp == "low" & TSS.ori.df$eff == "stochastic"),]$DSB.count) # p-value < 2.2e-16
ks.test(TSS.ori.df[which(TSS.ori.df$exp == "medium" & TSS.ori.df$eff == "core"),]$DSB.count,
        TSS.ori.df[which(TSS.ori.df$exp == "medium" & TSS.ori.df$eff == "stochastic"),]$DSB.count) # p-value < 2.2e-16
ks.test(TSS.ori.df[which(TSS.ori.df$exp == "high" & TSS.ori.df$eff == "core"),]$DSB.count,
        TSS.ori.df[which(TSS.ori.df$exp == "high" & TSS.ori.df$eff == "stochastic"),]$DSB.count) # p-value < 2.2e-16
ks.test(splice.sites.ori.df[which(splice.sites.ori.df$exp == "non" & splice.sites.ori.df$eff == "core"),]$DSB.count,
        splice.sites.ori.df[which(splice.sites.ori.df$exp == "non" & splice.sites.ori.df$eff == "stochastic"),]$DSB.count) # p-value < 2.2e-16
ks.test(splice.sites.ori.df[which(splice.sites.ori.df$exp == "low" & splice.sites.ori.df$eff == "core"),]$DSB.count,
        splice.sites.ori.df[which(splice.sites.ori.df$exp == "low" & splice.sites.ori.df$eff == "stochastic"),]$DSB.count) # p-value < 2.2e-16
ks.test(splice.sites.ori.df[which(splice.sites.ori.df$exp == "medium" & splice.sites.ori.df$eff == "core"),]$DSB.count,
        splice.sites.ori.df[which(splice.sites.ori.df$exp == "medium" & splice.sites.ori.df$eff == "stochastic"),]$DSB.count) # p-value < 2.2e-16
ks.test(splice.sites.ori.df[which(splice.sites.ori.df$exp == "high" & splice.sites.ori.df$eff == "core"),]$DSB.count,
        splice.sites.ori.df[which(splice.sites.ori.df$exp == "high" & splice.sites.ori.df$eff == "stochastic"),]$DSB.count) # p-value < 2.2e-16
ks.test(TSS.ori.df[which(TSS.ori.df$exp == "medium" & TSS.ori.df$eff == "high"),]$DSB.count,
        TSS.ori.df[which(TSS.ori.df$exp == "high" & TSS.ori.df$eff == "high"),]$DSB.count) # p-value = 0.001606


























