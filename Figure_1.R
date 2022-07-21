
# P. Murat, MRC Laboratory of Molecular Biology, July 2022
# Code for generating Figure 1 with associated supporting Figures
# in Murat et al. DNA replication initiation shapes the mutational landscape and expression of the human genome

library(dplyr)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
library(ggplot2)
library(ggpubr)
library(forcats)

setwd("/Volumes/pmurat/OriMut")

############################################################################
############################################################################
###                Replication origins classification                    ###
############################################################################
############################################################################

############################################################################
# Load position of mapped origins (SNS and IniSeq)

# SNS-seq data from
# Akerman et al. Nature Communication, 2020, 11, 1, 11-15
# Downloaded from GSE128477
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE128477
# GSE128477_Core_origins_hg38.bed       64,148 regions, termed “core origins”, as replication initiation hotspots, irrespective of the cell type
# GSE128477_Stochastic_origins_hg38.bed    The remaining 80% of IS (Q3–Q10, 256,600 regions), hereby termed “stochastic origins”, had low mean activity across 19 samples and only hosted ~15–30% of total SNS-seq signal in each cell type
# Genomic location are in hg38

SNS.core <- read.table("./Dataset/GSE128477_Core_origins_hg38.bed") # 64,148 origins
colnames(SNS.core) <- c("chr", "start", "end", "ID")

SNS.stochastic <- read.table("./Dataset/GSE128477_Stochastic_origins_hg38.bed") # 256,600 origins
colnames(SNS.stochastic) <- c("chr", "start", "end", "ID")

# Ini-seq data from
# Guilbaud et al.  Nucleic Acids Research, 2022, https://doi.org/10.1093/nar/gkac555
# GSM5658908_Ini-seq2.called.replication.origins.bed
# Genomic location are in hg38

IniSeq.ori <- read.table("./Dataset/GSM5658908_Ini-seq2.called.replication.origins.bed", header = F) # 23,905 origins
colnames(IniSeq.ori) <- c("chr", "start", "end", "EFF")

############################################################################
# Define center of origins and assign ori id

SNS.core.center <- SNS.core %>% mutate(start=(as.integer(round(start+end)/2)), end = start+1, class = "core")
SNS.stochastic.center <- SNS.stochastic %>% mutate(start=(as.integer(round(start+end)/2)), end = start+1, class = "stochastic")
IniSeq.center <- IniSeq.ori %>% mutate(start=(as.integer(round(start+end)/2)), end = start+1, class = "iniseq", chr = paste("chr", chr, sep = "")) 

# Combine all SNS origins

SNS.ori.center <- rbind(SNS.core.center, SNS.stochastic.center)

# Save bed files

write.table(SNS.core.center, "./Figure_1/SNS.core.bed", sep="\t", col.names = F, row.names = F, quote = F)
write.table(SNS.stochastic.center, "./Figure_1/SNS.stochastic.bed", sep="\t", col.names = F, row.names = F, quote = F)
write.table(SNS.ori.center, "./Figure_1/SNS.all.bed", sep="\t", col.names = F, row.names = F, quote = F)
write.table(IniSeq.center, "./Figure_1/IniSeq.bed", sep="\t", col.names = F, row.names = F, quote = F)

# Compute distances between origins with bedtools closest

# sort -k1,1 -k2,2n SNS.core.bed > SNS.core.sorted.bed
# sort -k1,1 -k2,2n SNS.stochastic.bed > SNS.stochastic.sorted.bed
# sort -k1,1 -k2,2n SNS.all.bed > SNS.all.sorted.bed
# sort -k1,1 -k2,2n IniSeq.bed > IniSeq.sorted.bed
# 
# bedtools closest -a IniSeq.sorted.bed -b SNS.all.sorted.bed -d > Closest.IniSeq.SNS.all.bed
# bedtools closest -a IniSeq.sorted.bed -b SNS.core.sorted.bed -d > Closest.IniSeq.SNS.core.bed
# bedtools closest -a IniSeq.sorted.bed -b SNS.stochastic.sorted.bed -d > Closest.IniSeq.SNS.stochastic.bed

# Compute overlap

# IniSeq ori vs all SNS origins

Closest.IniSeq.SNS.all <- read.table("./Figure_1/Closest.IniSeq.SNS.all.bed") %>% # 23,918 entries
  dplyr::dplyr::select(V1, V2, V9, V10, V11)
colnames(Closest.IniSeq.SNS.all) <- c("chr", "ori.iniseq.pos", "closest.SNS.id", "SNS.class", "dist")
Closest.IniSeq.SNS.all <- Closest.IniSeq.SNS.all %>% mutate(iniseq.id = paste(chr, ori.iniseq.pos, sep = "_")) %>%
  distinct(iniseq.id, .keep_all = TRUE) # 23,905 entries

Closest.IniSeq.SNS.all.5kb <- Closest.IniSeq.SNS.all %>% filter(dist <= 5000) # 22,050 origins

# IniSeq ori vs SNS core

Closest.IniSeq.SNS.core <- read.table("./Figure_1/Closest.IniSeq.SNS.core.bed") %>% # 23,912 entries
  dplyr::dplyr::select(V1, V2, V9, V10, V11)
colnames(Closest.IniSeq.SNS.core) <- c("chr", "ori.iniseq.pos", "closest.SNS.id", "SNS.class", "dist")
Closest.IniSeq.SNS.core <- Closest.IniSeq.SNS.core %>% mutate(iniseq.id = paste(chr, ori.iniseq.pos, sep = "_")) %>%
  distinct(iniseq.id, .keep_all = TRUE) # 23,905 entries

Closest.IniSeq.SNS.core.5kb <- Closest.IniSeq.SNS.core %>% filter(dist <= 5000) # 17,406 origins --> 46,742 core origins remaining

# IniSeq ori vs SNS stochastic

Closest.IniSeq.SNS.stochastic <- read.table("./Figure_1/Closest.IniSeq.SNS.stochastic.bed") %>% # 23,911 entries
  dplyr::dplyr::select(V1, V2, V9, V10, V11)
colnames(Closest.IniSeq.SNS.stochastic) <- c("chr", "ori.iniseq.pos", "closest.SNS.id", "SNS.class", "dist")
Closest.IniSeq.SNS.stochastic <- Closest.IniSeq.SNS.stochastic %>% mutate(iniseq.id = paste(chr, ori.iniseq.pos, sep = "_")) %>%
  distinct(iniseq.id, .keep_all = TRUE) %>% # 23,905 entries, dplyr::select only origins that do not overlap with SNS core origins
  filter(iniseq.id %in% Closest.IniSeq.SNS.core.more.5kb$iniseq.id) # 6,499 origins

Closest.IniSeq.SNS.stochastic.5kb <- Closest.IniSeq.SNS.stochastic %>% filter(dist <= 5000) # 4,644 origins --> 251956 stochastic origins remaining

# Prepare donut plot

donut.df <- cbind.data.frame(category = c("stochastic", "core", "constituive"), count = c(251956,46742,23905))
donut.df$fraction = donut.df$count / sum(donut.df$count)
donut.df$ymax = cumsum(donut.df$fraction)
donut.df$ymin = c(0, head(donut.df$ymax, n=-1))
donut.plot <- donut.df %>% ggplot(aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) + geom_rect() + coord_polar(theta="y") + xlim(c(2, 4))

pdf("./Rplot/Fig_1A.pdf", width=6, height=4, useDingbats=FALSE)
donut.plot
dev.off()

############################################################################
# Define origin classes
# Remove Iniseq origins from SNS core
# and filter out origins with an inter-origins distance smaller than 20kb

IniSeq.ori.bedtools <- IniSeq.ori %>% mutate(chr = paste0("chr", chr))
write.table(IniSeq.ori.bedtools, "./Figure_1/IniSeq_ori.bed", sep="\t", col.names = F, row.names = F, quote = F)
# bedtools closest -a ./Dataset/GSE128477_Core_origins_hg38.bed -b ./Figure_1/IniSeq_ori.bed -d > ./Figure_1/Closest.SNS.core.IniSeq.bed
SNS.core.closest <- read.table("./Figure_1/Closest.SNS.core.IniSeq.bed") %>% filter(V9 <= 5000)
'%ni%' <- Negate('%in%')
SNS.core.unique <- SNS.core %>% filter(ID %ni% SNS.core.closest$V4)
write.table(SNS.core.unique, "./Figure_1/Core_origins_unique.bed", sep="\t", col.names = F, row.names = F, quote = F)
SNS.core.center <- SNS.core.unique %>% mutate(center = round((start+end)/2)) %>% 
  mutate(start=as.integer(center), end=as.integer(center+1)) %>% dplyr::select(chr, start, end, ID)
SNS.core.inter <- SNS.core.center %>% mutate(prev = start-lag(start), nxt = lead(start)-start) %>% 
  mutate(min = pmin(prev, nxt, na.rm = T)) %>% mutate(min = ifelse(min > 0, min, NA)) %>% 
  filter(min >= 20000) %>% dplyr::select(-prev, -nxt, -min)
write.table(SNS.core.inter, "./Figure_1/SNS.core.center.bed", sep="\t", col.names = F, row.names = F, quote = F)
# 12,224 SNS core origins

SNS.stochastic.center <- SNS.stochastic %>% mutate(center = round((start+end)/2)) %>% 
  mutate(start=as.integer(center), end=as.integer(center+1)) %>% dplyr::select(chr, start, end, ID)
SNS.stochastic.inter <- SNS.stochastic.center %>% mutate(prev = start-lag(start), nxt = lead(start)-start) %>% 
  mutate(min = pmin(prev, nxt, na.rm = T)) %>% mutate(min = ifelse(min > 0, min, NA)) %>% 
  filter(min >= 20000) %>% dplyr::select(-prev, -nxt, -min)
write.table(SNS.stochastic.inter, "./Figure_1/SNS.stochastic.center.bed", sep="\t", col.names = F, row.names = F, quote = F)
# 10,923 SNS stochastic origins

IniSeq.center <- IniSeq.ori %>% mutate(chr = paste0("chr", chr)) %>% mutate(center = round((start+end)/2)) %>% 
  mutate(start=as.integer(center), end=as.integer(center+1)) %>% dplyr::select(chr, start, end, EFF)
IniSeq.inter <- IniSeq.center %>% mutate(prev = start-lag(start), nxt = lead(start)-start) %>% 
  mutate(min = pmin(prev, nxt, na.rm = T)) %>% mutate(min = ifelse(min > 0, min, NA)) %>% 
  filter(min >= 20000) %>% dplyr::select(-prev, -nxt, -min)
write.table(IniSeq.inter, "./Figure_1/IniSeq.center.bed", sep="\t", col.names = F, row.names = F, quote = F)
# 9,351 IniSeq origins

# sort -k1,1 -k2,2n SNS.core.center.bed > SNS.core.center.isol.sorted.bed
# sort -k1,1 -k2,2n SNS.stochastic.center.bed > SNS.stochastic.center.isol.sorted.bed
# sort -k1,1 -k2,2n IniSeq.center.bed > IniSeq.center.isol.sorted.bed

############################################################################
############################################################################
###                   Replication timing analysis                        ###
############################################################################
############################################################################

############################################################################
# Repli-seq dataset from
# Sequencing newly replicated DNA reveals widespread plasticity in human replication timing
# Hansen et al. PNAS, 2010, 107, 139-144
# Data are available from Encode
# Recovered data related to four cell lines: GM12812 (LCLs), BG02 (hESCs), IMR-90 (lung fibroblast), MCF7 (breast adenocarcinoma)
# bigwig files corresponding to G1, S1, S4 and G2 phases of the cell cycle were downloaded for each cell lines
# bigwig files are in hg19 coordinates

# Replication timing score was computed as RT(score) = log2((G1+S1)/(S4+G2))
# The data used in Extended Figure 1a is considering the average of RT(score) over the four cell lines

# For BG02 cell line
BG02.G1 <- as.data.frame(import.bw("./Dataset/Repli-seq/BG02/BG02_G1B_ENCFF001GHF.bigWig"))
BG02.S1 <- as.data.frame(import.bw("./Dataset/Repli-seq/BG02/BG02_S1_ENCFF001GHM.bigWig"))
BG02.S4 <- as.data.frame(import.bw("./Dataset/Repli-seq/BG02/BG02_S4_ENCFF001GHV.bigWig"))
BG02.G2 <- as.data.frame(import.bw("./Dataset/Repli-seq/BG02/BG02_G2_ENCFF001GHH.bigWig"))
# Transform df
BG02.G1 <- BG02.G1 %>% mutate(G1 = score) %>% unite("bin", c(seqnames, start, end), sep = "_") %>% dplyr::select(bin, G1)
BG02.S1 <- BG02.S1 %>% mutate(S1 = score) %>% unite("bin", c(seqnames, start, end), sep = "_") %>% dplyr::select(bin, S1)
BG02.S4 <- BG02.S4 %>% mutate(S4 = score) %>% unite("bin", c(seqnames, start, end), sep = "_") %>% dplyr::select(bin, S4)
BG02.G2 <- BG02.G2 %>% mutate(G2 = score) %>% unite("bin", c(seqnames, start, end), sep = "_") %>% dplyr::select(bin, G2)
BG02.df <- BG02.G1 %>% right_join(BG02.S1, by = "bin") %>% right_join(BG02.S4, by = "bin") %>% right_join(BG02.G2, by = "bin")
BG02.df2 <- BG02.df %>% mutate(score = log2((G1+S1)/(S4+G2))) %>% replace(is.na(.), 0) %>%
  separate(bin, c("seqnames", "start", "end"), sep = "_") %>% dplyr::select(seqnames, start, end, score)
# Normalise score
BG02.df2.norm <- BG02.df2 %>% mutate(score = 2*score/((max(score)-min(score))))
# Create Grange object
BG02.gr <- makeGRangesFromDataFrame(BG02.df2.norm, keep.extra.columns=TRUE)
# Save bigwig in hg19
seqinfo.hg19 <- Seqinfo(genome="hg19")
seqinfo(BG02.gr) <- seqinfo.hg19
export(BG02.gr, con = "./Dataset/Repli-seq/BG02/BG02_replication_timing_hg19.bigWig", format = "bigWig")

# For GM12812 cell line
GM12812.G1 <- as.data.frame(import.bw("./Dataset/Repli-seq/GM12812/GM12812_G1B_ENCFF001GLS.bigWig"))
GM12812.S1 <- as.data.frame(import.bw("./Dataset/Repli-seq/GM12812/GM12812_S1_ENCFF001GLY.bigWig"))
GM12812.S4 <- as.data.frame(import.bw("./Dataset/Repli-seq/GM12812/GM12812_S4_ENCFF001GMG.bigWig"))
GM12812.G2 <- as.data.frame(import.bw("./Dataset/Repli-seq/GM12812/GM12812_G2_ENCFF001GLV.bigWig"))
# Transform df
GM12812.G1 <- GM12812.G1 %>% mutate(G1 = score) %>% unite("bin", c(seqnames, start, end), sep = "_") %>% dplyr::select(bin, G1)
GM12812.S1 <- GM12812.S1 %>% mutate(S1 = score) %>% unite("bin", c(seqnames, start, end), sep = "_") %>% dplyr::select(bin, S1)
GM12812.S4 <- GM12812.S4 %>% mutate(S4 = score) %>% unite("bin", c(seqnames, start, end), sep = "_") %>% dplyr::select(bin, S4)
GM12812.G2 <- GM12812.G2 %>% mutate(G2 = score) %>% unite("bin", c(seqnames, start, end), sep = "_") %>% dplyr::select(bin, G2)
GM12812.df <- GM12812.G1 %>% right_join(GM12812.S1, by = "bin") %>% right_join(GM12812.S4, by = "bin") %>% right_join(GM12812.G2, by = "bin")
GM12812.df2 <- GM12812.df %>% mutate(score = log2((G1+S1)/(S4+G2))) %>% replace(is.na(.), 0) %>%
  separate(bin, c("seqnames", "start", "end"), sep = "_") %>% dplyr::select(seqnames, start, end, score)
# Normalise score
GM12812.df2.norm <- GM12812.df2 %>% mutate(score = 2*score/((max(score)-min(score))))
# Create Grange object
GM12812.gr <- makeGRangesFromDataFrame(GM12812.df2.norm, keep.extra.columns=TRUE)
# Save bigwig in hg19
seqinfo.hg19 <- Seqinfo(genome="hg19")
seqinfo(GM12812.gr) <- seqinfo.hg19
export(GM12812.gr, con = "./Dataset/Repli-seq/GM12812/GM12812_replication_timing_hg19.bigWig", format = "bigWig")

# For IMR90 cell line
IMR90.G1 <- as.data.frame(import.bw("./Dataset/Repli-seq/IMR90/IMR90_G1B_ENCFF001GRA.bigWig"))
IMR90.S1 <- as.data.frame(import.bw("./Dataset/Repli-seq/IMR90/IMR90_S1_ENCFF001GRG.bigWig"))
IMR90.S4 <- as.data.frame(import.bw("./Dataset/Repli-seq/IMR90/IMR90_S4_ENCFF001GRQ.bigWig"))
IMR90.G2 <- as.data.frame(import.bw("./Dataset/Repli-seq/IMR90/IMR90_G2_ENCFF001GRD.bigWig"))
# Transform df
IMR90.G1 <- IMR90.G1 %>% mutate(G1 = score) %>% unite("bin", c(seqnames, start, end), sep = "_") %>% dplyr::select(bin, G1)
IMR90.S1 <- IMR90.S1 %>% mutate(S1 = score) %>% unite("bin", c(seqnames, start, end), sep = "_") %>% dplyr::select(bin, S1)
IMR90.S4 <- IMR90.S4 %>% mutate(S4 = score) %>% unite("bin", c(seqnames, start, end), sep = "_") %>% dplyr::select(bin, S4)
IMR90.G2 <- IMR90.G2 %>% mutate(G2 = score) %>% unite("bin", c(seqnames, start, end), sep = "_") %>% dplyr::select(bin, G2)
IMR90.df <- IMR90.G1 %>% right_join(IMR90.S1, by = "bin") %>% right_join(IMR90.S4, by = "bin") %>% right_join(IMR90.G2, by = "bin")
IMR90.df2 <- IMR90.df %>% mutate(score = log2((G1+S1)/(S4+G2))) %>% replace(is.na(.), 0) %>%
  separate(bin, c("seqnames", "start", "end"), sep = "_") %>% dplyr::select(seqnames, start, end, score)
# Normalise score
IMR90.df2.norm <- IMR90.df2 %>% mutate(score = 2*score/((max(score)-min(score))))
# Create Grange object
IMR90.gr <- makeGRangesFromDataFrame(IMR90.df2.norm, keep.extra.columns=TRUE)
# Save bigwig in hg19
seqinfo.hg19 <- Seqinfo(genome="hg19")
seqinfo(IMR90.gr) <- seqinfo.hg19
export(IMR90.gr, con = "./Dataset/Repli-seq/IMR90/IMR90_replication_timing_hg19.bigWig", format = "bigWig")

# For MCF7 cell line
MCF7.G1 <- as.data.frame(import.bw("./Dataset/Repli-seq/MCF7/MCF7_G1B_ENCFF001GSV.bigWig"))
MCF7.S1 <- as.data.frame(import.bw("./Dataset/Repli-seq/MCF7/MCF7_S1_ENCFF001GTD.bigWig"))
MCF7.S4 <- as.data.frame(import.bw("./Dataset/Repli-seq/MCF7/MCF7_S4_ENCFF001GTK.bigWig"))
MCF7.G2 <- as.data.frame(import.bw("./Dataset/Repli-seq/MCF7/MCF7_G2_ENCFF001GSX.bigWig"))
# Transform df
MCF7.G1 <- MCF7.G1 %>% mutate(G1 = score) %>% unite("bin", c(seqnames, start, end), sep = "_") %>% dplyr::select(bin, G1)
MCF7.S1 <- MCF7.S1 %>% mutate(S1 = score) %>% unite("bin", c(seqnames, start, end), sep = "_") %>% dplyr::select(bin, S1)
MCF7.S4 <- MCF7.S4 %>% mutate(S4 = score) %>% unite("bin", c(seqnames, start, end), sep = "_") %>% dplyr::select(bin, S4)
MCF7.G2 <- MCF7.G2 %>% mutate(G2 = score) %>% unite("bin", c(seqnames, start, end), sep = "_") %>% dplyr::select(bin, G2)
MCF7.df <- MCF7.G1 %>% right_join(MCF7.S1, by = "bin") %>% right_join(MCF7.S4, by = "bin") %>% right_join(MCF7.G2, by = "bin")
MCF7.df2 <- MCF7.df %>% mutate(score = log2((G1+S1)/(S4+G2))) %>% replace(is.na(.), 0) %>%
  separate(bin, c("seqnames", "start", "end"), sep = "_") %>% dplyr::select(seqnames, start, end, score)
# Normalise score
MCF7.df2.norm <- MCF7.df2 %>% mutate(score = 2*score/((max(score)-min(score))))
# Create Grange object
MCF7.gr <- makeGRangesFromDataFrame(MCF7.df2.norm, keep.extra.columns=TRUE)
# Save bigwig in hg19
seqinfo.hg19 <- Seqinfo(genome="hg19")
seqinfo(MCF7.gr) <- seqinfo.hg19
export(MCF7.gr, con = "./Dataset/Repli-seq/MCF7/MCF7_replication_timing_hg19.bigWig", format = "bigWig")

# Average Repli-seq score over the four cell lines
BG02.score <- as.data.frame(BG02.gr) %>% unite("bin", c(seqnames, start, end), sep = "_") %>% mutate(score_BG02 = score) %>% dplyr::select(bin, score_BG02)
GM12812.score <- as.data.frame(GM12812.gr) %>% unite("bin", c(seqnames, start, end), sep = "_") %>% mutate(score_GM12812 = score) %>% dplyr::select(bin, score_GM12812)
IMR90.score <- as.data.frame(IMR90.gr) %>% unite("bin", c(seqnames, start, end), sep = "_") %>% mutate(score_IMR90 = score) %>% dplyr::select(bin, score_IMR90)
MCF7.score <- as.data.frame(MCF7.gr) %>% unite("bin", c(seqnames, start, end), sep = "_") %>% mutate(score_MCF7 = score) %>% dplyr::select(bin, score_MCF7)
Repli.score <- BG02.score %>% right_join(GM12812.score) %>% right_join(IMR90.score) %>% right_join(MCF7.score)
Repli.score.avg <- Repli.score %>% mutate(score = (score_BG02+score_GM12812+score_IMR90+score_MCF7)/4) %>% replace(is.na(.), 0) %>% 
  separate(bin, c("seqnames", "start", "end"), sep = "_") %>% dplyr::select(seqnames, start, end, score)
# Normalise score
Repli.score.avg.norm <- Repli.score.avg %>% mutate(score = (2*score/((max(score)-min(score))))-median(score))
# Create Grange object
Repli.score.gr <- makeGRangesFromDataFrame(Repli.score.avg.norm, keep.extra.columns=TRUE)
# Save bigwig in hg19
seqinfo.hg19 <- Seqinfo(genome="hg19")
seqinfo(Repli.score.gr) <- seqinfo.hg19
export(Repli.score.gr, con = "./Dataset/Repli-seq/Repli_seq_score_hg19.bigWig", format = "bigWig")
# Modify file for LiftOver
Repli.score.for.liftover <- Repli.score.avg.norm %>% mutate(end = as.numeric(start) + 1)
# Save bed file
write.table(Repli.score.for.liftover, "./Dataset/Repli-seq/Repli_seq_score_for_liftover.bed", sep="\t", col.names = F, row.names = F, quote = F)
# Submit to UCSC liftover for generating Repli_seq_score_liftover_hg38.bed
# Sort
# sort -k1,1 -k2,2n Repli_seq_score_liftover_hg38.bed > Repli_seq_score_liftover_hg38.sorted.bed
# Open result
Repli.score.hg38 <- read.table("./Dataset/Repli-seq/Repli_seq_score_hg38.sorted.bed", header = F)
colnames(Repli.score.hg38) <- c("chr", "start", "end", "score")
Repli.score.hg38.bin <- cbind.data.frame(Repli.score.hg38, end2 = c(Repli.score.hg38$start[2:nrow(Repli.score.hg38)],0)) %>% 
  mutate(end = (end2-1)) %>% dplyr::select(chr, start, end, score)
# Save bedgraph
write.table(Repli.score.hg38.bin, "./Dataset/Repli-seq/Repli_seq_score_hg38.bedgraph", sep="\t", col.names = F, row.names = F, quote = F)
# Convert to bigWig
# bedGraphToBigWig Repli_seq_score_hg38.bedgraph hg38.chrom.sizes Repli_seq_score_hg38.bigWig
# manually clean bedgraph for entries with end smaller than start

############################################################################
# Distribution of origins within replication domain timing

# Compute replication timing in 100 kb zone using bedtools

# Split the hg38 genome into windows of 100 kb
# bedtools makewindows -g hg38.chrom.sizes -w 100000 > hg38_100kb.bed
# sort -k1,1 -k2,2n hg38_100kb.bed > hg38_100kb.sorted.bed
# Compute mean coverage for each window using bedGraph inputs
# bedtools map -a hg38_100kb.sorted.bed -b Repli_seq_score_hg38.bedgraph -c 4 -o mean -null 0 > Repli_seq_score_hg38_in_100kb_windows.bedgraph

# Open RepliSeq values in 100kb bins results
RepliSeq.mean.100kb <- read.table("./Dataset/Repli-seq/Repli_seq_score_hg38_in_100kb_windows.bedgraph")
colnames(RepliSeq.mean.100kb) <- c("chr", "start", "end", "timing")

# Compute the distribution of replication timing in bins
RepliSeq.mean.100kb <- RepliSeq.mean.100kb %>% mutate(chr = gsub("chr", "", chr)) %>% unite("bin", c(chr, start), sep = "_") %>%
  dplyr::dplyr::select(bin, timing) %>% filter(timing != 0)
RepliSeq.mean.100kb.bin <- RepliSeq.mean.100kb %>% mutate(timing.bin = (trunc(as.numeric(timing)*50))/50) %>% dplyr::dplyr::select(bin, timing.bin)

# Cont number of origins in bins of 100kb

SNS.core <- read.table("./Figure_1/Core_origins_unique.bed") # 64,148 origins
colnames(SNS.core) <- c("chr", "start", "end", "ID")

SNS.core.bin <- SNS.core %>% filter(ID %ni% SNS.core.closest$V4) %>% mutate(int = as.integer((trunc(start/100000))*100000)) %>% 
  unite("bin", chr, int, remove = FALSE)
SNS.core.bin.count <- as.data.frame(table(SNS.core.bin$bin)) 
colnames(SNS.core.bin.count) <- c("bin", "SNS.core")
SNS.core.bin.count <- SNS.core.bin.count %>% mutate(bin = gsub("chr", "", bin))

SNS.stochastic.bin <- SNS.stochastic %>% mutate(int = as.integer((trunc(start/100000))*100000)) %>% 
  unite("bin", chr, int, remove = FALSE)
SNS.stochastic.bin.count <- as.data.frame(table(SNS.stochastic.bin$bin))
colnames(SNS.stochastic.bin.count) <- c("bin", "SNS.stochastic")
SNS.stochastic.bin.count <- SNS.stochastic.bin.count %>% mutate(bin = gsub("chr", "", bin))

IniSeq.ori.bin <- as.data.frame(IniSeq.ori) %>% mutate(int = as.integer((trunc(start/100000))*100000)) %>% 
  unite("bin", chr, int, remove = FALSE)
IniSeq.ori.bin.count <- as.data.frame(table(IniSeq.ori.bin$bin))
colnames(IniSeq.ori.bin.count) <- c("bin", "IniSeq")

# Combine data

Ori.distribution <- RepliSeq.mean.100kb.bin %>% left_join(SNS.stochastic.bin.count) %>% left_join(SNS.core.bin.count) %>% left_join(IniSeq.ori.bin.count)
# 29514 bins

# Compute fraction of origins in each bins
all.sites <- vector()
fct.SNS.core <- vector()
fct.SNS.stochastic <- vector()
fct.IniSeq <- vector()
for (i in unique(Ori.distribution$timing.bin)) {
  Timing.df <- Ori.distribution %>% filter(timing.bin == i)
  p <- nrow(Timing.df)/nrow(Ori.distribution)
  q <- sum(Timing.df$SNS.stochastic, na.rm = TRUE)/256600
  r <- sum(Timing.df$SNS.core, na.rm = TRUE)/31523
  s <- sum(Timing.df$IniSeq, na.rm = TRUE)/23905
  all.sites <- append(all.sites, p)
  fct.SNS.core <- append(fct.SNS.core, q)
  fct.SNS.stochastic <- append(fct.SNS.stochastic, r)
  fct.IniSeq <- append(fct.IniSeq, s)
}
Ori.Timing.profile <- cbind.data.frame(timing.bin = unique(Ori.distribution$timing.bin), all.sites, fct.SNS.stochastic, fct.SNS.core, fct.IniSeq)
Ori.Timing.profile <- Ori.Timing.profile[!(Ori.Timing.profile$timing.bin == 0),]

# Prepare data to plot
df1 <- Ori.Timing.profile %>% dplyr::dplyr::select(timing.bin, fct = fct.SNS.stochastic)  %>% mutate(class = "stochastic")
df2 <- Ori.Timing.profile %>% dplyr::dplyr::select(timing.bin, fct = fct.SNS.core)  %>% mutate(class = "core")
df3 <- Ori.Timing.profile %>% dplyr::dplyr::select(timing.bin, fct = fct.IniSeq)  %>% mutate(class = "constitutive")
df4 <- Ori.Timing.profile %>% dplyr::dplyr::select(timing.bin, fct = all.sites)  %>% mutate(class = "all sites")
Plot.df <- rbind(df1, df2, df3, df4)

# Plot
pdf("./Rplot/Fig_S1A.pdf", width=6, height=4, useDingbats=FALSE)
Plot.df %>% ggplot(aes(x=timing.bin, y=fct, color = class)) +
  geom_point(size = 0.5) +
  geom_line() +
  xlab("Late <-                 Replication timing                 -> Early") + ylab("Fraction") +
  scale_color_manual(values=c("#999999", "#FC4E07", "#009E73", "#0072B2")) + xlim(-1,1) +
  theme_bw() + theme(aspect.ratio=0.4)
dev.off()

############################################################################
############################################################################
###                 Compute GC content of origins                        ###
############################################################################
############################################################################

############################################################################
# Analyse GC content of all origin classes (+- 1kb)

# Prepare genomic coordinates

SNS.core.1kb <- SNS.core %>% mutate(center = round((start+end)/2)) %>% mutate(start = as.integer(center-1000), end = as.integer(center+1000)) %>% 
  dplyr::dplyr::select(chr, start, end, ID)
SNS.stochastic.1kb <- SNS.stochastic %>% mutate(center = round((start+end)/2)) %>% mutate(start = as.integer(center-1000), end = as.integer(center+1000)) %>% 
  dplyr::dplyr::select(chr, start, end, ID)
IniSeq.ori.1kb <- IniSeq.ori %>% mutate(chr = paste0("chr", chr)) %>% mutate(center = round((start+end)/2)) %>% unite("ID", c(chr, start, end), sep = "_", remove = FALSE) %>% mutate(start = as.integer(center-1000), end = as.integer(center+1000)) %>% 
  dplyr::dplyr::select(chr, start, end, ID, EFF)

# Bin IniSeq origins by efficiency

quantile(IniSeq.ori$EFF, c(0.333,0.666))
#  33.3%  66.6% 
# 0.8418 0.9065 

IniSeq.ori.100bp <- IniSeq.ori %>% mutate(chr = paste0("chr", chr)) %>% mutate(center = round((start+end)/2)) %>% unite("ID", c(chr, start, end), sep = "_", remove = FALSE) %>% mutate(start = as.integer(center-100), end = as.integer(center+100)) %>% 
  dplyr::dplyr::select(chr, start, end, ID, EFF)

IniSeq.ori.low.100bp <- IniSeq.ori.100bp %>% filter(EFF < 0.8418)
IniSeq.ori.medium.100bp <- IniSeq.ori.100bp %>% filter(EFF >= 0.8418 & EFF < 0.9065)
IniSeq.ori.high.100bp <- IniSeq.ori.100bp %>% filter(EFF >= 0.9065)

# Load as Grange object
SNS.core.1kb.gr <- makeGRangesFromDataFrame(SNS.core.1kb, keep.extra.columns=T)
seqlevelsStyle(SNS.core.1kb.gr) <- "UCSC"
SNS.stochastic.1kb.gr <- makeGRangesFromDataFrame(SNS.stochastic.1kb, keep.extra.columns=T)
seqlevelsStyle(SNS.stochastic.1kb.gr) <- "UCSC"
IniSeq.ori.1kb.gr <- makeGRangesFromDataFrame(IniSeq.ori.1kb, keep.extra.columns=T)
seqlevelsStyle(IniSeq.ori.1kb.gr) <- "UCSC"
IniSeq.ori.low.100bp.gr <- makeGRangesFromDataFrame(IniSeq.ori.low.100bp, keep.extra.columns=T)
seqlevelsStyle(IniSeq.ori.low.100bp.gr) <- "UCSC"
IniSeq.ori.medium.100bp.gr <- makeGRangesFromDataFrame(IniSeq.ori.medium.100bp, keep.extra.columns=T)
seqlevelsStyle(IniSeq.ori.medium.100bp.gr) <- "UCSC"
IniSeq.ori.high.100bp.gr <- makeGRangesFromDataFrame(IniSeq.ori.high.100bp, keep.extra.columns=T)
seqlevelsStyle(IniSeq.ori.high.100bp.gr) <- "UCSC"

# Compute GC content

SNS.core.views <- Views(Hsapiens, SNS.core.1kb.gr)
gc.SNS.core <- letterFrequency(SNS.core.views, "GC", as.prob = TRUE)

SNS.stochastic.views <- Views(Hsapiens, SNS.stochastic.1kb.gr)
gc.SNS.stochastic <- letterFrequency(SNS.stochastic.views, "GC", as.prob = TRUE)

Ini.views <- Views(Hsapiens, IniSeq.ori.1kb.gr)
gc.Ini <- letterFrequency(Ini.views, "GC", as.prob = TRUE)

Ini.high.views <- Views(Hsapiens, IniSeq.ori.high.100bp.gr)
gc.Ini.high <- letterFrequency(Ini.high.views, "GC", as.prob = TRUE)

Ini.medium.views <- Views(Hsapiens, IniSeq.ori.medium.100bp.gr)
gc.Ini.medium <- letterFrequency(Ini.medium.views, "GC", as.prob = TRUE)

Ini.low.views <- Views(Hsapiens, IniSeq.ori.low.100bp.gr)
gc.Ini.low <- letterFrequency(Ini.low.views, "GC", as.prob = TRUE)

# Prepare plot

GC.df1 <- as.data.frame(gc.SNS.core) %>% mutate(class = "core")
GC.df2 <- as.data.frame(gc.SNS.stochastic) %>% mutate(class = "stochastic")
GC.df3 <- as.data.frame(gc.Ini) %>% mutate(class = "constitutive")
GC.plot <- rbind(GC.df1, GC.df2, GC.df3)

GC.plot.usage <- ggplot(GC.plot, aes(x=`G|C`, after_stat(density), fill=class, color=class)) +
  geom_histogram(alpha=0.5, position="identity", binwidth=0.01) +
  scale_color_manual(values=c("#FC4E07", "#009E73", "#0072B2")) + ggtitle("All p-value < 2.2e-16") +
  xlab("GC content (origin center +- 1kb)") +
  theme_bw() + theme(aspect.ratio=1)

ks.test(GC.df1$`G|C`, GC.df2$`G|C`)  # p-value < 2.2e-16
ks.test(GC.df1$`G|C`, GC.df3$`G|C`)  # p-value < 2.2e-16

GC.df3 <- as.data.frame(gc.Ini.low) %>% mutate(class = "low")
GC.df4 <- as.data.frame(gc.Ini.medium) %>% mutate(class = "medium")
GC.df5 <- as.data.frame(gc.Ini.high) %>% mutate(class = "high")
GC.plot.2 <- rbind(GC.df3, GC.df4, GC.df5)

GC.plot.eff <- ggplot(GC.plot.2, aes(x=`G|C`, after_stat(density), fill=class, color=class)) +
  geom_histogram(alpha=0.5, position="identity", binwidth=0.01) +
  scale_color_manual(values=c("#FC4E07", "#009E73", "#0072B2")) + ggtitle("All p-value < 2.2e-16") +
  xlab("GC content (origin center +- 100bp)") +
  theme_bw() + theme(aspect.ratio=1)

ks.test(GC.df3$`G|C`, GC.df4$`G|C`)  # p-value < 2.2e-16
ks.test(GC.df4$`G|C`, GC.df5$`G|C`)  # p-value < 2.2e-16

# Save plots

pdf("./Rplot/Fig_S1B_S1C.pdf", width=12, height=5, useDingbats=FALSE)
ggarrange(GC.plot.usage, GC.plot.eff, ncol = 2, nrow = 1)
dev.off()

############################################################################
############################################################################
###         Compute mutation rates at replication origins                ###
############################################################################
############################################################################

############################################################################
# DNM mutations are from the Gene4Denovo data-base
# All_De_novo_mutations_1.2.txt
# recovered from http://www.genemed.tech/gene4denovo/download
# data is in hg19
# Formatting

DNM <- fread("./Dataset/All_De_novo_mutations_1.2.txt", select = c(1:5))
DNM <- as.data.frame(DNM) # 670,082 mutations
colnames(DNM) <- c("chr", "start", "end", "REF", "ALT")
# Add mutation id
DNM <-  DNM %>% mutate(id = c(1:nrow(DNM))) %>% mutate(end = start+1) %>% mutate(chr = paste("chr", chr, sep = ""))
# Select SNVs
DNM.SNV <- DNM %>% filter(nchar(REF) == 1 & nchar(ALT) == 1 & REF != "-" & ALT != "-") # 613,230 SNVs
# Save bed files
write.table(DNM.SNV, "./Figure_1/DNM.SNV.hg19.bed", sep="\t", col.names = F, row.names = F, quote = F)
# Liftover to hg38
# sort -k1,1 -k2,2n DNM.SNV.hg19.bed > DNM.SNV.hg19.sorted.bed
# liftOver DNM.SNV.hg19.sorted.bed hg19ToHg38.over.chain DNM.SNV.hg38.bed DNM.SNV.unlifted.bed

# We show hereafter how we compute mutation rates from common SNPs
# The same procedure was used for DNMs

############################################################################
# All common germline variants from dbSNP build 151 on hg38
# common_all_20180418.vcf.gz
# recovered from ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/
# This dataset reports all variants representing alleles observed in the germline with 
# a minor allele frequency ≥ 0.01 in at least one 1000 Genomes Phase III major population, 
# with at least two individuals from different families having the same minor allele
# Use of vctools to select only SNPs
# vcftools --vcf common_all_20180418.vcf --remove-indels --recode --recode-INFO-all --out common_all_20180418.SNP.vcf
# SNPs together with the Variant Allele Origin (SAO) and allele frequency (CAF)
# with perl/bash
# perl -lane 'print join "\t",(@F[0,1,3,4], /(?:SAO)=[^;]+/g, /(?:CAF)=[^,]+/g)' common_all_20180418.SNP.vcf.recode.vcf > temp.txt
# sed 's/\CAF=//g' temp.txt > temp2.txt
# sed 's/\SAO=//g' temp2.txt > REDUCED_common_all_20180418.SNP.SAO.CAF.recode.vcf
# Use of vctools to dplyr::select only INDELs
# vcftools --vcf common_all_20180418.vcf --keep-only-indels --recode --recode-INFO-all --out common_all_20180418.INDEL.vcf
# awk '{print $1, $2, $4 ,$5}' common_all_20180418.INDEL.vcf.recode.vcf > REDUCED_common_all_20180418.INDEL.recode.vcf

population.SNP <- read.table("./Dataset/REDUCED_common_all_20180418.SNP.SAO.CAF.recode.vcf", header = F)
colnames(population.SNP) <- c("chr", "pos", "REF", "ALT", "SAO", "CAF")
nrow(population.SNP) # 33,629,539 SNPs

population.INDEL <- read.table("./Dataset/REDUCED_common_all_20180418.INDEL.recode.vcf", header = F)
colnames(population.INDEL) <- c("chr", "pos", "REF", "ALT")
nrow(population.INDEL) # 3,673,439 INDELs

# Prepare bed files of sequence variants for bedtools closest

population.SNP.bed <- population.SNP %>% mutate(chr = paste0("chr", chr), start = pos, end = as.integer(pos+1))  %>% dplyr::select(chr, start, end, REF, ALT, SAO, CAF)
population.INDEL.bed <- population.INDEL %>% mutate(chr = paste0("chr", chr), start = pos, end = as.integer(pos+1))  %>% dplyr::select(chr, start, end, REF, ALT)
write.table(population.SNP.bed, "./Figure_1/population_SNP.bed", sep="\t", col.names = F, row.names = F, quote = F)
write.table(population.INDEL.bed, "./Figure_1/population_INDEL.bed", sep="\t", col.names = F, row.names = F, quote = F)

# sort -k1,1 -k2,2n population_SNP.bed > population_SNP.sorted.bed
# sort -k1,1 -k2,2n population_INDEL.bed > population_INDEL.sorted.bed

############################################################################
# Compute the distance of each mutation to the closest origins

# bedtools closest -a population_SNP.sorted.bed -b SNS.core.center.isol.sorted.bed > population.SNP.SNS.core.closest.bed
# bedtools closest -a population_INDEL.sorted.bed -b SNS.core.center.isol.sorted.bed > population.INDEL.SNS.core.closest.bed
# 
# bedtools closest -a population_SNP.sorted.bed -b SNS.stochastic.center.isol.sorted.bed > population.SNP.SNS.stochastic.closest.bed
# bedtools closest -a population_INDEL.sorted.bed -b SNS.stochastic.center.isol.sorted.bed > population.INDEL.SNS.stochastic.closest.bed
# 
# bedtools closest -a population_SNP.sorted.bed -b IniSeq.center.isol.sorted.bed > population.SNP.IniSeq.closest.bed
# bedtools closest -a population_INDEL.sorted.bed -b IniSeq.center.isol.sorted.bed > population.INDEL.IniSeq.closest.bed

# Open results / dplyr::select mutations at a maximum of +-10kb from an origin and save results

population.SNP.SNS.core.closest <- read.table("./Figure_1/population.SNP.SNS.core.closest.bed") %>% 
  dplyr::select(chr = V1, pos = V2, REF = V4, ALT = V5, SAO = V6, CAF = V7, ori.pos = V9) %>% mutate(dist = pos - ori.pos) %>% 
  filter(dist >= -10000 & dist <= 10000)
population.SNP.SNS.stochastic.closest <- read.table("./Figure_1/population.SNP.SNS.stochastic.closest.bed") %>% 
  dplyr::select(chr = V1, pos = V2, REF = V4, ALT = V5, SAO = V6, CAF = V7, ori.pos = V9) %>% mutate(dist = pos - ori.pos)  %>% 
  filter(dist >= -10000 & dist <= 10000)
population.SNP.IniSeq.closest <- read.table("./Figure_1/population.SNP.IniSeq.closest.bed") %>% 
  dplyr::select(chr = V1, pos = V2, REF = V4, ALT = V5, SAO = V6, CAF = V7, ori.pos = V9, EFF = V11) %>% mutate(dist = pos - ori.pos)  %>% 
  filter(dist >= -10000 & dist <= 10000)

population.INDEL.SNS.core.closest <- read.table("./Figure_1/population.INDEL.SNS.core.closest.bed") %>% 
  dplyr::select(chr = V1, pos = V2, REF = V4, ALT = V5, ori.pos = V7) %>% mutate(dist = pos - ori.pos) %>% 
  filter(dist >= -10000 & dist <= 10000)
population.INDEL.SNS.stochastic.closest <- read.table("./Figure_1/population.INDEL.SNS.stochastic.closest.bed") %>% 
  dplyr::select(chr = V1, pos = V2, REF = V4, ALT = V5, ori.pos = V7) %>% mutate(dist = pos - ori.pos)  %>% 
  filter(dist >= -10000 & dist <= 10000)
population.INDEL.IniSeq.closest <- read.table("./Figure_1/population.INDEL.IniSeq.closest.bed") %>% 
  dplyr::select(chr = V1, pos = V2, REF = V4, ALT = V5, ori.pos = V7, EFF = V9) %>% mutate(dist = pos - ori.pos)  %>% 
  filter(dist >= -10000 & dist <= 10000)

# Save tables

write.table(population.SNP.SNS.core.closest, "./Figure_1/population.SNP.SNS.core.closest.10kb.bed", sep="\t", col.names = F, row.names = F, quote = F)
write.table(population.SNP.SNS.stochastic.closest, "./Figure_1/population.SNP.SNS.stochastic.closest.10kb.bed", sep="\t", col.names = F, row.names = F, quote = F)
write.table(population.SNP.IniSeq.closest, "./Figure_1/population.SNP.IniSeq.closest.10kb.bed", sep="\t", col.names = F, row.names = F, quote = F)
write.table(population.INDEL.SNS.core.closest, "./Figure_1/population.INDEL.SNS.core.closest.10kb.bed", sep="\t", col.names = F, row.names = F, quote = F)
write.table(population.INDEL.SNS.stochastic.closest, "./Figure_1/population.INDEL.SNS.stochastic.closest.10kb.bed", sep="\t", col.names = F, row.names = F, quote = F)
write.table(population.INDEL.IniSeq.closest, "./Figure_1/population.INDEL.IniSeq.closest.10kb.bed", sep="\t", col.names = F, row.names = F, quote = F)

# Open results

population.SNP.SNS.core.closest <- read.csv("./Figure_1/population.SNP.SNS.core.closest.10kb.bed", sep="\t")
population.SNP.SNS.stochastic.closest <- read.csv("./Figure_1/population.SNP.SNS.stochastic.closest.10kb.bed", sep="\t")
population.SNP.IniSeq.closest <- read.csv("./Figure_1/population.SNP.IniSeq.closest.10kb.bed", sep="\t")
population.INDEL.SNS.core.closest <- read.csv("./Figure_1/population.INDEL.SNS.core.closest.10kb.bed", sep="\t")
population.INDEL.SNS.stochastic.closest <- read.csv("./Figure_1/population.INDEL.SNS.stochastic.closest.10kb.bed", sep="\t")
population.INDEL.IniSeq.closest <- read.csv("./Figure_1/population.INDEL.IniSeq.closest.10kb.bed", sep="\t")
colnames(population.SNP.SNS.core.closest) <- c("chr", "pos", "REF", "ALT", "SAO", "CAF", "ori.pos", "dist")
colnames(population.SNP.SNS.stochastic.closest) <- c("chr", "pos", "REF", "ALT", "SAO", "CAF", "ori.pos", "dist")
colnames(population.SNP.IniSeq.closest) <- c("chr", "pos", "REF", "ALT", "SAO", "CAF", "ori.pos", "EFF", "dist")
colnames(population.INDEL.SNS.core.closest) <- c("chr", "pos", "REF", "ALT", "ori.pos", "dist")
colnames(population.INDEL.SNS.stochastic.closest) <- c("chr", "pos", "REF", "ALT", "ori.pos", "dist")
colnames(population.INDEL.IniSeq.closest) <- c("chr", "pos", "REF", "ALT", "ori.pos", "EFF", "dist")

############################################################################
# Compute base composition of selected origins for correcting mutation rates for local sequence composition effects

# Partition origin domains in 100 nucleotides bins

SNS.core.inter <- read.table("./Figure_1/SNS.core.center.bed")
SNS.stochastic.inter <- read.table("./Figure_1/SNS.stochastic.center.bed")
IniSeq.inter <- read.table("./Figure_1/IniSeq.center.bed")
colnames(SNS.core.inter) <- c("chr", "start", "end", "INFO")
colnames(SNS.stochastic.inter) <- c("chr", "start", "end", "INFO")
colnames(IniSeq.inter) <- c("chr", "start", "end", "INFO")

SNS.core.inter.10kb <- SNS.core.inter %>% mutate(start = as.integer(start-10000), end = as.integer(start+20000))
SNS.stochastic.inter.10kb <- SNS.stochastic.inter %>% mutate(start = as.integer(start-10000), end = as.integer(start+20000))
IniSeq.inter.10kb <- IniSeq.inter %>% mutate(start = as.integer(start-10000), end = as.integer(start+20000))
write.table(SNS.core.inter.10kb, "./Figure_1/SNS.core.inter.center.10kb.bed", sep="\t", col.names = F, row.names = F, quote = F)
write.table(SNS.stochastic.inter.10kb, "./Figure_1/SNS.stochastic.inter.center.10kb.bed", sep="\t", col.names = F, row.names = F, quote = F)
write.table(IniSeq.inter.10kb, "./Figure_1/IniSeq.inter.center.10kb.bed", sep="\t", col.names = F, row.names = F, quote = F)

# sort -k1,1 -k2,2n SNS.core.inter.center.10kb.bed > SNS.core.inter.center.10kb.sorted.bed
# sort -k1,1 -k2,2n SNS.stochastic.inter.center.10kb.bed > SNS.stochastic.inter.center.10kb.sorted.bed
# sort -k1,1 -k2,2n IniSeq.inter.center.10kb.bed > IniSeq.inter.center.10kb.sorted.bed
# 
# bedtools makewindows -b SNS.core.inter.center.10kb.sorted.bed -w 100 -i winnum > SNS.core.inter.10kb.100nt.split.bed
# bedtools makewindows -b SNS.stochastic.inter.center.10kb.sorted.bed -w 100 -i winnum > SNS.stochastic.inter.10kb.100nt.split.bed
# bedtools makewindows -b IniSeq.inter.center.10kb.sorted.bed -w 100 -i winnum > IniSeq.inter.10kb.100nt.split.bed

SNS.core.inter.10kb.gr <- import("./Figure_1/SNS.core.inter.10kb.100nt.split.bed")
SNS.stochastic.inter.10kb.gr <- import("./Figure_1/SNS.stochastic.inter.10kb.100nt.split.bed")
IniSeq.inter.10kb.gr <- import("./Figure_1/IniSeq.inter.10kb.100nt.split.bed")

SNS.core.inter.views <- Views(Hsapiens, SNS.core.inter.10kb.gr)
SNS.core.inter.bc <- letterFrequency(SNS.core.inter.views, c("A", "C", "G", "T"), as.prob = FALSE)
SNS.core.inter.bc.df <- cbind.data.frame(bin = SNS.core.inter.10kb.gr$name, SNS.core.inter.bc)
SNS.core.inter.summary <- SNS.core.inter.bc.df %>% group_by(bin) %>% summarise_at(c("A", "C", "G", "T"), sum) %>% mutate(dist = (as.numeric(bin)-100)*100, GC.skew = (G-C)/(G+C), AT.skew = (A-T)/(A+T)) %>% 
  mutate(A.freq = A/(A+C+T+G), C.freq = C/(A+C+T+G), G.freq = G/(A+C+T+G), T.freq = T/(A+C+T+G))

SNS.stochastic.inter.views <- Views(Hsapiens, SNS.stochastic.inter.10kb.gr)
SNS.stochastic.inter.bc <- letterFrequency(SNS.stochastic.inter.views, c("A", "C", "G", "T"), as.prob = FALSE)
SNS.stochastic.inter.bc.df <- cbind.data.frame(bin = SNS.stochastic.inter.10kb.gr$name, SNS.stochastic.inter.bc)
SNS.stochastic.inter.summary <- SNS.stochastic.inter.bc.df %>% group_by(bin) %>% summarise_at(c("A", "C", "G", "T"), sum) %>% mutate(dist = (as.numeric(bin)-100)*100, GC.skew = (G-C)/(G+C), AT.skew = (A-T)/(A+T)) %>% 
  mutate(A.freq = A/(A+C+T+G), C.freq = C/(A+C+T+G), G.freq = G/(A+C+T+G), T.freq = T/(A+C+T+G))

IniSeq.inter.views <- Views(Hsapiens, IniSeq.inter.10kb.gr)
IniSeq.inter.bc <- letterFrequency(IniSeq.inter.views, c("A", "C", "G", "T"), as.prob = FALSE)
IniSeq.inter.bc.df <- cbind.data.frame(bin = IniSeq.inter.10kb.gr$name, IniSeq.inter.bc)
IniSeq.inter.summary <- IniSeq.inter.bc.df %>% group_by(bin) %>% summarise_at(c("A", "C", "G", "T"), sum) %>% mutate(dist = (as.numeric(bin)-100)*100, GC.skew = (G-C)/(G+C), AT.skew = (A-T)/(A+T)) %>% 
  mutate(A.freq = A/(A+C+T+G), C.freq = C/(A+C+T+G), G.freq = G/(A+C+T+G), T.freq = T/(A+C+T+G))

write.csv(SNS.core.inter.summary, "./Figure_1/SNS.core.inter.10kb.base.composition.csv", quote = F, row.names = F)
write.csv(SNS.stochastic.inter.summary, "./Figure_1/SNS.stochastic.inter.10kb.base.composition.csv", quote = F, row.names = F)
write.csv(IniSeq.inter.summary, "./Figure_1/IniSeq.inter.10kb.base.composition.csv", quote = F, row.names = F)

SNS.core.inter.summary <- read.csv("./Figure_1/SNS.core.inter.10kb.base.composition.csv")
SNS.stochastic.inter.summary <- read.csv("./Figure_1/SNS.stochastic.inter.10kb.base.composition.csv")
IniSeq.inter.summary <- read.csv("./Figure_1/IniSeq.inter.10kb.base.composition.csv")

############################################################################
# Compute SNP mutation rates and mutation interdistances around origins
# Mutation rates are computed by counting the number of SNP divided by the number of considered origins
# Mutation interdistances are computed as the distance between two mutations and considering the mean values of the distance in each bin

# SNS core | 12224 origins

nrow(population.SNP.SNS.core.closest) # 2937612
population.SNP.SNS.core.closest.bin <- population.SNP.SNS.core.closest %>% 
  mutate(prev = pos-lag(pos), nxt = lead(pos)-pos) %>% mutate(dist.SNP = pmin(prev, nxt, na.rm = T)) %>% mutate(dist.SNP = ifelse(dist.SNP > 0 & dist.SNP < 100, dist.SNP, NA)) %>%  
  mutate(dist = (as.numeric(cut(dist, breaks = seq(-10000, 10000, 100)))-100)*100) %>% 
  mutate(mut.type = ifelse((REF == "G" | REF == "C"), "GC", "AT")) %>% group_by(dist) %>% 
  summarise(GC.mut = sum(mut.type == "GC"), AT.mut = sum(mut.type == "AT"), dist.SNP.mean = mean(dist.SNP, na.rm = T)) %>%
  right_join(SNS.core.inter.summary) %>% mutate(mut.rate = (GC.mut+AT.mut)/12224) %>%
  mutate(mut.rate.fold = log2(mut.rate/mean(mut.rate[c(1:20,180:200)])))

# SNS stochastic | 10923  origins

nrow(population.SNP.SNS.stochastic.closest) # 2505633
population.SNP.SNS.stochastic.closest.bin <- population.SNP.SNS.stochastic.closest %>%
  mutate(prev = pos-lag(pos), nxt = lead(pos)-pos) %>% mutate(dist.SNP = pmin(prev, nxt, na.rm = T)) %>% mutate(dist.SNP = ifelse(dist.SNP > 0 & dist.SNP < 100, dist.SNP, NA)) %>%  
  mutate(dist = (as.numeric(cut(dist, breaks = seq(-10000, 10000, 100)))-100)*100) %>% 
  mutate(mut.type = ifelse((REF == "G" | REF == "C"), "GC", "AT")) %>% group_by(dist) %>%
  summarise(GC.mut = sum(mut.type == "GC"), AT.mut = sum(mut.type == "AT"), dist.SNP.mean = mean(dist.SNP, na.rm = T)) %>%
  right_join(SNS.stochastic.inter.summary) %>% mutate(mut.rate = (GC.mut+AT.mut)/10923) %>%
  mutate(mut.rate.fold = log2(mut.rate/mean(mut.rate[c(1:20,180:200)])))

# IniSeq | 9351 origins

nrow(population.SNP.IniSeq.closest) # 2171614
population.SNP.IniSeq.closest.bin <- population.SNP.IniSeq.closest %>% 
  mutate(prev = pos-lag(pos), nxt = lead(pos)-pos) %>% mutate(dist.SNP = pmin(prev, nxt, na.rm = T)) %>% mutate(dist.SNP = ifelse(dist.SNP > 0 & dist.SNP < 100, dist.SNP, NA)) %>%  
  mutate(dist = (as.numeric(cut(dist, breaks = seq(-10000, 10000, 100)))-100)*100) %>% 
  mutate(mut.type = ifelse((REF == "G" | REF == "C"), "GC", "AT")) %>% group_by(dist) %>% 
  summarise(GC.mut = sum(mut.type == "GC"), AT.mut = sum(mut.type == "AT"), dist.SNP.mean = mean(dist.SNP, na.rm = T)) %>%
  right_join(IniSeq.inter.summary) %>% mutate(mut.rate = (GC.mut+AT.mut)/9351) %>%
  mutate(mut.rate.fold = log2(mut.rate/mean(mut.rate[c(1:20,180:200)])))

# Prepare plot for mutation rates

population.SNP.plot.df1 <- population.SNP.SNS.stochastic.closest.bin %>% dplyr::select(dist,mut.rate.fold) %>% mutate(class = "SNS stochastic")
population.SNP.plot.df2 <- population.SNP.SNS.core.closest.bin %>% dplyr::select(dist,mut.rate.fold) %>% mutate(class = "SNS core")
population.SNP.plot.df3 <- population.SNP.IniSeq.closest.bin %>% dplyr::select(dist,mut.rate.fold) %>% mutate(class = "IniSeq")
population.SNP.plot <- rbind(population.SNP.plot.df1, population.SNP.plot.df2, population.SNP.plot.df3)

population.SNP.mut.rate.plot <- ggplot(population.SNP.plot, aes(x=dist, y=mut.rate.fold, color=class)) +
  geom_line(size = 0.5) +
  scale_color_manual(values=c("#FC4E07", "#009E73", "#0072B2")) +
  xlab("Distance from origin (bp)") + ylab("Mutation rate (log2 FC)") + ylim(-0.11,0.2) +
  theme_bw() + theme(aspect.ratio=0.5)

pdf("./Rplot/Fig_1D.pdf", width=6, height=3, useDingbats=FALSE)
population.SNP.mut.rate.plot
dev.off()

population.SNP.plot.df4 <- population.SNP.SNS.stochastic.closest.bin %>% dplyr::select(dist,dist.SNP.mean) %>% mutate(class = "SNS stochastic")
population.SNP.plot.df5 <- population.SNP.SNS.core.closest.bin %>% dplyr::select(dist,dist.SNP.mean) %>% mutate(class = "SNS core")
population.SNP.plot.df6 <- population.SNP.IniSeq.closest.bin %>% dplyr::select(dist,dist.SNP.mean) %>% mutate(class = "IniSeq")
population.SNP.plot.2 <- rbind(population.SNP.plot.df4, population.SNP.plot.df5, population.SNP.plot.df6)

population.SNP.mut.interdistance.plot <- ggplot(population.SNP.plot.2, aes(x=dist, y=dist.SNP.mean, color=class)) +
  geom_line(size = 0.5) +
  scale_color_manual(values=c("#FC4E07", "#009E73", "#0072B2")) +
  xlab("Distance from origin (bp)") + ylab("Mutation interdistance (nt)") +
  theme_bw() + theme(aspect.ratio=0.5)

pdf("./Rplot/Figure_S1D.pdf", width=6, height=3, useDingbats=FALSE)
population.SNP.mut.interdistance.plot
dev.off()

# Statistics

# For SNP

ori.int <- c(-10:10)*100
background <- setdiff(c(-100:100)*100, ori.int)

population.SNP.IniSeq.closest.stat <- population.SNP.IniSeq.closest %>% mutate(bin = round(dist/100)*100) %>% mutate(stat = case_when(bin %in% ori.int ~ "in", TRUE ~ "out" ))
df <- table(population.SNP.IniSeq.closest.stat$stat)/c(length(ori.int),length(background))
barplot(df)
chisq.test(df)$p.value # 0.0002467

population.SNP.SNS.core.closest.stat <- population.SNP.SNS.core.closest %>% mutate(bin = round(dist/100)*100) %>% mutate(stat = case_when(bin %in% ori.int ~ "in", TRUE ~ "out" ))
df <- table(population.SNP.SNS.core.closest.stat$stat)/c(length(ori.int),length(background))
barplot(df)
chisq.test(df)$p.value # 3.49724e-08

population.SNP.SNS.stochastic.closest.stat <- population.SNP.SNS.stochastic.closest %>% mutate(bin = round(dist/100)*100) %>% mutate(stat = case_when(bin %in% ori.int ~ "in", TRUE ~ "out" ))
df <- table(population.SNP.SNS.stochastic.closest.stat$stat)/c(length(ori.int),length(background))
barplot(df)
chisq.test(df)$p.value # 0.01438375

# For mutation interdistances

population.SNP.IniSeq.closest.stat <- population.SNP.IniSeq.closest.bin %>% mutate(stat = case_when(dist %in% ori.int ~ "in", TRUE ~ "out" ), inter = round(dist.SNP.mean*10))
df <- table(population.SNP.IniSeq.closest.stat$stat,population.SNP.IniSeq.closest.stat$inter)
barplot(df)
chisq.test(df)$p.value # 1.679107e-08

population.SNP.SNS.core.closest.stat <- population.SNP.SNS.core.closest.bin %>% mutate(stat = case_when(dist %in% ori.int ~ "in", TRUE ~ "out" ), inter = round(dist.SNP.mean*10))
df <- table(population.SNP.SNS.core.closest.stat$stat,population.SNP.SNS.core.closest.stat$inter)
barplot(df)
chisq.test(df)$p.value # 7.795784e-24

population.SNP.SNS.stochastic.closest.stat <- population.SNP.SNS.stochastic.closest.bin %>% mutate(stat = case_when(dist %in% ori.int ~ "in", TRUE ~ "out" ), inter = round(dist.SNP.mean*10))
df <- table(population.SNP.SNS.stochastic.closest.stat$stat,population.SNP.SNS.stochastic.closest.stat$inter)
barplot(df)
chisq.test(df)$p.value # 1.145846e-05

############################################################################
# Compute INDEL mutation rates around origins
# Mutation rates are computed by counting the number of INDEL divided by the number of considered origins
# Consider short (= 1bp) and long (>= 2bp) individually

# IniSeq | 9351 origins
nrow(population.INDEL.IniSeq.closest) # 259535
population.INDEL.IniSeq.closest.bin <- population.INDEL.IniSeq.closest %>%
  mutate(mut.base = substring(REF, 1, 1), mut.type = ifelse(nchar(REF) == 1, "INS", "DEL")) %>%
  mutate(mut.base = ifelse((mut.base == "G" | mut.base == "C"), "GC", "AT")) %>%
  mutate(mut.length = ifelse((mut.type == "INS" & nchar(ALT) == 2) | (mut.type == "DEL" & nchar(REF) == 2), "short", "long")) %>%
  mutate(dist = (as.numeric(cut(dist, breaks = seq(-10000, 10000, 100)))-100)*100) 
# Compute statistics for long or short INDEL events
population.INDEL.IniSeq.closest.bin.short <- population.INDEL.IniSeq.closest.bin %>% filter(mut.length == "short") %>%
  group_by(dist) %>% summarise(GC.mut = sum(mut.base == "GC"), AT.mut = sum(mut.base == "AT")) %>% right_join(IniSeq.inter.summary, by = "dist") %>%
  mutate(mut.rate = (GC.mut+AT.mut)/9351) %>% 
  mutate(mut.rate.fold = log2(mut.rate/mean(mut.rate[c(1:20,180:200)])))
population.INDEL.IniSeq.closest.bin.long <- population.INDEL.IniSeq.closest.bin %>% filter(mut.length == "long") %>%
  group_by(dist) %>% summarise(GC.mut = sum(mut.base == "GC"), AT.mut = sum(mut.base == "AT")) %>% right_join(IniSeq.inter.summary, by = "dist") %>%
  mutate(mut.rate = (GC.mut+AT.mut)/9351) %>% 
  mutate(mut.rate.fold = log2(mut.rate/mean(mut.rate[c(1:20,180:200)])))

# SNS core | 12224 origins
nrow(population.INDEL.SNS.core.closest) # 321213
population.INDEL.SNS.core.closest.bin <- population.INDEL.SNS.core.closest %>%
  mutate(mut.base = substring(REF, 1, 1), mut.type = ifelse(nchar(REF) == 1, "INS", "DEL")) %>%
  mutate(mut.base = ifelse((mut.base == "G" | mut.base == "C"), "GC", "AT")) %>%
  mutate(mut.length = ifelse((mut.type == "INS" & nchar(ALT) == 2) | (mut.type == "DEL" & nchar(REF) == 2), "short", "long")) %>%
  mutate(dist = (as.numeric(cut(dist, breaks = seq(-10000, 10000, 100)))-100)*100) 
# Compute statistics for long or short INDEL events
population.INDEL.SNS.core.closest.bin.short <- population.INDEL.SNS.core.closest.bin %>% filter(mut.length == "short") %>%
  group_by(dist) %>% summarise(GC.mut = sum(mut.base == "GC"), AT.mut = sum(mut.base == "AT")) %>% right_join(SNS.core.inter.summary, by = "dist") %>%
  mutate(mut.rate = (GC.mut+AT.mut)/12224) %>% 
  mutate(mut.rate.fold = log2(mut.rate/mean(mut.rate[c(1:20,180:200)])))
population.INDEL.SNS.core.closest.bin.long <- population.INDEL.SNS.core.closest.bin %>% filter(mut.length == "long") %>%
  group_by(dist) %>% summarise(GC.mut = sum(mut.base == "GC"), AT.mut = sum(mut.base == "AT")) %>% right_join(SNS.core.inter.summary, by = "dist") %>%
  mutate(mut.rate = (GC.mut+AT.mut)/12224) %>% 
  mutate(mut.rate.fold = log2(mut.rate/mean(mut.rate[c(1:20,180:200)])))

# SNS stochastic | 10923 origins
nrow(population.INDEL.SNS.stochastic.closest) # 276499
population.INDEL.SNS.stochastic.closest.bin <- population.INDEL.SNS.stochastic.closest %>%
  mutate(mut.base = substring(REF, 1, 1), mut.type = ifelse(nchar(REF) == 1, "INS", "DEL")) %>%
  mutate(mut.base = ifelse((mut.base == "G" | mut.base == "C"), "GC", "AT")) %>%
  mutate(mut.length = ifelse((mut.type == "INS" & nchar(ALT) == 2) | (mut.type == "DEL" & nchar(REF) == 2), "short", "long")) %>%
  mutate(dist = (as.numeric(cut(dist, breaks = seq(-10000, 10000, 100)))-100)*100) 
# Compute statistics for long or short INDEL events
population.INDEL.SNS.stochastic.closest.bin.short <- population.INDEL.SNS.stochastic.closest.bin %>% filter(mut.length == "short") %>%
  group_by(dist) %>% summarise(GC.mut = sum(mut.base == "GC"), AT.mut = sum(mut.base == "AT")) %>% right_join(SNS.stochastic.inter.summary, by = "dist") %>%
  mutate(mut.rate = (GC.mut+AT.mut)/12224) %>% 
  mutate(mut.rate.fold = log2(mut.rate/mean(mut.rate[c(1:20,180:200)])))
population.INDEL.SNS.stochastic.closest.bin.long <- population.INDEL.SNS.stochastic.closest.bin %>% filter(mut.length == "long") %>%
  group_by(dist) %>% summarise(GC.mut = sum(mut.base == "GC"), AT.mut = sum(mut.base == "AT")) %>% right_join(SNS.stochastic.inter.summary, by = "dist") %>%
  mutate(mut.rate = (GC.mut+AT.mut)/12224) %>% 
  mutate(mut.rate.fold = log2(mut.rate/mean(mut.rate[c(1:20,180:200)])))

# Prepare plot for INDEL mutation rates

population.INDEL.Iniseq.plot.df1 <- population.INDEL.IniSeq.closest.bin.short %>% dplyr::select(dist, rate = mut.rate.fold) %>% mutate(class = "Short INDELs")
population.INDEL.Iniseq.plot.df2 <- population.INDEL.IniSeq.closest.bin.long %>% dplyr::select(dist, rate = mut.rate.fold) %>% mutate(class = "Long INDELs")
population.INDEL.Iniseq.plot <- rbind(population.INDEL.Iniseq.plot.df1, population.INDEL.Iniseq.plot.df2)

a <- ggplot(population.INDEL.Iniseq.plot, aes(x=dist, y=rate, color=class)) +
  geom_line(size = 0.5) + ylim(-0.5, 0.5) +
  scale_color_manual(values=c("#FC4E07", "#0072B2")) +
  xlab("Distance from origin (bp)") + ylab("Mutation rate (log2 FC)") + ggtitle("Population INDELs at IniSeq origins") +
  theme_bw() + theme(aspect.ratio=0.5)

population.INDEL.SNS.core.plot.df1 <- population.INDEL.SNS.core.closest.bin.short %>% dplyr::select(dist, rate = mut.rate.fold) %>% mutate(class = "Short INDELs")
population.INDEL.SNS.core.plot.df2 <- population.INDEL.SNS.core.closest.bin.long %>% dplyr::select(dist, rate = mut.rate.fold) %>% mutate(class = "Long INDELs")
population.INDEL.SNS.core.plot <- rbind(population.INDEL.SNS.core.plot.df1, population.INDEL.SNS.core.plot.df2)

b <- ggplot(population.INDEL.SNS.core.plot, aes(x=dist, y=rate, color=class)) +
  geom_line(size = 0.5) + ylim(-0.5, 0.5) +
  scale_color_manual(values=c("#FC4E07", "#0072B2")) +
  xlab("Distance from origin (bp)") + ylab("Mutation rate (log2 FC)") + ggtitle("Population INDELs at SNS core origins") +
  theme_bw() + theme(aspect.ratio=0.5)

population.INDEL.SNS.stochastic.plot.df1 <- population.INDEL.SNS.stochastic.closest.bin.short %>% dplyr::select(dist, rate = mut.rate.fold) %>% mutate(class = "Short INDELs")
population.INDEL.SNS.stochastic.plot.df2 <- population.INDEL.SNS.stochastic.closest.bin.long %>% dplyr::select(dist, rate = mut.rate.fold) %>% mutate(class = "Long INDELs")
population.INDEL.SNS.stochastic.plot <- rbind(population.INDEL.SNS.stochastic.plot.df1, population.INDEL.SNS.stochastic.plot.df2)

c <- ggplot(population.INDEL.SNS.stochastic.plot, aes(x=dist, y=rate, color=class)) +
  geom_line(size = 0.5) + ylim(-0.5, 0.5) +
  scale_color_manual(values=c("#FC4E07", "#0072B2")) +
  xlab("Distance from origin (bp)") + ylab("Mutation rate (log2 FC)") + ggtitle("Population INDELs at SNS stochastic origins") +
  theme_bw() + theme(aspect.ratio=0.5)

pdf("./Rplot/Figure_S1E.pdf", width=8, height=7, useDingbats=FALSE)
ggarrange(a, b, c, ncol = 1, nrow = 3)
dev.off()

# Statistics

ori.int <- c(-2:2)*100
background <- c(-100:-50,50:100)*100

population.INDEL.IniSeq.closest.stat <- population.INDEL.IniSeq.closest.bin %>% filter(mut.length == "long") %>% mutate(stat = case_when(dist %in% ori.int ~ "in", dist %in% background ~ "out", TRUE ~ "NA")) %>% filter(stat == "in" | stat == "out")
df <- table(population.INDEL.IniSeq.closest.stat$stat)/c(length(ori.int),length(background))
barplot(df)
chisq.test(df)$p.value # 2.816039e-05

population.INDEL.SNS.core.closest.stat <- population.INDEL.SNS.core.closest.bin %>% filter(mut.length == "long") %>% mutate(stat = case_when(dist %in% ori.int ~ "in", dist %in% background ~ "out", TRUE ~ "NA")) %>% filter(stat == "in" | stat == "out")
df <- table(population.INDEL.SNS.core.closest.stat$stat)/c(length(ori.int),length(background))
barplot(df)
chisq.test(df)$p.value # 0.0003726243

population.INDEL.SNS.stochastic.closest.stat <- population.INDEL.SNS.stochastic.closest.bin %>% filter(mut.length == "long") %>% mutate(stat = case_when(dist %in% ori.int ~ "in", dist %in% background ~ "out", TRUE ~ "NA")) %>% filter(stat == "in" | stat == "out")
df <- table(population.INDEL.SNS.stochastic.closest.stat$stat)/c(length(ori.int),length(background))
barplot(df)
chisq.test(df)$p.value # 0.001315007

############################################################################
# Analyse population SNP densities (corrected for base composition) around origins

# At IniSeq origins

population.SNP.IniSeq.closest.pyr <- population.SNP.IniSeq.closest %>% filter(nchar(ALT) == 1) %>% mutate(dist = (as.numeric(cut(dist, breaks = seq(-10000, 10000, 100)))-100)*100) %>%
  mutate(class = paste0(REF, ">", ALT)) %>%
  mutate(class = gsub("G>A", "C>T", class), class = gsub("G>C", "C>G", class), class = gsub("G>T", "C>A", class), class = gsub("A>G", "T>C", class),class = gsub("A>C", "T>G", class), class = gsub("A>T", "T>A", class)) %>% 
  group_by(dist) %>% summarise(CT.mut = sum(class == "C>T"), CG.mut = sum(class == "C>G") , CA.mut = sum(class == "C>A"), TC.mut = sum(class == "T>C"), TG.mut = sum(class == "T>G") , TA.mut = sum(class == "T>A")) %>% 
  right_join(IniSeq.inter.summary) %>% mutate(CT.mut.cor = CT.mut/(G+C), CG.mut.cor = CG.mut/(G+C), CA.mut.cor = CA.mut/(G+C), TC.mut.cor = TC.mut/(T+A), TG.mut.cor = TG.mut/(T+A), TA.mut.cor = TA.mut/(T+A)) %>% 
  mutate(CT.mut.fold = log2(CT.mut.cor/mean(CT.mut.cor[c(1:20,180:200)])), CG.mut.fold = log2(CG.mut.cor/mean(CG.mut.cor[c(1:20,180:200)])), CA.mut.fold = log2(CA.mut.cor/mean(CA.mut.cor[c(1:20,180:200)])),
         TC.mut.fold = log2(TC.mut.cor/mean(TC.mut.cor[c(1:20,180:200)])), TG.mut.fold = log2(TG.mut.cor/mean(TG.mut.cor[c(1:20,180:200)])), TA.mut.fold = log2(TA.mut.cor/mean(TA.mut.cor[c(1:20,180:200)])))

# Prepare plot
population.SNP.IniSeq.closest.CT <- population.SNP.IniSeq.closest.pyr %>% dplyr::select(dist, mut = CT.mut.fold) %>% mutate(group = "C>T")
population.SNP.IniSeq.closest.CG <- population.SNP.IniSeq.closest.pyr %>% dplyr::select(dist, mut = CG.mut.fold) %>% mutate(group = "C>G")
population.SNP.IniSeq.closest.CA <- population.SNP.IniSeq.closest.pyr %>% dplyr::select(dist, mut = CA.mut.fold) %>% mutate(group = "C>A")
population.SNP.IniSeq.closest.TC <- population.SNP.IniSeq.closest.pyr %>% dplyr::select(dist, mut = TC.mut.fold) %>% mutate(group = "T>C")
population.SNP.IniSeq.closest.TG <- population.SNP.IniSeq.closest.pyr %>% dplyr::select(dist, mut = TG.mut.fold) %>% mutate(group = "T>G")
population.SNP.IniSeq.closest.TA <- population.SNP.IniSeq.closest.pyr %>% dplyr::select(dist, mut = TA.mut.fold) %>% mutate(group = "T>A")
population.SNP.IniSeq.closest.df <- rbind(population.SNP.IniSeq.closest.CT, population.SNP.IniSeq.closest.CG, population.SNP.IniSeq.closest.CA,
                                          population.SNP.IniSeq.closest.TC, population.SNP.IniSeq.closest.TG, population.SNP.IniSeq.closest.TA)

d <- ggplot(population.SNP.IniSeq.closest.df, aes(x=dist, y=mut, color=group)) +
  geom_line(size = 0.5) + ylim(-0.6,0.6) +
  scale_color_manual(values=c("#E69F00", "#56B4E9", "#009E73","#CC79A7", "#0072B2", "#D55E00")) +
  xlab("Distance from origin (bp)") + ylab("Mutation density (log2 FC)") + ggtitle("Population SNP / IniSeq") +
  theme_bw() + theme(aspect.ratio=0.5)

# At SNS core origins

population.SNP.SNS.core.closest.pyr <- population.SNP.SNS.core.closest %>% filter(nchar(ALT) == 1) %>% mutate(dist = (as.numeric(cut(dist, breaks = seq(-10000, 10000, 100)))-100)*100) %>%
  mutate(class = paste0(REF, ">", ALT)) %>%
  mutate(class = gsub("G>A", "C>T", class), class = gsub("G>C", "C>G", class), class = gsub("G>T", "C>A", class), class = gsub("A>G", "T>C", class),class = gsub("A>C", "T>G", class), class = gsub("A>T", "T>A", class)) %>% 
  group_by(dist) %>% summarise(CT.mut = sum(class == "C>T"), CG.mut = sum(class == "C>G") , CA.mut = sum(class == "C>A"), TC.mut = sum(class == "T>C"), TG.mut = sum(class == "T>G") , TA.mut = sum(class == "T>A")) %>% 
  right_join(SNS.core.inter.summary) %>% mutate(CT.mut.cor = CT.mut/(G+C), CG.mut.cor = CG.mut/(G+C), CA.mut.cor = CA.mut/(G+C), TC.mut.cor = TC.mut/(T+A), TG.mut.cor = TG.mut/(T+A), TA.mut.cor = TA.mut/(T+A)) %>% 
  mutate(CT.mut.fold = log2(CT.mut.cor/mean(CT.mut.cor[c(1:20,180:200)])), CG.mut.fold = log2(CG.mut.cor/mean(CG.mut.cor[c(1:20,180:200)])), CA.mut.fold = log2(CA.mut.cor/mean(CA.mut.cor[c(1:20,180:200)])),
         TC.mut.fold = log2(TC.mut.cor/mean(TC.mut.cor[c(1:20,180:200)])), TG.mut.fold = log2(TG.mut.cor/mean(TG.mut.cor[c(1:20,180:200)])), TA.mut.fold = log2(TA.mut.cor/mean(TA.mut.cor[c(1:20,180:200)])))

# Prepare plot
population.SNP.SNS.core.closest.CT <- population.SNP.SNS.core.closest.pyr %>% dplyr::select(dist, mut = CT.mut.fold) %>% mutate(group = "C>T")
population.SNP.SNS.core.closest.CG <- population.SNP.SNS.core.closest.pyr %>% dplyr::select(dist, mut = CG.mut.fold) %>% mutate(group = "C>G")
population.SNP.SNS.core.closest.CA <- population.SNP.SNS.core.closest.pyr %>% dplyr::select(dist, mut = CA.mut.fold) %>% mutate(group = "C>A")
population.SNP.SNS.core.closest.TC <- population.SNP.SNS.core.closest.pyr %>% dplyr::select(dist, mut = TC.mut.fold) %>% mutate(group = "T>C")
population.SNP.SNS.core.closest.TG <- population.SNP.SNS.core.closest.pyr %>% dplyr::select(dist, mut = TG.mut.fold) %>% mutate(group = "T>G")
population.SNP.SNS.core.closest.TA <- population.SNP.SNS.core.closest.pyr %>% dplyr::select(dist, mut = TA.mut.fold) %>% mutate(group = "T>A")
population.SNP.SNS.core.closest.df <- rbind(population.SNP.SNS.core.closest.CT, population.SNP.SNS.core.closest.CG, population.SNP.SNS.core.closest.CA,
                                            population.SNP.SNS.core.closest.TC, population.SNP.SNS.core.closest.TG, population.SNP.SNS.core.closest.TA)

e <- ggplot(population.SNP.SNS.core.closest.df, aes(x=dist, y=mut, color=group)) +
  geom_line(size = 0.5) + ylim(-0.6,0.6) +
  scale_color_manual(values=c("#E69F00", "#56B4E9", "#009E73","#CC79A7", "#0072B2", "#D55E00")) +
  xlab("Distance from origin (bp)") + ylab("Mutation density (log2 FC)") + ggtitle("Population SNP / SNS core") +
  theme_bw() + theme(aspect.ratio=0.5)

# At SNS stochastic origins

population.SNP.SNS.stochastic.closest.pyr <- population.SNP.SNS.stochastic.closest %>% filter(nchar(ALT) == 1) %>% mutate(dist = (as.numeric(cut(dist, breaks = seq(-10000, 10000, 100)))-100)*100) %>%
  mutate(class = paste0(REF, ">", ALT)) %>%
  mutate(class = gsub("G>A", "C>T", class), class = gsub("G>C", "C>G", class), class = gsub("G>T", "C>A", class), class = gsub("A>G", "T>C", class),class = gsub("A>C", "T>G", class), class = gsub("A>T", "T>A", class)) %>% 
  group_by(dist) %>% summarise(CT.mut = sum(class == "C>T"), CG.mut = sum(class == "C>G") , CA.mut = sum(class == "C>A"), TC.mut = sum(class == "T>C"), TG.mut = sum(class == "T>G") , TA.mut = sum(class == "T>A")) %>% 
  right_join(SNS.stochastic.inter.summary) %>% mutate(CT.mut.cor = CT.mut/(G+C), CG.mut.cor = CG.mut/(G+C), CA.mut.cor = CA.mut/(G+C), TC.mut.cor = TC.mut/(T+A), TG.mut.cor = TG.mut/(T+A), TA.mut.cor = TA.mut/(T+A)) %>% 
  mutate(CT.mut.fold = log2(CT.mut.cor/mean(CT.mut.cor[c(1:20,180:200)])), CG.mut.fold = log2(CG.mut.cor/mean(CG.mut.cor[c(1:20,180:200)])), CA.mut.fold = log2(CA.mut.cor/mean(CA.mut.cor[c(1:20,180:200)])),
         TC.mut.fold = log2(TC.mut.cor/mean(TC.mut.cor[c(1:20,180:200)])), TG.mut.fold = log2(TG.mut.cor/mean(TG.mut.cor[c(1:20,180:200)])), TA.mut.fold = log2(TA.mut.cor/mean(TA.mut.cor[c(1:20,180:200)])))

# Prepare plot
population.SNP.SNS.stochastic.closest.CT <- population.SNP.SNS.stochastic.closest.pyr %>% dplyr::select(dist, mut = CT.mut.fold) %>% mutate(group = "C>T")
population.SNP.SNS.stochastic.closest.CG <- population.SNP.SNS.stochastic.closest.pyr %>% dplyr::select(dist, mut = CG.mut.fold) %>% mutate(group = "C>G")
population.SNP.SNS.stochastic.closest.CA <- population.SNP.SNS.stochastic.closest.pyr %>% dplyr::select(dist, mut = CA.mut.fold) %>% mutate(group = "C>A")
population.SNP.SNS.stochastic.closest.TC <- population.SNP.SNS.stochastic.closest.pyr %>% dplyr::select(dist, mut = TC.mut.fold) %>% mutate(group = "T>C")
population.SNP.SNS.stochastic.closest.TG <- population.SNP.SNS.stochastic.closest.pyr %>% dplyr::select(dist, mut = TG.mut.fold) %>% mutate(group = "T>G")
population.SNP.SNS.stochastic.closest.TA <- population.SNP.SNS.stochastic.closest.pyr %>% dplyr::select(dist, mut = TA.mut.fold) %>% mutate(group = "T>A")
population.SNP.SNS.stochastic.closest.df <- rbind(population.SNP.SNS.stochastic.closest.CT, population.SNP.SNS.stochastic.closest.CG, population.SNP.SNS.stochastic.closest.CA,
                                                  population.SNP.SNS.stochastic.closest.TC, population.SNP.SNS.stochastic.closest.TG, population.SNP.SNS.stochastic.closest.TA)

f <- ggplot(population.SNP.SNS.stochastic.closest.df, aes(x=dist, y=mut, color=group)) +
  geom_line(size = 0.5) + ylim(-0.6,0.6) +
  scale_color_manual(values=c("#E69F00", "#56B4E9", "#009E73","#CC79A7", "#0072B2", "#D55E00")) +
  xlab("Distance from origin (bp)") + ylab("Mutation density (log2 FC)") + ggtitle("Population SNP / SNS stochastic") +
  theme_bw() + theme(aspect.ratio=0.5)

pdf("./Rplot/Fig_1E.pdf", width=8, height=7, useDingbats=FALSE)
ggarrange(d, e, f, ncol = 1, nrow = 3)
dev.off()

############################################################################
############################################################################
###                     Replication strand bias                          ###
############################################################################
############################################################################

# Convert SNP distances in 100 bp bins

population.SNP.SNS.core.closest.100bp <- population.SNP.SNS.core.closest %>% mutate(dist = (as.numeric(cut(dist, breaks = seq(-10000, 10000, 100)))-100)*100)
population.SNP.SNS.stochastic.closest.100bp <- population.SNP.SNS.stochastic.closest %>% mutate(dist = (as.numeric(cut(dist, breaks = seq(-10000, 10000, 100)))-100)*100)
population.SNP.IniSeq.closest.100bp <- population.SNP.IniSeq.closest %>% mutate(dist = (as.numeric(cut(dist, breaks = seq(-10000, 10000, 100)))-100)*100)

############################################################################
# For IniSeq origins

# C mutations

# For C:G > A:T transversion
CG.AT.IniSeq.10kb.bin <- population.SNP.IniSeq.closest.100bp %>% filter(REF == "C" & ALT == "A")
GC.TA.IniSeq.10kb.bin <- population.SNP.IniSeq.closest.100bp %>% filter(REF == "G" & ALT == "T")
# Prepare frequency tables
CG.AT.IniSeq.10kb.freq <- as_data_frame(ftable(CG.AT.IniSeq.10kb.bin$dist))
colnames(CG.AT.IniSeq.10kb.freq) <- c("Distance", "Count.CG.AT")
#
GC.TA.IniSeq.10kb.freq <- as_data_frame(ftable(GC.TA.IniSeq.10kb.bin$dist))
colnames(GC.TA.IniSeq.10kb.freq) <- c("Distance", "Count.GC.TA")
#
CG.AT.IniSeq.assymetry <- CG.AT.IniSeq.10kb.freq %>% right_join(GC.TA.IniSeq.10kb.freq, by = "Distance") %>% mutate(bias = log2(Count.CG.AT/Count.GC.TA)) %>% 
  mutate(class = "C>A", Distance = (as.numeric(Distance)-100)*100) %>% dplyr::dplyr::select(Distance, bias, class)

# For C:G > G:C transversion
CG.GC.IniSeq.10kb.bin <- population.SNP.IniSeq.closest.100bp %>% filter(REF == "C" & ALT == "G")
GC.CG.IniSeq.10kb.bin <- population.SNP.IniSeq.closest.100bp %>% filter(REF == "G" & ALT == "C")
# Prepare frequency tables
CG.GC.IniSeq.10kb.freq <- as_data_frame(ftable(CG.GC.IniSeq.10kb.bin$dist))
colnames(CG.GC.IniSeq.10kb.freq) <- c("Distance", "Count.CG.GC")
#
GC.CG.IniSeq.10kb.freq <- as_data_frame(ftable(GC.CG.IniSeq.10kb.bin$dist))
colnames(GC.CG.IniSeq.10kb.freq) <- c("Distance", "Count.GC.CG")
#
CG.GC.IniSeq.assymetry <- CG.GC.IniSeq.10kb.freq %>% right_join(GC.CG.IniSeq.10kb.freq, by = "Distance") %>% mutate(bias = log2(Count.CG.GC/Count.GC.CG)) %>% 
  mutate(class = "C>G", Distance = (as.numeric(Distance)-100)*100) %>% dplyr::dplyr::select(Distance, bias, class)

# For C:G > T:A transition
CG.TA.IniSeq.10kb.bin <- population.SNP.IniSeq.closest.100bp %>% filter(REF == "C" & ALT == "T")
GC.AT.IniSeq.10kb.bin <- population.SNP.IniSeq.closest.100bp %>% filter(REF == "G" & ALT == "A")
# Prepare frequency tables
CG.TA.IniSeq.10kb.freq <- as_data_frame(ftable(CG.TA.IniSeq.10kb.bin$dist))
colnames(CG.TA.IniSeq.10kb.freq) <- c("Distance", "Count.CG.TA")
#
GC.AT.IniSeq.10kb.freq <- as_data_frame(ftable(GC.AT.IniSeq.10kb.bin$dist))
colnames(GC.AT.IniSeq.10kb.freq) <- c("Distance", "Count.GC.AT")
#
CG.TA.IniSeq.assymetry <- CG.TA.IniSeq.10kb.freq %>% right_join(GC.AT.IniSeq.10kb.freq, by = "Distance") %>% mutate(bias = log2(Count.CG.TA/Count.GC.AT)) %>% 
  mutate(class = "C>T", Distance = (as.numeric(Distance)-100)*100) %>% dplyr::dplyr::select(Distance, bias, class)

# T mutations

# For T:A > A:T transversion
TA.AT.IniSeq.10kb.bin <- population.SNP.IniSeq.closest.100bp %>% filter(REF == "T" & ALT == "A")
AT.TA.IniSeq.10kb.bin <- population.SNP.IniSeq.closest.100bp %>% filter(REF == "A" & ALT == "T")
# Prepare frequency tables
TA.AT.IniSeq.10kb.freq <- as_data_frame(ftable(TA.AT.IniSeq.10kb.bin$dist))
colnames(TA.AT.IniSeq.10kb.freq) <- c("Distance", "Count.TA.AT")
#
AT.TA.IniSeq.10kb.freq <- as_data_frame(ftable(AT.TA.IniSeq.10kb.bin$dist))
colnames(AT.TA.IniSeq.10kb.freq) <- c("Distance", "Count.AT.TA")
#
TA.AT.IniSeq.assymetry <- TA.AT.IniSeq.10kb.freq %>% right_join(AT.TA.IniSeq.10kb.freq, by = "Distance") %>% mutate(bias = log2(Count.TA.AT/Count.AT.TA)) %>% 
  mutate(class = "T>A", Distance = (as.numeric(Distance)-100)*100) %>% dplyr::dplyr::select(Distance, bias, class)

# For T:A > C:G transition
TA.CG.IniSeq.10kb.bin <- population.SNP.IniSeq.closest.100bp %>% filter(REF == "T" & ALT == "C")
AT.GC.IniSeq.10kb.bin <- population.SNP.IniSeq.closest.100bp %>% filter(REF == "A" & ALT == "G")
# Prepare frequency tables
TA.CG.IniSeq.10kb.freq <- as_data_frame(ftable(TA.CG.IniSeq.10kb.bin$dist))
colnames(TA.CG.IniSeq.10kb.freq) <- c("Distance", "Count.TA.CG")
#
AT.GC.IniSeq.10kb.freq <- as_data_frame(ftable(AT.GC.IniSeq.10kb.bin$dist))
colnames(AT.GC.IniSeq.10kb.freq) <- c("Distance", "Count.AT.GC")
#
TA.CG.IniSeq.assymetry <- TA.CG.IniSeq.10kb.freq %>% right_join(AT.GC.IniSeq.10kb.freq, by = "Distance") %>% mutate(bias = log2(Count.TA.CG/Count.AT.GC)) %>% 
  mutate(class = "T>C", Distance = (as.numeric(Distance)-100)*100) %>% dplyr::dplyr::select(Distance, bias, class)

# For T:A > G:C transversion
TA.GC.IniSeq.10kb.bin <- population.SNP.IniSeq.closest.100bp %>% filter(REF == "T" & ALT == "G")
AT.CG.IniSeq.10kb.bin <- population.SNP.IniSeq.closest.100bp %>% filter(REF == "A" & ALT == "C")
# Prepare frequency tables
TA.GC.IniSeq.10kb.freq <- as_data_frame(ftable(TA.GC.IniSeq.10kb.bin$dist))
colnames(TA.GC.IniSeq.10kb.freq) <- c("Distance", "Count.TA.GC")
#
AT.CG.IniSeq.10kb.freq <- as_data_frame(ftable(AT.CG.IniSeq.10kb.bin$dist))
colnames(AT.CG.IniSeq.10kb.freq) <- c("Distance", "Count.AT.CG")
#
TA.GC.IniSeq.assymetry <- TA.GC.IniSeq.10kb.freq %>% right_join(AT.CG.IniSeq.10kb.freq, by = "Distance") %>% mutate(bias = log2(Count.TA.GC/Count.AT.CG)) %>% 
  mutate(class = "T>G", Distance = (as.numeric(Distance)-100)*100) %>% dplyr::dplyr::select(Distance, bias, class)

# Prepare data for plot

C.bias.IniSeq <- rbind(CG.AT.IniSeq.assymetry, CG.GC.IniSeq.assymetry, CG.TA.IniSeq.assymetry)

C.bias.IniSeq.plot <- C.bias.IniSeq %>% ggplot(aes(x=as.numeric(Distance), y=bias, color = class)) +
  geom_point(cex = 0.3) +
  geom_smooth(span = 0.1, alpha = 0.25, size = 1) + ylim(-0.4,0.4) +
  scale_color_manual(values = c("#2EBAED", "#999999", "#DE1C14")) +
  labs(x = "Distance from origin (bp)", y = "Replication bias (log2(FC))", title = "C mutations / IniSeq origins") +
  theme_bw() + theme(aspect.ratio=0.5)

T.bias.IniSeq <- rbind(TA.AT.IniSeq.assymetry, TA.CG.IniSeq.assymetry, TA.GC.IniSeq.assymetry)

T.bias.IniSeq.plot <- T.bias.IniSeq %>% ggplot(aes(x=as.numeric(Distance), y=bias, color = class)) +
  geom_point(cex = 0.3) +
  geom_smooth(span = 0.1, alpha = 0.25, size = 1) + ylim(-0.4,0.4) +
  scale_color_manual(values = c("#0072B2", "#ADCC54", "#D55E00")) +
  labs(x = "Distance from origin (bp)", y = "Replication bias (log2(FC))", title = "T mutations / IniSeq origins") +
  theme_bw() + theme(aspect.ratio=0.5)

pdf("./Rplot/Fig_1G.pdf", width=8, height=7, useDingbats=FALSE)
ggarrange(C.bias.IniSeq.plot, T.bias.IniSeq.plot, ncol = 1, nrow = 2)
dev.off()

# Save analysis for later use

write.csv(CG.AT.IniSeq.assymetry, "./Figure_1/CG.AT.IniSeq.assymetry.csv")
write.csv(CG.GC.IniSeq.assymetry, "./Figure_1/CG.GC.IniSeq.assymetry.csv")
write.csv(TA.CG.IniSeq.assymetry, "./Figure_1/TA.CG.IniSeq.assymetry.csv")
write.csv(TA.GC.IniSeq.assymetry, "./Figure_1/TA.GC.IniSeq.assymetry.csv")

############################################################################
# For SNS core origins

# C mutations

# For C:G > A:T transversion
CG.AT.SNS.core.10kb.bin <- population.SNP.SNS.core.closest.100bp %>% filter(REF == "C" & ALT == "A")
GC.TA.SNS.core.10kb.bin <- population.SNP.SNS.core.closest.100bp %>% filter(REF == "G" & ALT == "T")
# Prepare frequency tables
CG.AT.SNS.core.10kb.freq <- as_data_frame(ftable(CG.AT.SNS.core.10kb.bin$dist))
colnames(CG.AT.SNS.core.10kb.freq) <- c("Distance", "Count.CG.AT")
#
GC.TA.SNS.core.10kb.freq <- as_data_frame(ftable(GC.TA.SNS.core.10kb.bin$dist))
colnames(GC.TA.SNS.core.10kb.freq) <- c("Distance", "Count.GC.TA")
#
CG.AT.SNS.core.assymetry <- CG.AT.SNS.core.10kb.freq %>% right_join(GC.TA.SNS.core.10kb.freq, by = "Distance") %>% mutate(bias = log2(Count.CG.AT/Count.GC.TA)) %>% 
  mutate(class = "C>A", Distance = (as.numeric(Distance)-100)*100) %>% dplyr::dplyr::select(Distance, bias, class)

# For C:G > G:C transversion
CG.GC.SNS.core.10kb.bin <- population.SNP.SNS.core.closest.100bp %>% filter(REF == "C" & ALT == "G")
GC.CG.SNS.core.10kb.bin <- population.SNP.SNS.core.closest.100bp %>% filter(REF == "G" & ALT == "C")
# Prepare frequency tables
CG.GC.SNS.core.10kb.freq <- as_data_frame(ftable(CG.GC.SNS.core.10kb.bin$dist))
colnames(CG.GC.SNS.core.10kb.freq) <- c("Distance", "Count.CG.GC")
#
GC.CG.SNS.core.10kb.freq <- as_data_frame(ftable(GC.CG.SNS.core.10kb.bin$dist))
colnames(GC.CG.SNS.core.10kb.freq) <- c("Distance", "Count.GC.CG")
#
CG.GC.SNS.core.assymetry <- CG.GC.SNS.core.10kb.freq %>% right_join(GC.CG.SNS.core.10kb.freq, by = "Distance") %>% mutate(bias = log2(Count.CG.GC/Count.GC.CG)) %>% 
  mutate(class = "C>G", Distance = (as.numeric(Distance)-100)*100) %>% dplyr::dplyr::select(Distance, bias, class)

# For C:G > T:A transition
CG.TA.SNS.core.10kb.bin <- population.SNP.SNS.core.closest.100bp %>% filter(REF == "C" & ALT == "T")
GC.AT.SNS.core.10kb.bin <- population.SNP.SNS.core.closest.100bp %>% filter(REF == "G" & ALT == "A")
# Prepare frequency tables
CG.TA.SNS.core.10kb.freq <- as_data_frame(ftable(CG.TA.SNS.core.10kb.bin$dist))
colnames(CG.TA.SNS.core.10kb.freq) <- c("Distance", "Count.CG.TA")
#
GC.AT.SNS.core.10kb.freq <- as_data_frame(ftable(GC.AT.SNS.core.10kb.bin$dist))
colnames(GC.AT.SNS.core.10kb.freq) <- c("Distance", "Count.GC.AT")
#
CG.TA.SNS.core.assymetry <- CG.TA.SNS.core.10kb.freq %>% right_join(GC.AT.SNS.core.10kb.freq, by = "Distance") %>% mutate(bias = log2(Count.CG.TA/Count.GC.AT)) %>% 
  mutate(class = "C>T", Distance = (as.numeric(Distance)-100)*100) %>% dplyr::dplyr::select(Distance, bias, class)

# T mutations

# For T:A > A:T transversion
TA.AT.SNS.core.10kb.bin <- population.SNP.SNS.core.closest.100bp %>% filter(REF == "T" & ALT == "A")
AT.TA.SNS.core.10kb.bin <- population.SNP.SNS.core.closest.100bp %>% filter(REF == "A" & ALT == "T")
# Prepare frequency tables
TA.AT.SNS.core.10kb.freq <- as_data_frame(ftable(TA.AT.SNS.core.10kb.bin$dist))
colnames(TA.AT.SNS.core.10kb.freq) <- c("Distance", "Count.TA.AT")
#
AT.TA.SNS.core.10kb.freq <- as_data_frame(ftable(AT.TA.SNS.core.10kb.bin$dist))
colnames(AT.TA.SNS.core.10kb.freq) <- c("Distance", "Count.AT.TA")
#
TA.AT.SNS.core.assymetry <- TA.AT.SNS.core.10kb.freq %>% right_join(AT.TA.SNS.core.10kb.freq, by = "Distance") %>% mutate(bias = log2(Count.TA.AT/Count.AT.TA)) %>% 
  mutate(class = "T>A", Distance = (as.numeric(Distance)-100)*100) %>% dplyr::dplyr::select(Distance, bias, class)

# For T:A > C:G transition
TA.CG.SNS.core.10kb.bin <- population.SNP.SNS.core.closest.100bp %>% filter(REF == "T" & ALT == "C")
AT.GC.SNS.core.10kb.bin <- population.SNP.SNS.core.closest.100bp %>% filter(REF == "A" & ALT == "G")
# Prepare frequency tables
TA.CG.SNS.core.10kb.freq <- as_data_frame(ftable(TA.CG.SNS.core.10kb.bin$dist))
colnames(TA.CG.SNS.core.10kb.freq) <- c("Distance", "Count.TA.CG")
#
AT.GC.SNS.core.10kb.freq <- as_data_frame(ftable(AT.GC.SNS.core.10kb.bin$dist))
colnames(AT.GC.SNS.core.10kb.freq) <- c("Distance", "Count.AT.GC")
#
TA.CG.SNS.core.assymetry <- TA.CG.SNS.core.10kb.freq %>% right_join(AT.GC.SNS.core.10kb.freq, by = "Distance") %>% mutate(bias = log2(Count.TA.CG/Count.AT.GC)) %>% 
  mutate(class = "T>C", Distance = (as.numeric(Distance)-100)*100) %>% dplyr::dplyr::select(Distance, bias, class)

# For T:A > G:C transversion
TA.GC.SNS.core.10kb.bin <- population.SNP.SNS.core.closest.100bp %>% filter(REF == "T" & ALT == "G")
AT.CG.SNS.core.10kb.bin <- population.SNP.SNS.core.closest.100bp %>% filter(REF == "A" & ALT == "C")
# Prepare frequency tables
TA.GC.SNS.core.10kb.freq <- as_data_frame(ftable(TA.GC.SNS.core.10kb.bin$dist))
colnames(TA.GC.SNS.core.10kb.freq) <- c("Distance", "Count.TA.GC")
#
AT.CG.SNS.core.10kb.freq <- as_data_frame(ftable(AT.CG.SNS.core.10kb.bin$dist))
colnames(AT.CG.SNS.core.10kb.freq) <- c("Distance", "Count.AT.CG")
#
TA.GC.SNS.core.assymetry <- TA.GC.SNS.core.10kb.freq %>% right_join(AT.CG.SNS.core.10kb.freq, by = "Distance") %>% mutate(bias = log2(Count.TA.GC/Count.AT.CG)) %>% 
  mutate(class = "T>G", Distance = (as.numeric(Distance)-100)*100) %>% dplyr::dplyr::select(Distance, bias, class)

# Prepare data for plot

C.bias.SNS.core <- rbind(CG.AT.SNS.core.assymetry, CG.GC.SNS.core.assymetry, CG.TA.SNS.core.assymetry)

C.bias.SNS.core.plot <- C.bias.SNS.core %>% ggplot(aes(x=as.numeric(Distance), y=bias, color = class)) +
  geom_point(cex = 0.3) +
  geom_smooth(span = 0.1, alpha = 0.25, size = 1) + ylim(-0.4,0.4) +
  scale_color_manual(values = c("#2EBAED", "#999999", "#DE1C14")) +
  labs(x = "Distance from origin (bp)", y = "Replication bias (log2(FC))", title = "C mutations / SNS core origins") +
  theme_bw() + theme(aspect.ratio=0.5)

T.bias.SNS.core <- rbind(TA.AT.SNS.core.assymetry, TA.CG.SNS.core.assymetry, TA.GC.SNS.core.assymetry)

T.bias.SNS.core.plot <- T.bias.SNS.core %>% ggplot(aes(x=as.numeric(Distance), y=bias, color = class)) +
  geom_point(cex = 0.3) +
  geom_smooth(span = 0.1, alpha = 0.25, size = 1) + ylim(-0.4,0.4) +
  scale_color_manual(values = c("#0072B2", "#ADCC54", "#D55E00")) +
  labs(x = "Distance from origin (bp)", y = "Replication bias (log2(FC))", title = "T mutations / SNS core origins") +
  theme_bw() + theme(aspect.ratio=0.5)

############################################################################
# For SNS stochastic origins

# C mutations

# For C:G > A:T transversion
CG.AT.SNS.stochastic.10kb.bin <- population.SNP.SNS.stochastic.closest.100bp %>% filter(REF == "C" & ALT == "A")
GC.TA.SNS.stochastic.10kb.bin <- population.SNP.SNS.stochastic.closest.100bp %>% filter(REF == "G" & ALT == "T")
# Prepare frequency tables
CG.AT.SNS.stochastic.10kb.freq <- as_data_frame(ftable(CG.AT.SNS.stochastic.10kb.bin$dist))
colnames(CG.AT.SNS.stochastic.10kb.freq) <- c("Distance", "Count.CG.AT")
#
GC.TA.SNS.stochastic.10kb.freq <- as_data_frame(ftable(GC.TA.SNS.stochastic.10kb.bin$dist))
colnames(GC.TA.SNS.stochastic.10kb.freq) <- c("Distance", "Count.GC.TA")
#
CG.AT.SNS.stochastic.assymetry <- CG.AT.SNS.stochastic.10kb.freq %>% right_join(GC.TA.SNS.stochastic.10kb.freq, by = "Distance") %>% mutate(bias = log2(Count.CG.AT/Count.GC.TA)) %>% 
  mutate(class = "C>A", Distance = (as.numeric(Distance)-100)*100) %>% dplyr::dplyr::select(Distance, bias, class)

# For C:G > G:C transversion
CG.GC.SNS.stochastic.10kb.bin <- population.SNP.SNS.stochastic.closest.100bp %>% filter(REF == "C" & ALT == "G")
GC.CG.SNS.stochastic.10kb.bin <- population.SNP.SNS.stochastic.closest.100bp %>% filter(REF == "G" & ALT == "C")
# Prepare frequency tables
CG.GC.SNS.stochastic.10kb.freq <- as_data_frame(ftable(CG.GC.SNS.stochastic.10kb.bin$dist))
colnames(CG.GC.SNS.stochastic.10kb.freq) <- c("Distance", "Count.CG.GC")
#
GC.CG.SNS.stochastic.10kb.freq <- as_data_frame(ftable(GC.CG.SNS.stochastic.10kb.bin$dist))
colnames(GC.CG.SNS.stochastic.10kb.freq) <- c("Distance", "Count.GC.CG")
#
CG.GC.SNS.stochastic.assymetry <- CG.GC.SNS.stochastic.10kb.freq %>% right_join(GC.CG.SNS.stochastic.10kb.freq, by = "Distance") %>% mutate(bias = log2(Count.CG.GC/Count.GC.CG)) %>% 
  mutate(class = "C>G", Distance = (as.numeric(Distance)-100)*100) %>% dplyr::dplyr::select(Distance, bias, class)

# For C:G > T:A transition
CG.TA.SNS.stochastic.10kb.bin <- population.SNP.SNS.stochastic.closest.100bp %>% filter(REF == "C" & ALT == "T")
GC.AT.SNS.stochastic.10kb.bin <- population.SNP.SNS.stochastic.closest.100bp %>% filter(REF == "G" & ALT == "A")
# Prepare frequency tables
CG.TA.SNS.stochastic.10kb.freq <- as_data_frame(ftable(CG.TA.SNS.stochastic.10kb.bin$dist))
colnames(CG.TA.SNS.stochastic.10kb.freq) <- c("Distance", "Count.CG.TA")
#
GC.AT.SNS.stochastic.10kb.freq <- as_data_frame(ftable(GC.AT.SNS.stochastic.10kb.bin$dist))
colnames(GC.AT.SNS.stochastic.10kb.freq) <- c("Distance", "Count.GC.AT")
#
CG.TA.SNS.stochastic.assymetry <- CG.TA.SNS.stochastic.10kb.freq %>% right_join(GC.AT.SNS.stochastic.10kb.freq, by = "Distance") %>% mutate(bias = log2(Count.CG.TA/Count.GC.AT)) %>% 
  mutate(class = "C>T", Distance = (as.numeric(Distance)-100)*100) %>% dplyr::dplyr::select(Distance, bias, class)

# T mutations

# For T:A > A:T transversion
TA.AT.SNS.stochastic.10kb.bin <- population.SNP.SNS.stochastic.closest.100bp %>% filter(REF == "T" & ALT == "A")
AT.TA.SNS.stochastic.10kb.bin <- population.SNP.SNS.stochastic.closest.100bp %>% filter(REF == "A" & ALT == "T")
# Prepare frequency tables
TA.AT.SNS.stochastic.10kb.freq <- as_data_frame(ftable(TA.AT.SNS.stochastic.10kb.bin$dist))
colnames(TA.AT.SNS.stochastic.10kb.freq) <- c("Distance", "Count.TA.AT")
#
AT.TA.SNS.stochastic.10kb.freq <- as_data_frame(ftable(AT.TA.SNS.stochastic.10kb.bin$dist))
colnames(AT.TA.SNS.stochastic.10kb.freq) <- c("Distance", "Count.AT.TA")
#
TA.AT.SNS.stochastic.assymetry <- TA.AT.SNS.stochastic.10kb.freq %>% right_join(AT.TA.SNS.stochastic.10kb.freq, by = "Distance") %>% mutate(bias = log2(Count.TA.AT/Count.AT.TA)) %>% 
  mutate(class = "T>A", Distance = (as.numeric(Distance)-100)*100) %>% dplyr::dplyr::select(Distance, bias, class)

# For T:A > C:G transition
TA.CG.SNS.stochastic.10kb.bin <- population.SNP.SNS.stochastic.closest.100bp %>% filter(REF == "T" & ALT == "C")
AT.GC.SNS.stochastic.10kb.bin <- population.SNP.SNS.stochastic.closest.100bp %>% filter(REF == "A" & ALT == "G")
# Prepare frequency tables
TA.CG.SNS.stochastic.10kb.freq <- as_data_frame(ftable(TA.CG.SNS.stochastic.10kb.bin$dist))
colnames(TA.CG.SNS.stochastic.10kb.freq) <- c("Distance", "Count.TA.CG")
#
AT.GC.SNS.stochastic.10kb.freq <- as_data_frame(ftable(AT.GC.SNS.stochastic.10kb.bin$dist))
colnames(AT.GC.SNS.stochastic.10kb.freq) <- c("Distance", "Count.AT.GC")
#
TA.CG.SNS.stochastic.assymetry <- TA.CG.SNS.stochastic.10kb.freq %>% right_join(AT.GC.SNS.stochastic.10kb.freq, by = "Distance") %>% mutate(bias = log2(Count.TA.CG/Count.AT.GC)) %>% 
  mutate(class = "T>C", Distance = (as.numeric(Distance)-100)*100) %>% dplyr::dplyr::select(Distance, bias, class)

# For T:A > G:C transversion
TA.GC.SNS.stochastic.10kb.bin <- population.SNP.SNS.stochastic.closest.100bp %>% filter(REF == "T" & ALT == "G")
AT.CG.SNS.stochastic.10kb.bin <- population.SNP.SNS.stochastic.closest.100bp %>% filter(REF == "A" & ALT == "C")
# Prepare frequency tables
TA.GC.SNS.stochastic.10kb.freq <- as_data_frame(ftable(TA.GC.SNS.stochastic.10kb.bin$dist))
colnames(TA.GC.SNS.stochastic.10kb.freq) <- c("Distance", "Count.TA.GC")
#
AT.CG.SNS.stochastic.10kb.freq <- as_data_frame(ftable(AT.CG.SNS.stochastic.10kb.bin$dist))
colnames(AT.CG.SNS.stochastic.10kb.freq) <- c("Distance", "Count.AT.CG")
#
TA.GC.SNS.stochastic.assymetry <- TA.GC.SNS.stochastic.10kb.freq %>% right_join(AT.CG.SNS.stochastic.10kb.freq, by = "Distance") %>% mutate(bias = log2(Count.TA.GC/Count.AT.CG)) %>% 
  mutate(class = "T>G", Distance = (as.numeric(Distance)-100)*100) %>% dplyr::dplyr::select(Distance, bias, class)

# Prepare data for plot

C.bias.SNS.stochastic <- rbind(CG.AT.SNS.stochastic.assymetry, CG.GC.SNS.stochastic.assymetry, CG.TA.SNS.stochastic.assymetry)

C.bias.SNS.stochastic.plot <- C.bias.SNS.stochastic %>% ggplot(aes(x=as.numeric(Distance), y=bias, color = class)) +
  geom_point(cex = 0.3) +
  geom_smooth(span = 0.1, alpha = 0.25, size = 1) + ylim(-0.4,0.4) +
  scale_color_manual(values = c("#2EBAED", "#999999", "#DE1C14")) +
  labs(x = "Distance from origin (bp)", y = "Replication bias (log2(FC))", title = "C mutations / SNS stochastic origins") +
  theme_bw() + theme(aspect.ratio=0.5)

T.bias.SNS.stochastic <- rbind(TA.AT.SNS.stochastic.assymetry, TA.CG.SNS.stochastic.assymetry, TA.GC.SNS.stochastic.assymetry)

T.bias.SNS.stochastic.plot <- T.bias.SNS.stochastic %>% ggplot(aes(x=as.numeric(Distance), y=bias, color = class)) +
  geom_point(cex = 0.3) +
  geom_smooth(span = 0.1, alpha = 0.25, size = 1) + ylim(-0.4,0.4) +
  scale_color_manual(values = c("#0072B2", "#ADCC54", "#D55E00")) +
  labs(x = "Distance from origin (bp)", y = "Replication bias (log2(FC))", title = "T mutations / SNS stochastic origins") +
  theme_bw() + theme(aspect.ratio=0.5)

############################################################################
# For IniSeq origins binned by efficiency

# Compute efficiency bins based on
#   33.3%     66.6%   quantiles

quantile(population.SNP.IniSeq.closest.100bp$EFF, c(0.333,0.666))
#    33.3%    66.6% 
# low    medium  high
#   0.8728     0.9297 

# For C:G > A:T transversion

CG.AT.low <- CG.AT.IniSeq.10kb.bin %>% filter(EFF <= 0.8728)
GC.TA.low <- GC.TA.IniSeq.10kb.bin %>% filter(EFF <= 0.8728)
CG.AT.low.freq <- as_data_frame(ftable(CG.AT.low$dist))
colnames(CG.AT.low.freq) <- c("Distance", "Count.CG.AT")
GC.TA.low.freq <- as_data_frame(ftable(GC.TA.low$dist))
colnames(GC.TA.low.freq) <- c("Distance", "Count.GC.TA")
CG.AT.low.assymetry <- CG.AT.low.freq %>% right_join(GC.TA.low.freq, by = "Distance") %>% mutate(bias = log2(Count.CG.AT/Count.GC.TA), EFF = "low", Distance = (as.numeric(Distance)-100)*100) %>% dplyr::dplyr::select(Distance, bias, EFF)

CG.AT.medium <- CG.AT.IniSeq.10kb.bin %>% filter(EFF >= 0.8728 & EFF <= 0.9297)
GC.TA.medium <- GC.TA.IniSeq.10kb.bin %>% filter(EFF >= 0.8728 & EFF <= 0.9297)
CG.AT.medium.freq <- as_data_frame(ftable(CG.AT.medium$dist))
colnames(CG.AT.medium.freq) <- c("Distance", "Count.CG.AT")
GC.TA.medium.freq <- as_data_frame(ftable(GC.TA.medium$dist))
colnames(GC.TA.medium.freq) <- c("Distance", "Count.GC.TA")
CG.AT.medium.assymetry <- CG.AT.medium.freq %>% right_join(GC.TA.medium.freq, by = "Distance") %>% mutate(bias = log2(Count.CG.AT/Count.GC.TA), EFF = "medium", Distance = (as.numeric(Distance)-100)*100) %>% dplyr::dplyr::select(Distance, bias, EFF)

CG.AT.high <- CG.AT.IniSeq.10kb.bin %>% filter(EFF > 0.9297)
GC.TA.high <- GC.TA.IniSeq.10kb.bin %>% filter(EFF > 0.9297)
CG.AT.high.freq <- as_data_frame(ftable(CG.AT.high$dist))
colnames(CG.AT.high.freq) <- c("Distance", "Count.CG.AT")
GC.TA.high.freq <- as_data_frame(ftable(GC.TA.high$dist))
colnames(GC.TA.high.freq) <- c("Distance", "Count.GC.TA")
CG.AT.high.assymetry <- CG.AT.high.freq %>% right_join(GC.TA.high.freq, by = "Distance") %>% mutate(bias = log2(Count.CG.AT/Count.GC.TA), EFF = "high", Distance = (as.numeric(Distance)-100)*100) %>% dplyr::dplyr::select(Distance, bias, EFF)

# For C:G > G:C transversion

CG.GC.low <- CG.GC.IniSeq.10kb.bin %>% filter(EFF <= 0.8728)
GC.CG.low <- GC.CG.IniSeq.10kb.bin %>% filter(EFF <= 0.8728)
CG.GC.low.freq <- as_data_frame(ftable(CG.GC.low$dist))
colnames(CG.GC.low.freq) <- c("Distance", "Count.CG.GC")
GC.CG.low.freq <- as_data_frame(ftable(GC.CG.low$dist))
colnames(GC.CG.low.freq) <- c("Distance", "Count.GC.CG")
CG.GC.low.assymetry <- CG.GC.low.freq %>% right_join(GC.CG.low.freq, by = "Distance") %>% mutate(bias = log2(Count.CG.GC/Count.GC.CG), EFF = "low", Distance = (as.numeric(Distance)-100)*100) %>% dplyr::dplyr::select(Distance, bias, EFF)

CG.GC.medium <- CG.GC.IniSeq.10kb.bin %>% filter(EFF >= 0.8728 & EFF <= 0.9297)
GC.CG.medium <- GC.CG.IniSeq.10kb.bin %>% filter(EFF >= 0.8728 & EFF <= 0.9297)
CG.GC.medium.freq <- as_data_frame(ftable(CG.GC.medium$dist))
colnames(CG.GC.medium.freq) <- c("Distance", "Count.CG.GC")
GC.CG.medium.freq <- as_data_frame(ftable(GC.CG.medium$dist))
colnames(GC.CG.medium.freq) <- c("Distance", "Count.GC.CG")
CG.GC.medium.assymetry <- CG.GC.medium.freq %>% right_join(GC.CG.medium.freq, by = "Distance") %>% mutate(bias = log2(Count.CG.GC/Count.GC.CG), EFF = "medium", Distance = (as.numeric(Distance)-100)*100) %>% dplyr::dplyr::select(Distance, bias, EFF)

CG.GC.high <- CG.GC.IniSeq.10kb.bin %>% filter(EFF > 0.9297)
GC.CG.high <- GC.CG.IniSeq.10kb.bin %>% filter(EFF > 0.9297)
CG.GC.high.freq <- as_data_frame(ftable(CG.GC.high$dist))
colnames(CG.GC.high.freq) <- c("Distance", "Count.CG.GC")
GC.CG.high.freq <- as_data_frame(ftable(GC.CG.high$dist))
colnames(GC.CG.high.freq) <- c("Distance", "Count.GC.CG")
CG.GC.high.assymetry <- CG.GC.high.freq %>% right_join(GC.CG.high.freq, by = "Distance") %>% mutate(bias = log2(Count.CG.GC/Count.GC.CG), EFF = "high", Distance = (as.numeric(Distance)-100)*100) %>% dplyr::dplyr::select(Distance, bias, EFF)

# For C:G > T:A transition

CG.TA.low <- CG.TA.IniSeq.10kb.bin %>% filter(EFF <= 0.8728)
GC.AT.low <- GC.AT.IniSeq.10kb.bin %>% filter(EFF <= 0.8728)
CG.TA.low.freq <- as_data_frame(ftable(CG.TA.low$dist))
colnames(CG.TA.low.freq) <- c("Distance", "Count.CG.TA")
GC.AT.low.freq <- as_data_frame(ftable(GC.AT.low$dist))
colnames(GC.AT.low.freq) <- c("Distance", "Count.GC.AT")
CG.TA.low.assymetry <- CG.TA.low.freq %>% right_join(GC.AT.low.freq, by = "Distance") %>% mutate(bias = log2(Count.CG.TA/Count.GC.AT), EFF = "low", Distance = (as.numeric(Distance)-100)*100) %>% dplyr::dplyr::select(Distance, bias, EFF)

CG.TA.medium <- CG.TA.IniSeq.10kb.bin %>% filter(EFF >= 0.8728 & EFF <= 0.9297)
GC.AT.medium <- GC.AT.IniSeq.10kb.bin %>% filter(EFF >= 0.8728 & EFF <= 0.9297)
CG.TA.medium.freq <- as_data_frame(ftable(CG.TA.medium$dist))
colnames(CG.TA.medium.freq) <- c("Distance", "Count.CG.TA")
GC.AT.medium.freq <- as_data_frame(ftable(GC.AT.medium$dist))
colnames(GC.AT.medium.freq) <- c("Distance", "Count.GC.AT")
CG.TA.medium.assymetry <- CG.TA.medium.freq %>% right_join(GC.AT.medium.freq, by = "Distance") %>% mutate(bias = log2(Count.CG.TA/Count.GC.AT), EFF = "medium", Distance = (as.numeric(Distance)-100)*100) %>% dplyr::dplyr::select(Distance, bias, EFF)

CG.TA.high <- CG.TA.IniSeq.10kb.bin %>% filter(EFF > 0.9297)
GC.AT.high <- GC.AT.IniSeq.10kb.bin %>% filter(EFF > 0.9297)
CG.TA.high.freq <- as_data_frame(ftable(CG.TA.high$dist))
colnames(CG.TA.high.freq) <- c("Distance", "Count.CG.TA")
GC.AT.high.freq <- as_data_frame(ftable(GC.AT.high$dist))
colnames(GC.AT.high.freq) <- c("Distance", "Count.GC.AT")
CG.TA.high.assymetry <- CG.TA.high.freq %>% right_join(GC.AT.high.freq, by = "Distance") %>% mutate(bias = log2(Count.CG.TA/Count.GC.AT), EFF = "high", Distance = (as.numeric(Distance)-100)*100) %>% dplyr::dplyr::select(Distance, bias, EFF)

# For T:A > A:T transversion

TA.AT.low <- TA.AT.IniSeq.10kb.bin %>% filter(EFF <= 0.8728)
AT.TA.low <- AT.TA.IniSeq.10kb.bin %>% filter(EFF <= 0.8728)
TA.AT.low.freq <- as_data_frame(ftable(TA.AT.low$dist))
colnames(TA.AT.low.freq) <- c("Distance", "Count.TA.AT")
AT.TA.low.freq <- as_data_frame(ftable(AT.TA.low$dist))
colnames(AT.TA.low.freq) <- c("Distance", "Count.AT.TA")
TA.AT.low.assymetry <- TA.AT.low.freq %>% right_join(AT.TA.low.freq, by = "Distance") %>% mutate(bias = log2(Count.TA.AT/Count.AT.TA), EFF = "low", Distance = (as.numeric(Distance)-100)*100) %>% dplyr::dplyr::select(Distance, bias, EFF)

TA.AT.medium <- TA.AT.IniSeq.10kb.bin %>% filter(EFF >= 0.8728 & EFF <= 0.9297)
AT.TA.medium <- AT.TA.IniSeq.10kb.bin %>% filter(EFF >= 0.8728 & EFF <= 0.9297)
TA.AT.medium.freq <- as_data_frame(ftable(TA.AT.medium$dist))
colnames(TA.AT.medium.freq) <- c("Distance", "Count.TA.AT")
AT.TA.medium.freq <- as_data_frame(ftable(AT.TA.medium$dist))
colnames(AT.TA.medium.freq) <- c("Distance", "Count.AT.TA")
TA.AT.medium.assymetry <- TA.AT.medium.freq %>% right_join(AT.TA.medium.freq, by = "Distance") %>% mutate(bias = log2(Count.TA.AT/Count.AT.TA), EFF = "medium", Distance = (as.numeric(Distance)-100)*100) %>% dplyr::dplyr::select(Distance, bias, EFF)

TA.AT.high <- TA.AT.IniSeq.10kb.bin %>% filter(EFF > 0.9297)
AT.TA.high <- AT.TA.IniSeq.10kb.bin %>% filter(EFF > 0.9297)
TA.AT.high.freq <- as_data_frame(ftable(TA.AT.high$dist))
colnames(TA.AT.high.freq) <- c("Distance", "Count.TA.AT")
AT.TA.high.freq <- as_data_frame(ftable(AT.TA.high$dist))
colnames(AT.TA.high.freq) <- c("Distance", "Count.AT.TA")
TA.AT.high.assymetry <- TA.AT.high.freq %>% right_join(AT.TA.high.freq, by = "Distance") %>% mutate(bias = log2(Count.TA.AT/Count.AT.TA), EFF = "high", Distance = (as.numeric(Distance)-100)*100) %>% dplyr::dplyr::select(Distance, bias, EFF)

# For T:A > C:G transition

TA.CG.low <- TA.CG.IniSeq.10kb.bin %>% filter(EFF <= 0.8728)
AT.GC.low <- AT.GC.IniSeq.10kb.bin %>% filter(EFF <= 0.8728)
TA.CG.low.freq <- as_data_frame(ftable(TA.CG.low$dist))
colnames(TA.CG.low.freq) <- c("Distance", "Count.TA.CG")
AT.GC.low.freq <- as_data_frame(ftable(AT.GC.low$dist))
colnames(AT.GC.low.freq) <- c("Distance", "Count.AT.GC")
TA.CG.low.assymetry <- TA.CG.low.freq %>% right_join(AT.GC.low.freq, by = "Distance") %>% mutate(bias = log2(Count.TA.CG/Count.AT.GC), EFF = "low", Distance = (as.numeric(Distance)-100)*100) %>% dplyr::dplyr::select(Distance, bias, EFF)

TA.CG.medium <- TA.CG.IniSeq.10kb.bin %>% filter(EFF >= 0.8728 & EFF <= 0.9297)
AT.GC.medium <- AT.GC.IniSeq.10kb.bin %>% filter(EFF >= 0.8728 & EFF <= 0.9297)
TA.CG.medium.freq <- as_data_frame(ftable(TA.CG.medium$dist))
colnames(TA.CG.medium.freq) <- c("Distance", "Count.TA.CG")
AT.GC.medium.freq <- as_data_frame(ftable(AT.GC.medium$dist))
colnames(AT.GC.medium.freq) <- c("Distance", "Count.AT.GC")
TA.CG.medium.assymetry <- TA.CG.medium.freq %>% right_join(AT.GC.medium.freq, by = "Distance") %>% mutate(bias = log2(Count.TA.CG/Count.AT.GC), EFF = "medium", Distance = (as.numeric(Distance)-100)*100) %>% dplyr::dplyr::select(Distance, bias, EFF)

TA.CG.high <- TA.CG.IniSeq.10kb.bin %>% filter(EFF > 0.9297)
AT.GC.high <- AT.GC.IniSeq.10kb.bin %>% filter(EFF > 0.9297)
TA.CG.high.freq <- as_data_frame(ftable(TA.CG.high$dist))
colnames(TA.CG.high.freq) <- c("Distance", "Count.TA.CG")
AT.GC.high.freq <- as_data_frame(ftable(AT.GC.high$dist))
colnames(AT.GC.high.freq) <- c("Distance", "Count.AT.GC")
TA.CG.high.assymetry <- TA.CG.high.freq %>% right_join(AT.GC.high.freq, by = "Distance") %>% mutate(bias = log2(Count.TA.CG/Count.AT.GC), EFF = "high", Distance = (as.numeric(Distance)-100)*100) %>% dplyr::dplyr::select(Distance, bias, EFF)

# For T:A > G:C transversion

TA.GC.low <- TA.GC.IniSeq.10kb.bin %>% filter(EFF <= 0.8728)
AT.CG.low <- AT.CG.IniSeq.10kb.bin %>% filter(EFF <= 0.8728)
TA.GC.low.freq <- as_data_frame(ftable(TA.GC.low$dist))
colnames(TA.GC.low.freq) <- c("Distance", "Count.TA.GC")
AT.CG.low.freq <- as_data_frame(ftable(AT.CG.low$dist))
colnames(AT.CG.low.freq) <- c("Distance", "Count.AT.CG")
TA.GC.low.assymetry <- TA.GC.low.freq %>% right_join(AT.CG.low.freq, by = "Distance") %>% mutate(bias = log2(Count.TA.GC/Count.AT.CG), EFF = "low", Distance = (as.numeric(Distance)-100)*100) %>% dplyr::dplyr::select(Distance, bias, EFF)

TA.GC.medium <- TA.GC.IniSeq.10kb.bin %>% filter(EFF >= 0.8728 & EFF <= 0.9297)
AT.CG.medium <- AT.CG.IniSeq.10kb.bin %>% filter(EFF >= 0.8728 & EFF <= 0.9297)
TA.GC.medium.freq <- as_data_frame(ftable(TA.GC.medium$dist))
colnames(TA.GC.medium.freq) <- c("Distance", "Count.TA.GC")
AT.CG.medium.freq <- as_data_frame(ftable(AT.CG.medium$dist))
colnames(AT.CG.medium.freq) <- c("Distance", "Count.AT.CG")
TA.GC.medium.assymetry <- TA.GC.medium.freq %>% right_join(AT.CG.medium.freq, by = "Distance") %>% mutate(bias = log2(Count.TA.GC/Count.AT.CG), EFF = "medium", Distance = (as.numeric(Distance)-100)*100) %>% dplyr::dplyr::select(Distance, bias, EFF)

TA.GC.high <- TA.GC.IniSeq.10kb.bin %>% filter(EFF > 0.9297)
AT.CG.high <- AT.CG.IniSeq.10kb.bin %>% filter(EFF > 0.9297)
TA.GC.high.freq <- as_data_frame(ftable(TA.GC.high$dist))
colnames(TA.GC.high.freq) <- c("Distance", "Count.TA.GC")
AT.CG.high.freq <- as_data_frame(ftable(AT.CG.high$dist))
colnames(AT.CG.high.freq) <- c("Distance", "Count.AT.CG")
TA.GC.high.assymetry <- TA.GC.high.freq %>% right_join(AT.CG.high.freq, by = "Distance") %>% mutate(bias = log2(Count.TA.GC/Count.AT.CG), EFF = "high", Distance = (as.numeric(Distance)-100)*100) %>% dplyr::dplyr::select(Distance, bias, EFF)

############################################################################
############################################################################
###       Replication strand bias quantification and statistics          ###
############################################################################
############################################################################

# Define standerd error to the mean function
stderr <- function(x, na.rm=FALSE) {
  if (na.rm) x <- na.omit(x)
  sqrt(var(x)/length(x))
}

# Replication strand bias for IniSeq and SNS origins

# IniSeq origins
CG.AT.IniSeq.bias <- CG.AT.IniSeq.assymetry %>% mutate(strand.bias = case_when(Distance < 0 ~ bias, TRUE ~ -bias))
CG.AT.IniSeq.bias.df <- cbind.data.frame(class = "C>A", mean = mean(CG.AT.IniSeq.bias$strand.bias), sem = stderr(CG.AT.IniSeq.bias$strand.bias), ori = "IniSeq")
#  
CG.GC.IniSeq.bias <- CG.GC.IniSeq.assymetry %>% mutate(strand.bias = case_when(Distance < 0 ~ bias, TRUE ~ -bias))
CG.GC.IniSeq.bias.df <- cbind.data.frame(class = "C>G", mean = mean(CG.GC.IniSeq.bias$strand.bias), sem = stderr(CG.GC.IniSeq.bias$strand.bias), ori = "IniSeq")
#  
CG.TA.IniSeq.bias <- CG.TA.IniSeq.assymetry %>% mutate(strand.bias = case_when(Distance < 0 ~ bias, TRUE ~ -bias))
CG.TA.IniSeq.bias.df <- cbind.data.frame(class = "C>T", mean = mean(CG.TA.IniSeq.bias$strand.bias), sem = stderr(CG.TA.IniSeq.bias$strand.bias), ori = "IniSeq")
#  
TA.AT.IniSeq.bias <- TA.AT.IniSeq.assymetry %>% mutate(strand.bias = case_when(Distance < 0 ~ bias, TRUE ~ -bias))
TA.AT.IniSeq.bias.df <- cbind.data.frame(class = "T>A", mean = mean(TA.AT.IniSeq.bias$strand.bias), sem = stderr(TA.AT.IniSeq.bias$strand.bias), ori = "IniSeq")
#  
TA.CG.IniSeq.bias <- TA.CG.IniSeq.assymetry %>% mutate(strand.bias = case_when(Distance < 0 ~ bias, TRUE ~ -bias))
TA.CG.IniSeq.bias.df <- cbind.data.frame(class = "T>C", mean = mean(TA.CG.IniSeq.bias$strand.bias), sem = stderr(TA.CG.IniSeq.bias$strand.bias), ori = "IniSeq")
#  
TA.GC.IniSeq.bias <- TA.GC.IniSeq.assymetry %>% mutate(strand.bias = case_when(Distance < 0 ~ bias, TRUE ~ -bias))
TA.GC.IniSeq.bias.df <- cbind.data.frame(class = "T>G", mean = mean(TA.GC.IniSeq.bias$strand.bias), sem = stderr(TA.GC.IniSeq.bias$strand.bias), ori = "IniSeq")
# Bind
IniSeq.bias <- rbind(CG.AT.IniSeq.bias.df, CG.GC.IniSeq.bias.df, CG.TA.IniSeq.bias.df, TA.AT.IniSeq.bias.df, TA.CG.IniSeq.bias.df, TA.GC.IniSeq.bias.df)

# SNS core origins
CG.AT.SNS.core.bias <- CG.AT.SNS.core.assymetry %>% mutate(strand.bias = case_when(Distance < 0 ~ bias, TRUE ~ -bias))
CG.AT.SNS.core.bias.df <- cbind.data.frame(class = "C>A", mean = mean(CG.AT.SNS.core.bias$strand.bias), sem = stderr(CG.AT.SNS.core.bias$strand.bias), ori = "SNS core")
#  
CG.GC.SNS.core.bias <- CG.GC.SNS.core.assymetry %>% mutate(strand.bias = case_when(Distance < 0 ~ bias, TRUE ~ -bias))
CG.GC.SNS.core.bias.df <- cbind.data.frame(class = "C>G", mean = mean(CG.GC.SNS.core.bias$strand.bias), sem = stderr(CG.GC.SNS.core.bias$strand.bias), ori = "SNS core")
#  
CG.TA.SNS.core.bias <- CG.TA.SNS.core.assymetry %>% mutate(strand.bias = case_when(Distance < 0 ~ bias, TRUE ~ -bias))
CG.TA.SNS.core.bias.df <- cbind.data.frame(class = "C>T", mean = mean(CG.TA.SNS.core.bias$strand.bias), sem = stderr(CG.TA.SNS.core.bias$strand.bias), ori = "SNS core")
#  
TA.AT.SNS.core.bias <- TA.AT.SNS.core.assymetry %>% mutate(strand.bias = case_when(Distance < 0 ~ bias, TRUE ~ -bias))
TA.AT.SNS.core.bias.df <- cbind.data.frame(class = "T>A", mean = mean(TA.AT.SNS.core.bias$strand.bias), sem = stderr(TA.AT.SNS.core.bias$strand.bias), ori = "SNS core")
#  
TA.CG.SNS.core.bias <- TA.CG.SNS.core.assymetry %>% mutate(strand.bias = case_when(Distance < 0 ~ bias, TRUE ~ -bias))
TA.CG.SNS.core.bias.df <- cbind.data.frame(class = "T>C", mean = mean(TA.CG.SNS.core.bias$strand.bias), sem = stderr(TA.CG.SNS.core.bias$strand.bias), ori = "SNS core")
#  
TA.GC.SNS.core.bias <- TA.GC.SNS.core.assymetry %>% mutate(strand.bias = case_when(Distance < 0 ~ bias, TRUE ~ -bias))
TA.GC.SNS.core.bias.df <- cbind.data.frame(class = "T>G", mean = mean(TA.GC.SNS.core.bias$strand.bias), sem = stderr(TA.GC.SNS.core.bias$strand.bias), ori = "SNS core")
# Bind
SNS.core.bias <- rbind(CG.AT.SNS.core.bias.df, CG.GC.SNS.core.bias.df, CG.TA.SNS.core.bias.df, TA.AT.SNS.core.bias.df, TA.CG.SNS.core.bias.df, TA.GC.SNS.core.bias.df)

# SNS stochastic origins
CG.AT.SNS.stochastic.bias <- CG.AT.SNS.stochastic.assymetry %>% mutate(strand.bias = case_when(Distance < 0 ~ bias, TRUE ~ -bias))
CG.AT.SNS.stochastic.bias.df <- cbind.data.frame(class = "C>A", mean = mean(CG.AT.SNS.stochastic.bias$strand.bias), sem = stderr(CG.AT.SNS.stochastic.bias$strand.bias), ori = "SNS stochastic")
#  
CG.GC.SNS.stochastic.bias <- CG.GC.SNS.stochastic.assymetry %>% mutate(strand.bias = case_when(Distance < 0 ~ bias, TRUE ~ -bias))
CG.GC.SNS.stochastic.bias.df <- cbind.data.frame(class = "C>G", mean = mean(CG.GC.SNS.stochastic.bias$strand.bias), sem = stderr(CG.GC.SNS.stochastic.bias$strand.bias), ori = "SNS stochastic")
#  
CG.TA.SNS.stochastic.bias <- CG.TA.SNS.stochastic.assymetry %>% mutate(strand.bias = case_when(Distance < 0 ~ bias, TRUE ~ -bias))
CG.TA.SNS.stochastic.bias.df <- cbind.data.frame(class = "C>T", mean = mean(CG.TA.SNS.stochastic.bias$strand.bias), sem = stderr(CG.TA.SNS.stochastic.bias$strand.bias), ori = "SNS stochastic")
#  
TA.AT.SNS.stochastic.bias <- TA.AT.SNS.stochastic.assymetry %>% mutate(strand.bias = case_when(Distance < 0 ~ bias, TRUE ~ -bias))
TA.AT.SNS.stochastic.bias.df <- cbind.data.frame(class = "T>A", mean = mean(TA.AT.SNS.stochastic.bias$strand.bias), sem = stderr(TA.AT.SNS.stochastic.bias$strand.bias), ori = "SNS stochastic")
#  
TA.CG.SNS.stochastic.bias <- TA.CG.SNS.stochastic.assymetry %>% mutate(strand.bias = case_when(Distance < 0 ~ bias, TRUE ~ -bias))
TA.CG.SNS.stochastic.bias.df <- cbind.data.frame(class = "T>C", mean = mean(TA.CG.SNS.stochastic.bias$strand.bias), sem = stderr(TA.CG.SNS.stochastic.bias$strand.bias), ori = "SNS stochastic")
#  
TA.GC.SNS.stochastic.bias <- TA.GC.SNS.stochastic.assymetry %>% mutate(strand.bias = case_when(Distance < 0 ~ bias, TRUE ~ -bias))
TA.GC.SNS.stochastic.bias.df <- cbind.data.frame(class = "T>G", mean = mean(TA.GC.SNS.stochastic.bias$strand.bias), sem = stderr(TA.GC.SNS.stochastic.bias$strand.bias), ori = "SNS stochastic")
# Bind
SNS.stochastic.bias <- rbind(CG.AT.SNS.stochastic.bias.df, CG.GC.SNS.stochastic.bias.df, CG.TA.SNS.stochastic.bias.df, TA.AT.SNS.stochastic.bias.df, TA.CG.SNS.stochastic.bias.df, TA.GC.SNS.stochastic.bias.df)

# Prepare plot
bias.plot <- rbind(IniSeq.bias, SNS.core.bias, SNS.stochastic.bias)

Origin.bias.plot <- bias.plot %>% mutate(Origins = fct_relevel(ori, "SNS stochastic", "SNS core", "IniSeq")) %>% 
  ggplot(aes(x=class, y=mean, fill=Origins)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem), width=.2,
                position=position_dodge(.9)) + 
  scale_fill_manual(values=c("#CC79A7", "#0072B2", "#FC4E07")) + ylim(-0.15, 0.28) +
  ylab("Replication strand bias\nlog2(leading/lagging)") + labs(fill = "Origins") +
  theme_bw() + theme(aspect.ratio=1, axis.title=element_text(size=12), axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1), axis.title.x=element_blank(), axis.text.y=element_text(size=12))

# Replication strand bias for IniSeq origins binned by efficiency

# Low efficiency
CG.AT.IniSeq.low.bias <- CG.AT.low.assymetry %>% mutate(strand.bias = case_when(Distance < 0 ~ bias, TRUE ~ -bias))
CG.AT.IniSeq.low.bias.df <- cbind.data.frame(class = "C>A", mean = mean(CG.AT.IniSeq.low.bias$strand.bias), sem = stderr(CG.AT.IniSeq.low.bias$strand.bias), eff = "low")
#
CG.GC.IniSeq.low.bias <- CG.GC.low.assymetry %>% mutate(strand.bias = case_when(Distance < 0 ~ bias, TRUE ~ -bias))
CG.GC.IniSeq.low.bias.df <- cbind.data.frame(class = "C>G", mean = mean(CG.GC.IniSeq.low.bias$strand.bias), sem = stderr(CG.GC.IniSeq.low.bias$strand.bias), eff = "low")
#
CG.TA.IniSeq.low.bias <- CG.TA.low.assymetry %>% mutate(strand.bias = case_when(Distance < 0 ~ bias, TRUE ~ -bias))
CG.TA.IniSeq.low.bias.df <- cbind.data.frame(class = "C>T", mean = mean(CG.TA.IniSeq.low.bias$strand.bias), sem = stderr(CG.TA.IniSeq.low.bias$strand.bias), eff = "low")
#
TA.AT.IniSeq.low.bias <- TA.AT.low.assymetry %>% mutate(strand.bias = case_when(Distance < 0 ~ bias, TRUE ~ -bias))
TA.AT.IniSeq.low.bias.df <- cbind.data.frame(class = "T>A", mean = mean(TA.AT.IniSeq.low.bias$strand.bias), sem = stderr(TA.AT.IniSeq.low.bias$strand.bias), eff = "low")
#
TA.CG.IniSeq.low.bias <- TA.CG.low.assymetry %>% mutate(strand.bias = case_when(Distance < 0 ~ bias, TRUE ~ -bias))
TA.CG.IniSeq.low.bias.df <- cbind.data.frame(class = "T>C", mean = mean(TA.CG.IniSeq.low.bias$strand.bias), sem = stderr(TA.CG.IniSeq.low.bias$strand.bias), eff = "low")
#
TA.GC.IniSeq.low.bias <- TA.GC.low.assymetry %>% mutate(strand.bias = case_when(Distance < 0 ~ bias, TRUE ~ -bias))
TA.GC.IniSeq.low.bias.df <- cbind.data.frame(class = "T>G", mean = mean(TA.GC.IniSeq.low.bias$strand.bias), sem = stderr(TA.GC.IniSeq.low.bias$strand.bias), eff = "low")
# Bind
IniSeq.low.bias <- rbind(CG.AT.IniSeq.low.bias.df, CG.GC.IniSeq.low.bias.df, CG.TA.IniSeq.low.bias.df, TA.AT.IniSeq.low.bias.df, TA.CG.IniSeq.low.bias.df, TA.GC.IniSeq.low.bias.df)
# Medium efficiency
CG.AT.IniSeq.medium.bias <- CG.AT.medium.assymetry %>% mutate(strand.bias = case_when(Distance < 0 ~ bias, TRUE ~ -bias))
CG.AT.IniSeq.medium.bias.df <- cbind.data.frame(class = "C>A", mean = mean(CG.AT.IniSeq.medium.bias$strand.bias), sem = stderr(CG.AT.IniSeq.medium.bias$strand.bias), eff = "medium")
#
CG.GC.IniSeq.medium.bias <- CG.GC.medium.assymetry %>% mutate(strand.bias = case_when(Distance < 0 ~ bias, TRUE ~ -bias))
CG.GC.IniSeq.medium.bias.df <- cbind.data.frame(class = "C>G", mean = mean(CG.GC.IniSeq.medium.bias$strand.bias), sem = stderr(CG.GC.IniSeq.medium.bias$strand.bias), eff = "medium")
#
CG.TA.IniSeq.medium.bias <- CG.TA.medium.assymetry %>% mutate(strand.bias = case_when(Distance < 0 ~ bias, TRUE ~ -bias))
CG.TA.IniSeq.medium.bias.df <- cbind.data.frame(class = "C>T", mean = mean(CG.TA.IniSeq.medium.bias$strand.bias), sem = stderr(CG.TA.IniSeq.medium.bias$strand.bias), eff = "medium")
#
TA.AT.IniSeq.medium.bias <- TA.AT.medium.assymetry %>% mutate(strand.bias = case_when(Distance < 0 ~ bias, TRUE ~ -bias))
TA.AT.IniSeq.medium.bias.df <- cbind.data.frame(class = "T>A", mean = mean(TA.AT.IniSeq.medium.bias$strand.bias), sem = stderr(TA.AT.IniSeq.medium.bias$strand.bias), eff = "medium")
#
TA.CG.IniSeq.medium.bias <- TA.CG.medium.assymetry %>% mutate(strand.bias = case_when(Distance < 0 ~ bias, TRUE ~ -bias))
TA.CG.IniSeq.medium.bias.df <- cbind.data.frame(class = "T>C", mean = mean(TA.CG.IniSeq.medium.bias$strand.bias), sem = stderr(TA.CG.IniSeq.medium.bias$strand.bias), eff = "medium")
#
TA.GC.IniSeq.medium.bias <- TA.GC.medium.assymetry %>% mutate(strand.bias = case_when(Distance < 0 ~ bias, TRUE ~ -bias))
TA.GC.IniSeq.medium.bias.df <- cbind.data.frame(class = "T>G", mean = mean(TA.GC.IniSeq.medium.bias$strand.bias), sem = stderr(TA.GC.IniSeq.medium.bias$strand.bias), eff = "medium")
# Bind
IniSeq.medium.bias <- rbind(CG.AT.IniSeq.medium.bias.df, CG.GC.IniSeq.medium.bias.df, CG.TA.IniSeq.medium.bias.df, TA.AT.IniSeq.medium.bias.df, TA.CG.IniSeq.medium.bias.df, TA.GC.IniSeq.medium.bias.df)
# High efficiency
CG.AT.IniSeq.high.bias <- CG.AT.high.assymetry %>% mutate(strand.bias = case_when(Distance < 0 ~ bias, TRUE ~ -bias))
CG.AT.IniSeq.high.bias.df <- cbind.data.frame(class = "C>A", mean = mean(CG.AT.IniSeq.high.bias$strand.bias), sem = stderr(CG.AT.IniSeq.high.bias$strand.bias), eff = "high")
#
CG.GC.IniSeq.high.bias <- CG.GC.high.assymetry %>% mutate(strand.bias = case_when(Distance < 0 ~ bias, TRUE ~ -bias))
CG.GC.IniSeq.high.bias.df <- cbind.data.frame(class = "C>G", mean = mean(CG.GC.IniSeq.high.bias$strand.bias), sem = stderr(CG.GC.IniSeq.high.bias$strand.bias), eff = "high")
#
CG.TA.IniSeq.high.bias <- CG.TA.high.assymetry %>% mutate(strand.bias = case_when(Distance < 0 ~ bias, TRUE ~ -bias))
CG.TA.IniSeq.high.bias.df <- cbind.data.frame(class = "C>T", mean = mean(CG.TA.IniSeq.high.bias$strand.bias), sem = stderr(CG.TA.IniSeq.high.bias$strand.bias), eff = "high")
#
TA.AT.IniSeq.high.bias <- TA.AT.high.assymetry %>% mutate(strand.bias = case_when(Distance < 0 ~ bias, TRUE ~ -bias))
TA.AT.IniSeq.high.bias.df <- cbind.data.frame(class = "T>A", mean = mean(TA.AT.IniSeq.high.bias$strand.bias), sem = stderr(TA.AT.IniSeq.high.bias$strand.bias), eff = "high")
#
TA.CG.IniSeq.high.bias <- TA.CG.high.assymetry %>% mutate(strand.bias = case_when(Distance < 0 ~ bias, TRUE ~ -bias))
TA.CG.IniSeq.high.bias.df <- cbind.data.frame(class = "T>C", mean = mean(TA.CG.IniSeq.high.bias$strand.bias), sem = stderr(TA.CG.IniSeq.high.bias$strand.bias), eff = "high")
#
TA.GC.IniSeq.high.bias <- TA.GC.high.assymetry %>% mutate(strand.bias = case_when(Distance < 0 ~ bias, TRUE ~ -bias))
TA.GC.IniSeq.high.bias.df <- cbind.data.frame(class = "T>G", mean = mean(TA.GC.IniSeq.high.bias$strand.bias), sem = stderr(TA.GC.IniSeq.high.bias$strand.bias), eff = "high")
# Bind
IniSeq.high.bias <- rbind(CG.AT.IniSeq.high.bias.df, CG.GC.IniSeq.high.bias.df, CG.TA.IniSeq.high.bias.df, TA.AT.IniSeq.high.bias.df, TA.CG.IniSeq.high.bias.df, TA.GC.IniSeq.high.bias.df)
# Bind all efficiency
IniSeq.eff.bias <- rbind(IniSeq.low.bias, IniSeq.medium.bias, IniSeq.high.bias)

IniSeq.eff.bias.plot <- IniSeq.eff.bias %>% mutate(Origins = fct_relevel(eff, "low", "medium", "high")) %>% 
  ggplot(aes(x=class, y=mean, fill=Origins)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem), width=.2,
                position=position_dodge(.9)) + ylim(-0.22,0.39) +
  scale_fill_manual(values=c("#CC79A7", "#0072B2", "#FC4E07")) +
  ylab("Replication strand bias\nlog2(leading/lagging)") + labs(fill = "Efficiency") +
  theme_bw() + theme(aspect.ratio=1, axis.title=element_text(size=12), axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1), axis.title.x=element_blank(), axis.text.y=element_text(size=12))

# Save plots

pdf("./Rplot/Fig_1F_1G.pdf", width=12, height=5, useDingbats=FALSE)
ggarrange(Origin.bias.plot, IniSeq.eff.bias.plot, ncol = 2, nrow = 1)
dev.off()

# Compute stats

ks.test(CG.AT.IniSeq.assymetry[which(CG.AT.IniSeq.assymetry$Distance < 0),]$bias, CG.AT.IniSeq.assymetry[which(CG.AT.IniSeq.assymetry$Distance > 0),]$bias) # p-value = 1.11e-15
ks.test(CG.GC.IniSeq.assymetry[which(CG.GC.IniSeq.assymetry$Distance < 0),]$bias, CG.GC.IniSeq.assymetry[which(CG.GC.IniSeq.assymetry$Distance > 0),]$bias) # p-value < 2.2e-16
ks.test(CG.TA.IniSeq.assymetry[which(CG.TA.IniSeq.assymetry$Distance < 0),]$bias, CG.TA.IniSeq.assymetry[which(CG.TA.IniSeq.assymetry$Distance > 0),]$bias) # p-value = 0.002202
ks.test(TA.AT.IniSeq.assymetry[which(TA.AT.IniSeq.assymetry$Distance < 0),]$bias, TA.AT.IniSeq.assymetry[which(TA.AT.IniSeq.assymetry$Distance > 0),]$bias) # p-value = 0.007142
ks.test(TA.CG.IniSeq.assymetry[which(TA.CG.IniSeq.assymetry$Distance < 0),]$bias, TA.CG.IniSeq.assymetry[which(TA.CG.IniSeq.assymetry$Distance > 0),]$bias) # p-value = 1.11e-15
ks.test(TA.GC.IniSeq.assymetry[which(TA.GC.IniSeq.assymetry$Distance < 0),]$bias, TA.GC.IniSeq.assymetry[which(TA.GC.IniSeq.assymetry$Distance > 0),]$bias) # p-value < 2.2e-16
#
ks.test(CG.AT.SNS.core.assymetry[which(CG.AT.SNS.core.assymetry$Distance < 0),]$bias, CG.AT.SNS.core.assymetry[which(CG.AT.SNS.core.assymetry$Distance > 0),]$bias) # p-value = 0.2234
ks.test(CG.GC.SNS.core.assymetry[which(CG.GC.SNS.core.assymetry$Distance < 0),]$bias, CG.GC.SNS.core.assymetry[which(CG.GC.SNS.core.assymetry$Distance > 0),]$bias) # p-value = 0.4221
ks.test(CG.TA.SNS.core.assymetry[which(CG.TA.SNS.core.assymetry$Distance < 0),]$bias, CG.TA.SNS.core.assymetry[which(CG.TA.SNS.core.assymetry$Distance > 0),]$bias) # p-value = 0.06057
ks.test(TA.AT.SNS.core.assymetry[which(TA.AT.SNS.core.assymetry$Distance < 0),]$bias, TA.AT.SNS.core.assymetry[which(TA.AT.SNS.core.assymetry$Distance > 0),]$bias) # p-value = 0.09153
ks.test(TA.CG.SNS.core.assymetry[which(TA.CG.SNS.core.assymetry$Distance < 0),]$bias, TA.CG.SNS.core.assymetry[which(TA.CG.SNS.core.assymetry$Distance > 0),]$bias) # p-value = 0.3314
ks.test(TA.GC.SNS.core.assymetry[which(TA.GC.SNS.core.assymetry$Distance < 0),]$bias, TA.GC.SNS.core.assymetry[which(TA.GC.SNS.core.assymetry$Distance > 0),]$bias) # p-value = 0.3098
#
ks.test(CG.AT.SNS.stochastic.assymetry[which(CG.AT.SNS.stochastic.assymetry$Distance < 0),]$bias, CG.AT.SNS.stochastic.assymetry[which(CG.AT.SNS.stochastic.assymetry$Distance > 0),]$bias) # p-value = 0.2247
ks.test(CG.GC.SNS.stochastic.assymetry[which(CG.GC.SNS.stochastic.assymetry$Distance < 0),]$bias, CG.GC.SNS.stochastic.assymetry[which(CG.GC.SNS.stochastic.assymetry$Distance > 0),]$bias) # p-value = 0.2397
ks.test(CG.TA.SNS.stochastic.assymetry[which(CG.TA.SNS.stochastic.assymetry$Distance < 0),]$bias, CG.TA.SNS.stochastic.assymetry[which(CG.TA.SNS.stochastic.assymetry$Distance > 0),]$bias) # p-value = 0.334
ks.test(TA.AT.SNS.stochastic.assymetry[which(TA.AT.SNS.stochastic.assymetry$Distance < 0),]$bias, TA.AT.SNS.stochastic.assymetry[which(TA.AT.SNS.stochastic.assymetry$Distance > 0),]$bias) # p-value = 0.04331
ks.test(TA.CG.SNS.stochastic.assymetry[which(TA.CG.SNS.stochastic.assymetry$Distance < 0),]$bias, TA.CG.SNS.stochastic.assymetry[which(TA.CG.SNS.stochastic.assymetry$Distance > 0),]$bias) # p-value = 0.01584
ks.test(TA.GC.SNS.stochastic.assymetry[which(TA.GC.SNS.stochastic.assymetry$Distance < 0),]$bias, TA.GC.SNS.stochastic.assymetry[which(TA.GC.SNS.stochastic.assymetry$Distance > 0),]$bias) # p-value = 0.6077
#
ks.test(CG.AT.IniSeq.low.bias$strand.bias, CG.AT.IniSeq.medium.bias$strand.bias) # p-value = 0.0006709
ks.test(CG.AT.IniSeq.high.bias$strand.bias, CG.AT.IniSeq.medium.bias$strand.bias) # p-value = 0.003068
#
ks.test(CG.GC.IniSeq.low.bias$strand.bias, CG.GC.IniSeq.medium.bias$strand.bias) # p-value = 0.0001932
ks.test(CG.GC.IniSeq.high.bias$strand.bias, CG.GC.IniSeq.medium.bias$strand.bias) # p-value = 0.01638
#
ks.test(CG.TA.IniSeq.low.bias$strand.bias, CG.TA.IniSeq.medium.bias$strand.bias) # p-value = 0.3927
ks.test(CG.TA.IniSeq.high.bias$strand.bias, CG.TA.IniSeq.medium.bias$strand.bias) # p-value = 0.06809
#
ks.test(TA.AT.IniSeq.low.bias$strand.bias, TA.AT.IniSeq.medium.bias$strand.bias) #  p-value = 0.1777
ks.test(TA.AT.IniSeq.high.bias$strand.bias, TA.AT.IniSeq.medium.bias$strand.bias) # p-value = 0.27
#
ks.test(TA.CG.IniSeq.low.bias$strand.bias, TA.CG.IniSeq.medium.bias$strand.bias) #  p-value < 2.2e-16
ks.test(TA.CG.IniSeq.high.bias$strand.bias, TA.CG.IniSeq.medium.bias$strand.bias) # p-value = 2.531e-14
#
ks.test(TA.GC.IniSeq.low.bias$strand.bias, TA.GC.IniSeq.medium.bias$strand.bias) #  p-value = 1.223e-05
ks.test(TA.GC.IniSeq.high.bias$strand.bias, TA.GC.IniSeq.medium.bias$strand.bias) # p-value = 0.3275
