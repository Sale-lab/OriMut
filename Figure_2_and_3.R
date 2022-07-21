
# P. Murat, MRC Laboratory of Molecular Biology, July 2022
# Code for generating Figures 2 and 3 with associated supporting Figures
# in Murat et al. DNA replication initiation shapes the mutational landscape and expression of the human genome

library(BSgenome.Hsapiens.UCSC.hg38)
library(stringr)
library(MutationalPatterns)
library(NMF)
library(ggpubr)
library(reshape2)
library(lsa)
library(RColorBrewer)
library(drc)
library(tidyr)

setwd("/Volumes/pmurat/OriMut")

############################################################################
############################################################################
###               De novo mutational signature extraction                ###
############################################################################
############################################################################

# We show hereafter how we extract mutational signatures from SNP count matrices
# The same procedure was used for DNMs

############################################################################
# Prepare vcf files reporting population SNP around origins
# SNP at +- 10kb of isolated IniSeq origins
# bin by distance from the origins

ref_genome <- BSgenome.Hsapiens.UCSC.hg38

population.SNP.IniSeq.closest <- read.csv("./Figure_1/population.SNP.IniSeq.closest.10kb.bed", sep="\t")
colnames(population.SNP.IniSeq.closest) <- c("chr", "pos", "REF", "ALT", "SAO", "CAF", "ori.pos", "EFF", "dist")

pop.SNP.IniSeq.bin <- population.SNP.IniSeq.closest %>% filter(nchar(ALT) == 1) %>% mutate(dist = (as.numeric(cut(dist, breaks = seq(-10000, 10000, 100)))-100)*100)

# Load base composition information

IniSeq.inter.summary <- read.csv("./Figure_1/IniSeq.inter.10kb.base.composition.csv")

# Prepare vcf files corresponding to each dist bin from the center of orgins

for (i in unique(pop.SNP.IniSeq.bin$dist)) {
  vcf <- pop.SNP.IniSeq.bin %>% filter(dist == i) %>% mutate(chr = str_replace(chr, 'chr', '')) %>%
    mutate(ID = paste("rs", pos, sep = ""), QUAL = ".", FILTER = ".", INFO = ".", FORMAT = ".") %>%
    dplyr::select(chr, pos, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT)
  colnames(vcf) <- c('#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT')
  vcffile <- paste0("./Figure_2/VCF/Population_SNP_IniSeq/Bin_",i,".bed")
  write.table(vcf, vcffile, sep="\t", col.names = T, row.names = F, quote = F)
}

# manually remove bed called NA

# Add header to all bed files

# for file in *.bed
# do
# awk 'BEGIN{print"##fileformat=VCFv4.1"}1' $file > $file.vcf
# done

# Import and prepare data from multiple VCF files

vcfFiles.pop.SNP.IniSeq.bin <- grep("*.vcf$", dir("./Figure_2/VCF/Population_SNP_IniSeq/"), value = TRUE)
vcf_files <- paste("./Figure_2/VCF/Population_SNP_IniSeq/", vcfFiles.pop.SNP.IniSeq.bin, sep = "")
sample_names <- vcfFiles.pop.SNP.IniSeq.bin %>% str_replace(".bed.vcf", "") %>% str_replace("Bin_", "")

SNP_grl <- read_vcfs_as_granges(vcf_files, sample_names, ref_genome)
# Save grange object
saveRDS(SNP_grl, file = "./Figure_2/IniSeq_SNP_grange.rds")
SNP_grl <- readRDS("./Figure_2/IniSeq_SNP_grange.rds")

############################################################################
# Compute 96 trinucleotide mutation count matrix

mut.mat <- mut_matrix(vcf_list = SNP_grl, ref_genome = ref_genome)
# Save
saveRDS(mut.mat, file = "./Figure_2/IniSeq_SNP_mutation_matrix.rds")
mut.mat <- readRDS("./Figure_2/IniSeq_SNP_mutation_matrix.rds")

# add a small psuedocount to the mutation count matrix
mut.mat.2 <- mut.mat + 0.0001

############################################################################
# Extract mutational signatures

# Use nmf function to determine the number of signature to extract

estimate <- nmf(mut.mat.2, rank = 2:8, method = "brunet", nrun = 200, seed = 123456, .opt = "v-p")

# The cophenetic correlation coefficient starts decreasing for a rank of 3
# We then extract 3 signatures

# Extract signatures

nmf_res <- extract_signatures(mut.mat.2, rank = 3, nrun = 1000, single_core = TRUE)
colnames(nmf_res$signatures) <- c("Signature C", "Signature B", "Signature A")
rownames(nmf_res$contribution) <- c("Signature C", "Signature B", "Signature A")
# Save results
saveRDS(nmf_res, file = "./Figure_2/Signatures_IniSeq_SNP.rds")
nmf_res <- readRDS("./Figure_2/Signatures_IniSeq_SNP.rds")

pdf("./Rplot/Fig_2C.pdf", width=8, height=5, useDingbats=FALSE)
plot_96_profile(nmf_res$signatures, condensed = TRUE, ymax = 0.1)
dev.off()

# Compute signature exposures

exp.ori <- as.data.frame(nmf_res$contribution)
exp.ori <- as.data.frame(t(exp.ori)) %>% add_rownames(var = "dist")
exp.ori$dist <- as.integer(as.character(exp.ori$dist))

Sig.A <- exp.ori %>% mutate(fraction = `Signature A`) %>% mutate(signature = "SBS A")  %>% dplyr::select(dist, fraction, signature) 
Sig.B <- exp.ori %>% mutate(fraction = `Signature B`) %>% mutate(signature = "SBS B")  %>% dplyr::select(dist, fraction, signature) 
Sig.C <- exp.ori %>% mutate(fraction = `Signature C`) %>% mutate(signature = "SBS C")  %>% dplyr::select(dist, fraction, signature) 
dat.ori <- rbind(Sig.A, Sig.B, Sig.C)

exposure <- ggplot(dat.ori, aes(x = dist, y = fraction)) +
  geom_col(aes(fill = signature), width = 100, alpha = 0.8) + xlim(-5000,5000) + ylim(0,0.525) +
  scale_fill_manual(values = c("#D4D2D2", "#2EBAED", "#DE1C14")) + ggtitle("Population SNP / IniSeq origins") + ylab("Relative contribution") + xlab("Distance from origins (bp)") +
  theme_classic() + theme(aspect.ratio=1)

pdf("./Rplot/Fig_2D.pdf", width=8, height=5, useDingbats=FALSE)
exposure
dev.off()

############################################################################
# Assess exposure of SBS A, B and C to constitutive origins when binned by efficiency

quantile(population.SNP.IniSeq.closest$EFF, c(0.333, 0.666))
#    33.3%    66.6% 
# low    medium  high
#   0.8728     0.9297 

population.SNP.IniSeq.low <- population.SNP.IniSeq.closest %>% filter(EFF < 0.8728)  %>% filter(nchar(ALT) == 1) %>% mutate(dist = (as.numeric(cut(dist, breaks = seq(-10000, 10000, 100)))-100)*100)
population.SNP.IniSeq.medium <- population.SNP.IniSeq.closest %>% filter(EFF >= 0.8728 & EFF < 0.9297) %>% filter(nchar(ALT) == 1) %>% mutate(dist = (as.numeric(cut(dist, breaks = seq(-10000, 10000, 100)))-100)*100)
population.SNP.IniSeq.high <- population.SNP.IniSeq.closest %>% filter(EFF >= 0.9297) %>% filter(nchar(ALT) == 1) %>% mutate(dist = (as.numeric(cut(dist, breaks = seq(-10000, 10000, 100)))-100)*100)

# Prepare vcf files

for (i in unique(population.SNP.IniSeq.low$dist)) {
  vcf <- population.SNP.IniSeq.low %>% filter(dist == i) %>% mutate(chr = str_replace(chr, 'chr', '')) %>%
    mutate(ID = paste("rs", pos, sep = ""), QUAL = ".", FILTER = ".", INFO = ".", FORMAT = ".") %>%
    dplyr::select(chr, pos, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT)
  colnames(vcf) <- c('#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT')
  vcffile <- paste0("./Figure_2/VCF/Population_SNP_IniSeq_low/Bin_",i,".bed")
  write.table(vcf, vcffile, sep="\t", col.names = T, row.names = F, quote = F)
}

for (i in unique(population.SNP.IniSeq.medium$dist)) {
  vcf <- population.SNP.IniSeq.medium %>% filter(dist == i) %>% mutate(chr = str_replace(chr, 'chr', '')) %>%
    mutate(ID = paste("rs", pos, sep = ""), QUAL = ".", FILTER = ".", INFO = ".", FORMAT = ".") %>%
    dplyr::select(chr, pos, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT)
  colnames(vcf) <- c('#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT')
  vcffile <- paste0("./Figure_2/VCF/Population_SNP_IniSeq_medium/Bin_",i,".bed")
  write.table(vcf, vcffile, sep="\t", col.names = T, row.names = F, quote = F)
}

for (i in unique(population.SNP.IniSeq.high$dist)) {
  vcf <- population.SNP.IniSeq.high %>% filter(dist == i) %>% mutate(chr = str_replace(chr, 'chr', '')) %>%
    mutate(ID = paste("rs", pos, sep = ""), QUAL = ".", FILTER = ".", INFO = ".", FORMAT = ".") %>%
    dplyr::select(chr, pos, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT)
  colnames(vcf) <- c('#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT')
  vcffile <- paste0("./Figure_2/VCF/Population_SNP_IniSeq_high/Bin_",i,".bed")
  write.table(vcf, vcffile, sep="\t", col.names = T, row.names = F, quote = F)
}

# Remove NA file and add header to all bed files

for file in *.bed
do
awk 'BEGIN{print"##fileformat=VCFv4.1"}1' $file > $file.vcf
done

# Import and prepare data from multiple VCF files

vcfFiles.pop.SNP.IniSeq.low <- grep("*.vcf$", dir("./Figure_2/VCF/Population_SNP_IniSeq_low/"), value = TRUE)
vcf_files <- paste("./Figure_2/VCF/Population_SNP_IniSeq_low/", vcfFiles.pop.SNP.IniSeq.low, sep = "")
sample_names <- vcfFiles.pop.SNP.IniSeq.low %>% str_replace(".bed.vcf", "") %>% str_replace("Bin_", "")
SNP_grl_IniSeq_low <- read_vcfs_as_granges(vcf_files, sample_names, ref_genome)
# Save grange object
saveRDS(SNP_grl_IniSeq_low, file = "./Figure_2/IniSeq_low_SNP_grange.rds")

vcfFiles.pop.SNP.IniSeq.medium <- grep("*.vcf$", dir("./Figure_2/VCF/Population_SNP_IniSeq_medium/"), value = TRUE)
vcf_files <- paste("./Figure_2/VCF/Population_SNP_IniSeq_medium/", vcfFiles.pop.SNP.IniSeq.medium, sep = "")
sample_names <- vcfFiles.pop.SNP.IniSeq.medium %>% str_replace(".bed.vcf", "") %>% str_replace("Bin_", "")
SNP_grl_IniSeq_medium <- read_vcfs_as_granges(vcf_files, sample_names, ref_genome)
# Save grange object
saveRDS(SNP_grl_IniSeq_medium, file = "./Figure_2/IniSeq_medium_SNP_grange.rds")

vcfFiles.pop.SNP.IniSeq.high <- grep("*.vcf$", dir("./Figure_2/VCF/Population_SNP_IniSeq_high/"), value = TRUE)
vcf_files <- paste("./Figure_2/VCF/Population_SNP_IniSeq_high/", vcfFiles.pop.SNP.IniSeq.high, sep = "")
sample_names <- vcfFiles.pop.SNP.IniSeq.high %>% str_replace(".bed.vcf", "") %>% str_replace("Bin_", "")
SNP_grl_IniSeq_high <- read_vcfs_as_granges(vcf_files, sample_names, ref_genome)
# Save grange object
saveRDS(SNP_grl_IniSeq_high, file = "./Figure_2/IniSeq_high_SNP_grange.rds")

# Open grange objects
SNP_grl_IniSeq_low <- readRDS("./Figure_2/IniSeq_low_SNP_grange.rds")
SNP_grl_IniSeq_medium <- readRDS("./Figure_2/IniSeq_medium_SNP_grange.rds")
SNP_grl_IniSeq_high <- readRDS("./Figure_2/IniSeq_high_SNP_grange.rds")

# Exposure to IniSeq origins of low efficiency

IniSeq.exp.signatures <- nmf_res$signatures

mut.mat.Iniseq.low <- mut_matrix(vcf_list = SNP_grl_IniSeq_low, ref_genome = ref_genome)
fit.sig.Iniseq.low <- fit_to_signatures(mut.mat.Iniseq.low, IniSeq.exp.signatures)
#
Iniseq.low.exposure <- as.data.frame(fit.sig.Iniseq.low$contribution)
Iniseq.low.exposure <- as.data.frame(t(Iniseq.low.exposure)) %>% add_rownames(var = "dist")
Iniseq.low.exposure$dist <- as.integer(as.character(Iniseq.low.exposure$dist))
#
Sig.A.low <- Iniseq.low.exposure %>% mutate(fraction = `Signature A`) %>% mutate(signature = "Signature A")  %>% dplyr::select(dist, fraction, signature) 
Sig.B.low <- Iniseq.low.exposure %>% mutate(fraction = `Signature B`) %>% mutate(signature = "Signature B")  %>% dplyr::select(dist, fraction, signature) 
Sig.C.low <- Iniseq.low.exposure %>% mutate(fraction = `Signature C`) %>% mutate(signature = "Signature C")  %>% dplyr::select(dist, fraction, signature) 
dat.IniSeq.low <- rbind(Sig.A.low, Sig.B.low, Sig.C.low)
#
a <- ggplot(dat.IniSeq.low, aes(x = dist, y = fraction)) +
  geom_col(aes(fill = signature), width = 100, alpha = 0.8) + ylim(0, 0.20) + xlim(-5000,5000) +
  scale_fill_manual(values = c("#D4D2D2", "#2EBAED", "#DE1C14")) + ggtitle("Low") +
  ylab("Relative contribution") + xlab("Distance from origins (bp)") +
  theme_classic() + theme(aspect.ratio=1)

# Exposure to IniSeq origins of medium efficiency

mut.mat.Iniseq.medium <- mut_matrix(vcf_list = SNP_grl_IniSeq_medium, ref_genome = ref_genome)
fit.sig.Iniseq.medium <- fit_to_signatures(mut.mat.Iniseq.medium, IniSeq.exp.signatures)
#
Iniseq.medium.exposure <- as.data.frame(fit.sig.Iniseq.medium$contribution)
Iniseq.medium.exposure <- as.data.frame(t(Iniseq.medium.exposure)) %>% add_rownames(var = "dist")
Iniseq.medium.exposure$dist <- as.integer(as.character(Iniseq.medium.exposure$dist))
#
Sig.A.medium <- Iniseq.medium.exposure %>% mutate(fraction = `Signature A`) %>% mutate(signature = "Signature A")  %>% dplyr::select(dist, fraction, signature) 
Sig.B.medium <- Iniseq.medium.exposure %>% mutate(fraction = `Signature B`) %>% mutate(signature = "Signature B")  %>% dplyr::select(dist, fraction, signature) 
Sig.C.medium <- Iniseq.medium.exposure %>% mutate(fraction = `Signature C`) %>% mutate(signature = "Signature C")  %>% dplyr::select(dist, fraction, signature) 
dat.IniSeq.medium <- rbind(Sig.A.medium, Sig.B.medium, Sig.C.medium)
#
b <- ggplot(dat.IniSeq.medium, aes(x = dist, y = fraction)) +
  geom_col(aes(fill = signature), width = 100, alpha = 0.8) + ylim(0, 0.20) + xlim(-5000,5000) +
  scale_fill_manual(values = c("#D4D2D2", "#2EBAED", "#DE1C14")) + ggtitle("Medium") +
  ylab("Relative contribution") + xlab("Distance from origins (bp)") +
  theme_classic() + theme(aspect.ratio=1)

# Exposure to IniSeq origins of high efficiency
mut.mat.Iniseq.high <- mut_matrix(vcf_list = SNP_grl_IniSeq_high, ref_genome = ref_genome)
fit.sig.Iniseq.high <- fit_to_signatures(mut.mat.Iniseq.high, IniSeq.exp.signatures)
#
Iniseq.high.exposure <- as.data.frame(fit.sig.Iniseq.high$contribution)
Iniseq.high.exposure <- as.data.frame(t(Iniseq.high.exposure)) %>% add_rownames(var = "dist")
Iniseq.high.exposure$dist <- as.integer(as.character(Iniseq.high.exposure$dist))
#
Sig.A.high <- Iniseq.high.exposure %>% mutate(fraction = `Signature A`) %>% mutate(signature = "Signature A")  %>% dplyr::select(dist, fraction, signature) 
Sig.B.high <- Iniseq.high.exposure %>% mutate(fraction = `Signature B`) %>% mutate(signature = "Signature B")  %>% dplyr::select(dist, fraction, signature) 
Sig.C.high <- Iniseq.high.exposure %>% mutate(fraction = `Signature C`) %>% mutate(signature = "Signature C")  %>% dplyr::select(dist, fraction, signature) 
dat.IniSeq.high <- rbind(Sig.A.high, Sig.B.high, Sig.C.high)
#
c <- ggplot(dat.IniSeq.high, aes(x = dist, y = fraction)) +
  geom_col(aes(fill = signature), width = 100, alpha = 0.8) + ylim(0, 0.20) + xlim(-5000,5000) +
  scale_fill_manual(values = c("#D4D2D2", "#2EBAED", "#DE1C14")) + ggtitle("High") +
  ylab("Relative contribution") + xlab("Distance from origins (bp)") +
  theme_classic() + theme(aspect.ratio=1)

# Save mutation matrixes

saveRDS(mut.mat.Iniseq.low, file = "./Figure_2/IniSeq_low_SNP_mutation_matrix.rds")
saveRDS(mut.mat.Iniseq.medium, file = "./Figure_2/IniSeq_medium_SNP_mutation_matrix.rds")
saveRDS(mut.mat.Iniseq.high, file = "./Figure_2/IniSeq_high_SNP_mutation_matrix.rds")
mut.mat.Iniseq.low <- readRDS("./Figure_2/IniSeq_low_SNP_mutation_matrix.rds")
mut.mat.Iniseq.medium <- readRDS("./Figure_2/IniSeq_medium_SNP_mutation_matrix.rds")
mut.mat.Iniseq.high <- readRDS("./Figure_2/IniSeq_high_SNP_mutation_matrix.rds")

# Prepare and save plots

pdf("./Rplot/Fig_S2A.pdf", width=15, height=6, useDingbats=FALSE)
ggarrange(a, b, c, ncol = 3, nrow = 1)
dev.off()

############################################################################
# Assess initiation efficiency dependency of identified mutational signatures

Sig.A.low <- Iniseq.low.exposure %>% melt(id.vars = "dist") %>% filter(dist <= -2000 | dist >= 2000) %>% filter(variable == "Signature A") %>% mutate(eff = "low")
Sig.A.medium <- Iniseq.medium.exposure %>% melt(id.vars = "dist") %>% filter(dist <= -2000 | dist >= 2000) %>% filter(variable == "Signature A") %>% mutate(eff = "medium")
Sig.A.high <- Iniseq.high.exposure %>% melt(id.vars = "dist") %>% filter(dist <= -2000 | dist >= 2000) %>% filter(variable == "Signature A") %>% mutate(eff = "high")
Sig.A.plot <- rbind(Sig.A.low, Sig.A.medium, Sig.A.high) %>% mutate(chi = round(value*10))
df.table <- table(Sig.A.plot$eff, Sig.A.plot$chi)
barplot(df.table)
chisq.test(df.table)$p.value # 0.3655874
Sig.A.eff <- Sig.A.plot %>% mutate(Efficiency = fct_relevel(eff, "low", "medium", "high")) %>% 
  ggplot(aes(x = Efficiency, y = value, colour = Efficiency)) +
  stat_summary(fun.min = function(z) { quantile(z,0.25) }, fun.max = function(z) { quantile(z,0.75) }, fun = median, size = 1.5) +
  scale_color_manual(values = c("#2EBAED", "#ADCC54", "#D6604D")) +
  theme_classic() + theme(aspect.ratio=2, legend.position = "none", axis.title=element_text(size=14), axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1), axis.text.y = element_text(size = 12), axis.title.x=element_blank(), axis.title.y=element_text(vjust=4)) +
  ggtitle("SBS A - P = 0.3655874") + ylab("Relative contribution") + ylim(0.025,0.06)

Sig.B.low <- Iniseq.low.exposure %>% melt(id.vars = "dist") %>% filter(dist >= -500 & dist <= 500) %>% filter(variable == "Signature B") %>% mutate(eff = "low")
Sig.B.medium <- Iniseq.medium.exposure %>% melt(id.vars = "dist") %>% filter(dist >= -500 & dist <= 500) %>% filter(variable == "Signature B") %>% mutate(eff = "medium")
Sig.B.high <- Iniseq.high.exposure %>% melt(id.vars = "dist") %>% filter(dist >= -500 & dist <= 500) %>% filter(variable == "Signature B") %>% mutate(eff = "high")
Sig.B.plot <- rbind(Sig.B.low, Sig.B.medium, Sig.B.high) %>% mutate(chi = round(value*10))
df.table <- table(Sig.B.plot$eff, Sig.B.plot$chi)
barplot(df.table)
chisq.test(df.table)$p.value # 0.02913794
Sig.B.eff <- Sig.B.plot %>% mutate(Efficiency = fct_relevel(eff, "low", "medium", "high")) %>% 
  ggplot(aes(x = Efficiency, y = value, colour = Efficiency)) +
  stat_summary(fun.min = function(z) { quantile(z,0.25) }, fun.max = function(z) { quantile(z,0.75) }, fun = median, size = 1.5) +
  scale_color_manual(values = c("#2EBAED", "#ADCC54", "#D6604D")) +
  theme_classic() + theme(aspect.ratio=2, legend.position = "none", axis.title=element_text(size=14), axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1), axis.text.y = element_text(size = 12), axis.title.x=element_blank(), axis.title.y=element_text(vjust=4)) +
  ggtitle("SBS B - P = 0.02913794") + ylab("Relative contribution") + ylim(0.125,0.19)

Sig.C.low <- Iniseq.low.exposure %>% melt(id.vars = "dist") %>% filter((dist >= -2000 & dist <= -500) | (dist >= 500 & dist <= 2000)) %>% filter(variable == "Signature C") %>% mutate(eff = "low")
Sig.C.medium <- Iniseq.medium.exposure %>% melt(id.vars = "dist") %>% filter((dist >= -2000 & dist <= -500) | (dist >= 500 & dist <= 2000)) %>% filter(variable == "Signature C") %>% mutate(eff = "medium")
Sig.C.high <- Iniseq.high.exposure %>% melt(id.vars = "dist") %>% filter((dist >= -2000 & dist <= -500) | (dist >= 500 & dist <= 2000)) %>% filter(variable == "Signature C") %>% mutate(eff = "high")
Sig.C.plot <- rbind(Sig.C.low, Sig.C.medium, Sig.C.high) %>% mutate(chi = round(value*10))
df.table <- table(Sig.C.plot$eff, Sig.C.plot$chi)
barplot(df.table)
chisq.test(df.table)$p.value # 7.588665e-11
Sig.C.eff <- Sig.C.plot %>% mutate(Efficiency = fct_relevel(eff, "low", "medium", "high")) %>% 
  ggplot(aes(x = Efficiency, y = value, colour = Efficiency)) +
  stat_summary(fun.min = function(z) { quantile(z,0.25) }, fun.max = function(z) { quantile(z,0.75) }, fun = median, size = 1.5) +
  scale_color_manual(values = c("#2EBAED", "#ADCC54", "#D6604D")) +
  theme_classic() + theme(aspect.ratio=2, legend.position = "none", axis.title=element_text(size=14), axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1), axis.text.y = element_text(size = 12), axis.title.x=element_blank(), axis.title.y=element_text(vjust=4)) +
  ggtitle("SBS C - P = 7.588665e-11") + ylab("Relative contribution") + ylim(0,0.1)

pdf("./Rplot/Figure_S2B.pdf", width=12, height=4, useDingbats=FALSE)
ggarrange(Sig.A.eff, Sig.B.eff, Sig.C.eff, ncol = 3, nrow = 1)
dev.off()

############################################################################
# Compute exposure of identified signatures to SNS origins

population.SNP.SNS.core.closest <- read.csv("./Figure_1/population.SNP.SNS.core.closest.10kb.bed", sep="\t")
colnames(population.SNP.SNS.core.closest) <- c("chr", "pos", "REF", "ALT", "SAO", "CAF", "ori.pos", "dist")

population.SNP.SNS.stochastic.closest <- read.csv("./Figure_1/population.SNP.SNS.stochastic.closest.10kb.bed", sep="\t")
colnames(population.SNP.SNS.stochastic.closest) <- c("chr", "pos", "REF", "ALT", "SAO", "CAF", "ori.pos", "dist")

pop.SNP.SNS.core.bin <- population.SNP.SNS.core.closest %>% filter(nchar(ALT) == 1) %>% mutate(dist = (as.numeric(cut(dist, breaks = seq(-10000, 10000, 100)))-100)*100)
pop.SNP.SNS.stochastic.bin <- population.SNP.SNS.stochastic.closest %>% filter(nchar(ALT) == 1) %>% mutate(dist = (as.numeric(cut(dist, breaks = seq(-10000, 10000, 100)))-100)*100)

# Prepare vcf files

for (i in unique(pop.SNP.SNS.core.bin$dist)) {
  vcf <- pop.SNP.SNS.core.bin %>% filter(dist == i) %>% mutate(chr = str_replace(chr, 'chr', '')) %>%
    mutate(ID = paste("rs", pos, sep = ""), QUAL = ".", FILTER = ".", INFO = ".", FORMAT = ".") %>%
    dplyr::select(chr, pos, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT)
  colnames(vcf) <- c('#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT')
  vcffile <- paste0("./Figure_2/VCF/Population_SNP_SNS_core/Bin_",i,".bed")
  write.table(vcf, vcffile, sep="\t", col.names = T, row.names = F, quote = F)
}

for (i in unique(pop.SNP.SNS.stochastic.bin$dist)) {
  vcf <- pop.SNP.SNS.stochastic.bin %>% filter(dist == i) %>% mutate(chr = str_replace(chr, 'chr', '')) %>%
    mutate(ID = paste("rs", pos, sep = ""), QUAL = ".", FILTER = ".", INFO = ".", FORMAT = ".") %>%
    dplyr::select(chr, pos, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT)
  colnames(vcf) <- c('#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT')
  vcffile <- paste0("./Figure_2/VCF/Population_SNP_SNS_stochastic/Bin_",i,".bed")
  write.table(vcf, vcffile, sep="\t", col.names = T, row.names = F, quote = F)
}

# Remove NA file and add header to all bed files

# for file in *.bed
# do
# awk 'BEGIN{print"##fileformat=VCFv4.1"}1' $file > $file.vcf
# done

# Prepare grange objects

vcfFiles.pop.SNP.SNS.core <- grep("*.vcf$", dir("./Figure_2/VCF/Population_SNP_SNS_core/"), value = TRUE)
vcf_files <- paste("./Figure_2/VCF/Population_SNP_SNS_core/", vcfFiles.pop.SNP.SNS.core, sep = "")
sample_names <- vcfFiles.pop.SNP.SNS.core %>% str_replace(".bed.vcf", "") %>% str_replace("Bin_", "")
SNP_grl_SNS_core <- read_vcfs_as_granges(vcf_files, sample_names, ref_genome)
# Save grange object
saveRDS(SNP_grl_SNS_core, file = "./Figure_2/SNS_core_SNP_grange.rds")
SNP_grl_SNS_core <- readRDS("./Figure_2/SNS_core_SNP_grange.rds")

vcfFiles.pop.SNP.SNS.stochastic <- grep("*.vcf$", dir("./Figure_2/VCF/Population_SNP_SNS_stochastic/"), value = TRUE)
vcf_files <- paste("./Figure_2/VCF/Population_SNP_SNS_stochastic/", vcfFiles.pop.SNP.SNS.stochastic, sep = "")
sample_names <- vcfFiles.pop.SNP.SNS.stochastic %>% str_replace(".bed.vcf", "") %>% str_replace("Bin_", "")
SNP_grl_SNS_stochastic <- read_vcfs_as_granges(vcf_files, sample_names, ref_genome)
# Save grange object
saveRDS(SNP_grl_SNS_stochastic, file = "./Figure_2/SNS_stochastic_SNP_grange.rds")
SNP_grl_SNS_stochastic <- readRDS("./Figure_2/SNS_stochastic_SNP_grange.rds")

# Compute exposure to IniSeq mutational signatures

# Exposure to SNS core origins
mut.mat.SNS.core <- mut_matrix(vcf_list = SNP_grl_SNS_core, ref_genome = ref_genome)
fit.sig.SNS.core <- fit_to_signatures(mut.mat.SNS.core, IniSeq.exp.signatures)
#
SNS.core.exposure <- as.data.frame(fit.sig.SNS.core$contribution)
SNS.core.exposure <- as.data.frame(t(SNS.core.exposure)) %>% add_rownames(var = "dist")
SNS.core.exposure$dist <- as.integer(as.character(SNS.core.exposure$dist))
#
Sig.A.SNS.core <- SNS.core.exposure %>% mutate(fraction = `Signature A`) %>% mutate(signature = "Signature A")  %>% dplyr::select(dist, fraction, signature) 
Sig.B.SNS.core <- SNS.core.exposure %>% mutate(fraction = `Signature B`) %>% mutate(signature = "Signature B")  %>% dplyr::select(dist, fraction, signature) 
Sig.C.SNS.core <- SNS.core.exposure %>% mutate(fraction = `Signature C`) %>% mutate(signature = "Signature C")  %>% dplyr::select(dist, fraction, signature) 
dat.SNS.core <- rbind(Sig.A.SNS.core, Sig.B.SNS.core, Sig.C.SNS.core)
#
d <- ggplot(dat.SNS.core, aes(x = dist, y = fraction)) +
  geom_col(aes(fill = signature), width = 100, alpha = 0.8) +
  scale_fill_manual(values = c("#D4D2D2", "#2EBAED", "#DE1C14")) + ggtitle("Population SNP / SNS core") +
  ylab("Relative contribution") + xlab("Distance form origins (bp)") + xlim(-5000,5000) + ylim(0,0.525) +
  theme_classic() + theme(aspect.ratio=1)

# Exposure to SNS stochastic origins
mut.mat.SNS.stochastic <- mut_matrix(vcf_list = SNP_grl_SNS_stochastic, ref_genome = ref_genome)
fit.sig.SNS.stochastic <- fit_to_signatures(mut.mat.SNS.stochastic, IniSeq.exp.signatures)
#
SNS.stochastic.exposure <- as.data.frame(fit.sig.SNS.stochastic$contribution)
SNS.stochastic.exposure <- as.data.frame(t(SNS.stochastic.exposure)) %>% add_rownames(var = "dist")
SNS.stochastic.exposure$dist <- as.integer(as.character(SNS.stochastic.exposure$dist))
#
Sig.A.SNS.stochastic <- SNS.stochastic.exposure %>% mutate(fraction = `Signature A`) %>% mutate(signature = "Signature A")  %>% dplyr::select(dist, fraction, signature) 
Sig.B.SNS.stochastic <- SNS.stochastic.exposure %>% mutate(fraction = `Signature B`) %>% mutate(signature = "Signature B")  %>% dplyr::select(dist, fraction, signature) 
Sig.C.SNS.stochastic <- SNS.stochastic.exposure %>% mutate(fraction = `Signature C`) %>% mutate(signature = "Signature C")  %>% dplyr::select(dist, fraction, signature) 
dat.SNS.stochastic <- rbind(Sig.A.SNS.stochastic, Sig.B.SNS.stochastic, Sig.C.SNS.stochastic)
#
e <- ggplot(dat.SNS.stochastic, aes(x = dist, y = fraction)) +
  geom_col(aes(fill = signature), width = 100, alpha = 0.8) +
  scale_fill_manual(values = c("#D4D2D2", "#2EBAED", "#DE1C14")) + ggtitle("Population SNP / SNS stochastic") +
  ylab("Relative contribution") + xlab("Distance form origins (bp)") + xlim(-5000,5000) + ylim(0,0.525) +
  theme_classic() + theme(aspect.ratio=1)

# Save mutational matrices
saveRDS(mut.mat.SNS.core, file = "./Figure_2/SNS_core_SNP_mutation_matrix.rds")
saveRDS(mut.mat.SNS.stochastic, file = "./Figure_2/SNS_stochastic_SNP_mutation_matrix.rds")
mut.mat.SNS.core <- readRDS("./Figure_2/SNS_core_SNP_mutation_matrix.rds")
mut.mat.SNS.stochastic <- readRDS("./Figure_2/SNS_stochastic_SNP_mutation_matrix.rds")

pdf("./Rplot/Figure_S2C.pdf", width=15, height=4, useDingbats=FALSE)
ggarrange(exposure, d, e, ncol = 3, nrow = 1)
dev.off()

############################################################################
# De novo extraction of indel signatures

# Prepare vcf files reporting population INDELs around origins
# INDELs at +- 10kb of isolated IniSeq origins
# bin by distance from the origins

population.INDEL.IniSeq.closest <- read.csv("./Figure_1/population.INDEL.IniSeq.closest.10kb.bed", sep="\t", header = F)
colnames(population.INDEL.IniSeq.closest) <- c("chr", "pos", "REF", "ALT", "ori.pos", "EFF", "dist")

pop.INDEL.IniSeq.bin <- population.INDEL.IniSeq.closest %>% mutate(dist = (as.numeric(cut(dist, breaks = seq(-10000, 10000, 100)))-100)*100)

# Save vcf files corresponding to each dist bin from the center of orgins

for (i in unique(pop.INDEL.IniSeq.bin$dist)) {
  vcf <- pop.INDEL.IniSeq.bin %>% filter(dist == i) %>% mutate(chr = str_replace(chr, 'chr', '')) %>%
    mutate(ID = paste("rs", pos, sep = ""), QUAL = ".", FILTER = ".", INFO = ".", FORMAT = ".") %>%
    dplyr::select(chr, pos, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT)
  colnames(vcf) <- c('#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT')
  vcffile <- paste0("./Figure_2/VCF/Population_INDEL_IniSeq/Bin_",i,".bed")
  write.table(vcf, vcffile, sep="\t", col.names = T, row.names = F, quote = F)
}

# manually remove bed called NA

# Add header to all bed files

for file in *.bed
do
awk 'BEGIN{print"##fileformat=VCFv4.1"}1' $file > $file.vcf
done

# Import and prepare data from multiple VCF files

vcfFiles.pop.INDEL.IniSeq.bin <- grep("*.vcf$", dir("./Figure_2/VCF/Population_INDEL_IniSeq/"), value = TRUE)
vcf_files <- paste("./Figure_2/VCF/Population_INDEL_IniSeq/", vcfFiles.pop.INDEL.IniSeq.bin, sep = "")
sample_names <- vcfFiles.pop.INDEL.IniSeq.bin %>% str_replace(".bed.vcf", "") %>% str_replace("Bin_", "")
grl_pop_indel <- read_vcfs_as_granges(vcf_files, sample_names, ref_genome, type = "indel")

# Get INDEL context
indel_grl <- get_indel_context(grl_pop_indel, ref_genome)
# Save grange object
saveRDS(indel_grl, file = "./Figure_2/IniSeq_INDEL_grange.rds")
indel_grl <- readRDS("./Figure_2/IniSeq_INDEL_grange.rds")

# Count the number of INDEL per type
indel_counts <- count_indel_contexts(indel_grl)

# add a small pseudocount to the mutation count matrix

indel_counts_2 <- indel_counts + 0.0001

# Use nmf function to determine the number of signature to extract

estimate.INDEL <- nmf(indel_counts_2, rank = 2:8, method = "brunet", nrun = 200, seed = 123456, .opt = "v-p")
# suggests 3 individual INDEL signatures

# Extract signatures

nmf_res_INDEL <- extract_signatures(indel_counts_2, rank = 3, nrun = 1000, single_core = TRUE)
colnames(nmf_res_INDEL$signatures) <- c("Signature IDB", "Signature IDA", "Signature IDC")
rownames(nmf_res_INDEL$contribution) <- c("Signature IDB", "Signature IDA", "Signature IDC")
# Save results
saveRDS(nmf_res_INDEL, file = "./Figure_2/Signatures_IniSeq_INDEL.rds")
nmf_res_INDEL <- readRDS("./Figure_2/Signatures_IniSeq_INDEL.rds")

pdf("./Rplot/Figure_S2D.pdf", width=14, height=7, useDingbats=FALSE)
plot_indel_contexts(nmf_res_INDEL$signatures, condensed = TRUE)
dev.off()

# Compute signature exposures

exp.ori.INDEL <- as.data.frame(nmf_res_INDEL$contribution)
exp.ori.INDEL <- as.data.frame(t(exp.ori.INDEL)) %>% add_rownames(var = "dist")
exp.ori.INDEL$dist <- as.integer(as.character(exp.ori.INDEL$dist))

Sig.A <- exp.ori.INDEL %>% mutate(fraction = `Signature IDA`) %>% mutate(signature = "ID A")  %>% dplyr::select(dist, fraction, signature) 
Sig.B <- exp.ori.INDEL %>% mutate(fraction = `Signature IDB`) %>% mutate(signature = "ID B")  %>% dplyr::select(dist, fraction, signature) 
Sig.C <- exp.ori.INDEL %>% mutate(fraction = `Signature IDC`) %>% mutate(signature = "ID C")  %>% dplyr::select(dist, fraction, signature) 
dat.ori.INDEL <- rbind(Sig.A, Sig.B, Sig.C)

pdf("./Rplot/Fig_S2E.pdf", width=8, height=5, useDingbats=FALSE)
dat.ori.INDEL %>% mutate(signature = fct_relevel(signature, "Signature IDA", "Signature IDB", "Signature IDC")) %>% 
  ggplot(aes(x = dist, y = fraction)) + 
  geom_col(aes(fill = signature), width = 100, alpha = 0.8) + xlim(-5000,5000) +
  scale_fill_manual(values = c("#D4D2D2", "#2EBAED", "#DE1C14")) + ggtitle("Population INDEL / IniSeq origins") + ylab("Relative contribution") + xlab("Distance from origins (bp)") +
  theme_classic() + theme(aspect.ratio=1)
dev.off()

############################################################################
############################################################################
###        Contribution of known mutational signatures (COSMIC)          ###
###                      to mutagenesis at origins                       ###
############################################################################
############################################################################

############################################################################
# Find optimal contribution of COSMIC signatures to mutational landscape of IniSeq origins
# using strict refitting

COSMIC.signatures <- get_known_signatures()
COSMIC.fit.IniSeq <- fit_to_signatures_strict(mut.mat.2, COSMIC.signatures, max_delta = 0.004)

# Identify operating signatures

COSMIC.fit.IniSeq <- as.data.frame(COSMIC.fit.IniSeq$fit_res$contribution)
COSMIC.fit.IniSeq <- as.data.frame(t(COSMIC.fit.IniSeq)) %>% dplyr::select_if(colSums(.) != 0) %>% add_rownames(var = "dist")
COSMIC.fit.IniSeq$dist <- as.numeric(as.character(COSMIC.fit.IniSeq$dist))
COSMIC.exposure.IniSeq <- melt(COSMIC.fit.IniSeq, id.vars = "dist")

# Contributing signatures
# SBS1  SBS3  SBS5  SBS6  SBS7b SBS15 SBS19 SBS23 SBS24 SBS26 SBS30 SBS31 SBS39 SBS87

# Compute cosine similarity between extracted SBS A, B and C and contributing COSMIC signatures

Sig <- c("SBS1", "SBS3", "SBS5", "SBS6", "SBS7b", "SBS15", "SBS19", "SBS23", "SBS24", "SBS26", "SBS30", "SBS31", "SBS39", "SBS87")
COSMIC.cont <- as.data.frame(COSMIC.signatures) %>% dplyr::select(all_of(Sig))

Signature.A <- as.data.frame(nmf_res$signatures)$'Signature A'
Signature.B <- as.data.frame(nmf_res$signatures)$'Signature B'
Signature.C <- as.data.frame(nmf_res$signatures)$'Signature C'
# Normalise vectors
scale1 <- function(x) {x / sum(x)}
Signature.A <- scale1(Signature.A)
Signature.B <- scale1(Signature.B)
Signature.C <- scale1(Signature.C)

Sig.A.cos <- vector()
Sig.B.cos <- vector()
Sig.C.cos <- vector()
for (i in colnames(COSMIC.cont)) {
  cos <- COSMIC.cont %>% dplyr::select(all_of(i))
  p <- cosine(Signature.A, cos[,1])
  q <- cosine(Signature.B, cos[,1])
  r <- cosine(Signature.C, cos[,1])
  Sig.A.cos <- append(Sig.A.cos, p)
  Sig.B.cos <- append(Sig.B.cos, q)
  Sig.C.cos <- append(Sig.C.cos, r)
}
cos.mat <- cbind.data.frame(colnames(COSMIC.cont), Sig.C.cos, Sig.B.cos, Sig.A.cos)
colnames(cos.mat) <- c("COSMIC", "Signature C", "Signature B", "Signature A")
cos.mat.melt <- melt(cos.mat)

pdf("./Rplot/Fig_S3A.pdf", width=8, height=5, useDingbats=FALSE)
ggplot(data = cos.mat.melt, aes(x=fct_inorder(COSMIC), y=variable, fill=value)) + geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0.5, limit = c(0,1), space = "Lab", name="Cosine\nSimilarity") +
  theme_minimal() + theme(aspect.ratio=0.1, axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1), axis.text.y = element_text(size = 12), axis.title.x=element_blank(), axis.title.y=element_blank())
dev.off()

# Define set of COSMIC signatures according to their similarity to SBS A, B and C

cos.sim.match <- cbind.data.frame(COSMIC = cos.mat$COSMIC, match = colnames(cos.mat[,2:4])[apply(cos.mat[,2:4],1,which.max)])

setA <- (cos.sim.match %>% filter(match == "Signature A"))$COSMIC # Signatures SBS1, SBS6, SBS87 
setB <- (cos.sim.match %>% filter(match == "Signature B"))$COSMIC # Signatures SBS7b, SBS15, SBS19, SBS23, SBS24, SBS31
setC <- (cos.sim.match %>% filter(match == "Signature C"))$COSMIC # Signatures SBS3, SBS5, SBS26, SBS30, SBS39
# Two signatures are misassigned
# SBS6 shows similarity with Signature A and B but present a profile more representative of Signature C --> manually assign to Signature C
# SBS87 shows a higher similarity to Signature A than B but very close --> manually assign to Signature B 

# Assess exposure of signature set A
COSMIC.exposure.IniSeq.setA <- COSMIC.exposure.IniSeq %>% filter(variable %in% setA[!setA %in% c("SBS6", "SBS87")])

setA.plot <- COSMIC.exposure.IniSeq.setA %>% mutate(signature = fct_relevel(variable, "SBS1")) %>%
  ggplot(aes(x = dist, y = value)) +
  geom_col(aes(fill = signature), width = 100, alpha = 0.8) + ylim(0, 1500) +  xlim(-5000,5000) +
  scale_fill_manual(values = c("#D4D2D2")) +
  ggtitle("set A") + xlab("Distance from origins (bp)") + ylab("Absolute contribution\n(number of mutations") + labs(fill = "COSMIC\nsignatures") +
  theme_classic() + theme(aspect.ratio=1)

# Assess exposure of signature set B
COSMIC.exposure.IniSeq.setB <- COSMIC.exposure.IniSeq %>% filter(variable %in% c(setB, "SBS87"))

display.brewer.pal(n = 9, name = 'Blues')
brewer.pal(n = 9, name = "Blues")
setB.plot <- COSMIC.exposure.IniSeq.setB %>%
  ggplot(aes(x = dist, y = value)) + ylim(0, 8500) + xlim(-5000,5000) +
  geom_col(aes(fill = variable), width = 100, alpha = 0.8) +
  scale_fill_manual(values = c("#000000", "#08306B","#08519C", "#2171B5", "#4292C6", "#6BAED6", "#9ECAE1")) +
  ggtitle("set B") + xlab("Distance from origins (bp)") + ylab("Absolute contribution\n(number of mutations") + labs(fill = "COSMIC\nsignatures") +
  theme_classic() + theme(aspect.ratio=1)

# Assess exposure of signature set C
COSMIC.exposure.IniSeq.setC <- COSMIC.exposure.IniSeq %>% filter(variable %in% c(setC, "SBS6"))# %>% filter(value >= 200)

display.brewer.pal(n = 9, name = 'YlOrRd')
brewer.pal(n = 9, name = "YlOrRd")
setC.plot <- COSMIC.exposure.IniSeq.setC %>%
  ggplot(aes(x = dist, y = value)) +
  geom_col(aes(fill = variable), width = 100, alpha = 0.8) + ylim(0, 8500) + xlim(-5000,5000) +
  scale_fill_manual(values = c("#999999", "#D4D2D2", "#800026", "#FC4E2A", "#FD8D3C", "#FED976")) +
  ggtitle("set C") + xlab("Distance from origins (bp)") + ylab("Absolute contribution\n(number of mutations") + labs(fill = "COSMIC\nsignatures") +
  theme_classic() + theme(aspect.ratio=1)

pdf("./Rplot/Fig_2E_2F_2G.pdf", width=15, height=6, useDingbats=FALSE)
ggarrange(setA.plot, setB.plot, setC.plot, ncol = 3, nrow = 1)
dev.off()

############################################################################
# Contribution of COSMIC signatures to IniSeq origins binned by efficiency

# Fit COSMIC signatures to IniSeq origins binned by efficiency
COSMIC.fit.IniSeq.low <- fit_to_signatures(mut.mat.Iniseq.low, as.matrix(COSMIC.cont))
COSMIC.fit.IniSeq.medium <- fit_to_signatures(mut.mat.Iniseq.medium, as.matrix(COSMIC.cont))
COSMIC.fit.IniSeq.high <- fit_to_signatures(mut.mat.Iniseq.high, as.matrix(COSMIC.cont))

# Assess contribution of each set of COSMIC signatures to the mutational lanscape of IniSeq origins
# Low efficiency
COSMIC.exposure.IniSeq.low <- as.data.frame(COSMIC.fit.IniSeq.low$contribution)
COSMIC.exposure.IniSeq.low <- as.data.frame(t(COSMIC.exposure.IniSeq.low)) %>% dplyr::select_if(colSums(.) != 0) %>% add_rownames(var = "dist")
COSMIC.exposure.IniSeq.low$dist <- as.numeric(as.character(COSMIC.exposure.IniSeq.low$dist))
COSMIC.exposure.IniSeq.low <- melt(COSMIC.exposure.IniSeq.low, id.vars = "dist")
# Medium efficiency
COSMIC.exposure.IniSeq.medium <- as.data.frame(COSMIC.fit.IniSeq.medium$contribution)
COSMIC.exposure.IniSeq.medium <- as.data.frame(t(COSMIC.exposure.IniSeq.medium)) %>% dplyr::select_if(colSums(.) != 0) %>% add_rownames(var = "dist")
COSMIC.exposure.IniSeq.medium$dist <- as.numeric(as.character(COSMIC.exposure.IniSeq.medium$dist))
COSMIC.exposure.IniSeq.medium <- melt(COSMIC.exposure.IniSeq.medium, id.vars = "dist")
# High efficiency
COSMIC.exposure.IniSeq.high <- as.data.frame(COSMIC.fit.IniSeq.high$contribution)
COSMIC.exposure.IniSeq.high <- as.data.frame(t(COSMIC.exposure.IniSeq.high)) %>% dplyr::select_if(colSums(.) != 0) %>% add_rownames(var = "dist")
COSMIC.exposure.IniSeq.high$dist <- as.numeric(as.character(COSMIC.exposure.IniSeq.high$dist))
COSMIC.exposure.IniSeq.high <- melt(COSMIC.exposure.IniSeq.high, id.vars = "dist")

# For COSMIC signatures Set A
COSMIC.exposure.IniSeq.low.setA <- COSMIC.exposure.IniSeq.low %>% filter(variable %in% c("SBS1"))
setA.low.plot <- COSMIC.exposure.IniSeq.low.setA %>%
  ggplot(aes(x = dist, y = value)) + ylim(0,600) +
  geom_col(aes(fill = variable), width = 100, alpha = 0.8) +
  scale_fill_manual(values = c("#D4D2D2", "#ADCC54", "#2EBAED")) + xlim(-2500,2500) +
  ggtitle("COSMIC signatures set A around low IniSeq origins") + xlab("Distance from origins (bp)") + ylab("Absolute contribution\n(number of mutations") + labs(fill = "COSMIC\nsignatures") +
  theme_classic() + theme(aspect.ratio=1.5)
COSMIC.exposure.IniSeq.medium.setA <- COSMIC.exposure.IniSeq.medium %>% filter(variable %in% c("SBS1"))
setA.medium.plot <- COSMIC.exposure.IniSeq.medium.setA %>%
  ggplot(aes(x = dist, y = value)) + ylim(0,600) +
  geom_col(aes(fill = variable), width = 100, alpha = 0.8) +
  scale_fill_manual(values = c("#D4D2D2", "#ADCC54", "#2EBAED")) + xlim(-2500,2500) +
  ggtitle("COSMIC signatures set A around medium IniSeq origins") + xlab("Distance from origins (bp)") + ylab("Absolute contribution\n(number of mutations") + labs(fill = "COSMIC\nsignatures") +
  theme_classic() + theme(aspect.ratio=1.5)
COSMIC.exposure.IniSeq.high.setA <- COSMIC.exposure.IniSeq.high %>% filter(variable %in% c("SBS1"))
setA.high.plot <- COSMIC.exposure.IniSeq.high.setA %>%
  ggplot(aes(x = dist, y = value)) + ylim(0,600) +
  geom_col(aes(fill = variable), width = 100, alpha = 0.8) +
  scale_fill_manual(values = c("#D4D2D2", "#ADCC54", "#2EBAED")) + xlim(-2500,2500) +
  ggtitle("COSMIC signatures set A around high IniSeq origins") + xlab("Distance from origins (bp)") + ylab("Absolute contribution\n(number of mutations") + labs(fill = "COSMIC\nsignatures") +
  theme_classic() + theme(aspect.ratio=1.5)

# For COSMIC signatures Set B
COSMIC.exposure.IniSeq.low.setB <- COSMIC.exposure.IniSeq.low %>% filter(variable %in% c(setB, "SBS87"))
setB.low.plot <- COSMIC.exposure.IniSeq.low.setB %>%
  ggplot(aes(x = dist, y = value)) + ylim(0,3000) +
  geom_col(aes(fill = variable), width = 100, alpha = 0.8) + xlim(-2500,2500) +
  scale_fill_manual(values = c("#000000", "#08306B","#08519C", "#2171B5", "#4292C6", "#6BAED6", "#9ECAE1")) +
  ggtitle("COSMIC signatures set B around low IniSeq origins") + xlab("Distance from origins (bp)") + ylab("Absolute contribution\n(number of mutations") + labs(fill = "COSMIC\nsignatures") +
  theme_classic() + theme(aspect.ratio=1.5)
COSMIC.exposure.IniSeq.medium.setB <- COSMIC.exposure.IniSeq.medium %>% filter(variable %in% c(setB, "SBS87"))
setB.medium.plot <- COSMIC.exposure.IniSeq.medium.setB %>%
  ggplot(aes(x = dist, y = value)) + ylim(0,3000) +
  geom_col(aes(fill = variable), width = 100, alpha = 0.8) + xlim(-2500,2500) +
  scale_fill_manual(values = c("#000000", "#08306B","#08519C", "#2171B5", "#4292C6", "#6BAED6", "#9ECAE1")) +
  ggtitle("COSMIC signatures set B around medium IniSeq origins") + xlab("Distance from origins (bp)") + ylab("Absolute contribution\n(number of mutations") + labs(fill = "COSMIC\nsignatures") +
  theme_classic() + theme(aspect.ratio=1.5)
COSMIC.exposure.IniSeq.high.setB <- COSMIC.exposure.IniSeq.high %>% filter(variable %in% c(setB, "SBS87"))
setB.high.plot <- COSMIC.exposure.IniSeq.high.setB %>%
  ggplot(aes(x = dist, y = value)) + ylim(0,3000) +
  geom_col(aes(fill = variable), width = 100, alpha = 0.8) + xlim(-2500,2500) +
  scale_fill_manual(values = c("#000000", "#08306B","#08519C", "#2171B5", "#4292C6", "#6BAED6", "#9ECAE1")) +
  ggtitle("COSMIC signatures set B around high IniSeq origins") + xlab("Distance from origins (bp)") + ylab("Absolute contribution\n(number of mutations") + labs(fill = "COSMIC\nsignatures") +
  theme_classic() + theme(aspect.ratio=1.5)

# For COSMIC signatures Set C
COSMIC.exposure.IniSeq.low.setC <- COSMIC.exposure.IniSeq.low %>% filter(variable %in% c(setC, "SBS6")) 
setC.low.plot <- COSMIC.exposure.IniSeq.low.setC %>%
  ggplot(aes(x = dist, y = value)) + ylim(0,3000) + xlim(-2500,2500) +
  geom_col(aes(fill = variable), width = 100, alpha = 0.8) + 
  scale_fill_manual(values = c("#999999", "#D4D2D2", "#800026", "#FC4E2A", "#FD8D3C", "#FED976")) +
  ggtitle("COSMIC signatures set C around low IniSeq origins") + xlab("Distance from origins (bp)") + ylab("Absolute contribution\n(number of mutations") + labs(fill = "COSMIC\nsignatures") +
  theme_classic() + theme(aspect.ratio=1.5)
COSMIC.exposure.IniSeq.medium.setC <- COSMIC.exposure.IniSeq.medium %>% filter(variable %in% c(setC, "SBS6")) 
setC.medium.plot <- COSMIC.exposure.IniSeq.medium.setC %>%
  ggplot(aes(x = dist, y = value)) + ylim(0,3000) + xlim(-2500,2500) +
  geom_col(aes(fill = variable), width = 100, alpha = 0.8) +
  scale_fill_manual(values = c("#999999", "#D4D2D2", "#800026", "#FC4E2A", "#FD8D3C", "#FED976")) +
  ggtitle("COSMIC signatures set C around medium IniSeq origins") + xlab("Distance from origins (bp)") + ylab("Absolute contribution\n(number of mutations") + labs(fill = "COSMIC\nsignatures") +
  theme_classic() + theme(aspect.ratio=1.5)
COSMIC.exposure.IniSeq.high.setC <- COSMIC.exposure.IniSeq.high %>% filter(variable %in% c(setC, "SBS6")) 
setC.high.plot <- COSMIC.exposure.IniSeq.high.setC %>%
  ggplot(aes(x = dist, y = value)) + ylim(0,3000) + xlim(-2500,2500) +
  geom_col(aes(fill = variable), width = 100, alpha = 0.8) +
  scale_fill_manual(values = c("#999999", "#D4D2D2", "#800026", "#FC4E2A", "#FD8D3C", "#FED976")) +
  ggtitle("COSMIC signatures set C around high IniSeq origins") + xlab("Distance from origins (bp)") + ylab("Absolute contribution\n(number of mutations") + labs(fill = "COSMIC\nsignatures") +
  theme_classic() + theme(aspect.ratio=1.5)

# Save plot

pdf("./Rplot/Fig_S3B.pdf", width=12, height=12, useDingbats=FALSE)
ggarrange(setA.high.plot, setB.high.plot, setC.high.plot, 
          setA.medium.plot, setB.medium.plot, setC.medium.plot, 
          setA.low.plot, setB.low.plot, setC.low.plot, 
          ncol = 3, nrow = 3)
dev.off()

############################################################################
# Contribution of COSMIC signatures to SNS core and stochastic origins

COSMIC.fit.SNS.core <- fit_to_signatures(mut.mat.SNS.core, as.matrix(COSMIC.cont))
COSMIC.fit.SNS.stochastic <- fit_to_signatures(mut.mat.SNS.stochastic, as.matrix(COSMIC.cont))
# Assess contribution of each set of COSMIC signatures to the mutational landscape of SNS origins
# SNS core
COSMIC.exposure.SNS.core <- as.data.frame(COSMIC.fit.SNS.core$contribution)
COSMIC.exposure.SNS.core <- as.data.frame(t(COSMIC.exposure.SNS.core)) %>% dplyr::select_if(colSums(.) != 0) %>% add_rownames(var = "dist")
COSMIC.exposure.SNS.core$dist <- as.numeric(as.character(COSMIC.exposure.SNS.core$dist))
COSMIC.exposure.SNS.core <- melt(COSMIC.exposure.SNS.core, id.vars = "dist")
# SNS stochastic
COSMIC.exposure.SNS.stochastic <- as.data.frame(COSMIC.fit.SNS.stochastic$contribution)
COSMIC.exposure.SNS.stochastic <- as.data.frame(t(COSMIC.exposure.SNS.stochastic)) %>% dplyr::select_if(colSums(.) != 0) %>% add_rownames(var = "dist")
COSMIC.exposure.SNS.stochastic$dist <- as.numeric(as.character(COSMIC.exposure.SNS.stochastic$dist))
COSMIC.exposure.SNS.stochastic <- melt(COSMIC.exposure.SNS.stochastic, id.vars = "dist")

# For COSMIC signatures Set A
COSMIC.exposure.SNS.core.setA <- COSMIC.exposure.SNS.core %>% filter(variable %in% c("SBS1"))
setA.SNS.core.plot <- COSMIC.exposure.SNS.core.setA %>%
  ggplot(aes(x = dist, y = value)) + ylim(0,2000) + xlim(-2500,2500) +
  geom_col(aes(fill = variable), width = 100, alpha = 0.8) +
  scale_fill_manual(values = c("#D4D2D2", "#ADCC54", "#2EBAED")) +
  ggtitle("COSMIC signatures set A around SNS core origins") + xlab("Distance from origins (bp)") + ylab("Absolute contribution\n(number of mutations") + labs(fill = "COSMIC\nsignatures") +
  theme_classic() + theme(aspect.ratio=1.5)
COSMIC.exposure.SNS.stochastic.setA <- COSMIC.exposure.SNS.stochastic %>% filter(variable %in% c("SBS1"))
setA.SNS.stochastic.plot <- COSMIC.exposure.SNS.stochastic.setA %>%
  ggplot(aes(x = dist, y = value)) + ylim(0,2000) + xlim(-2500,2500) +
  geom_col(aes(fill = variable), width = 100, alpha = 0.8) +
  scale_fill_manual(values = c("#D4D2D2", "#ADCC54", "#2EBAED")) +
  ggtitle("COSMIC signatures set A around SNS stochastic origins") + xlab("Distance from origins (bp)") + ylab("Absolute contribution\n(number of mutations") + labs(fill = "COSMIC\nsignatures") +
  theme_classic() + theme(aspect.ratio=1.5)

# For COSMIC signatures Set B
COSMIC.exposure.SNS.core.setB <- COSMIC.exposure.SNS.core %>% filter(variable %in% c(setB, "SBS87"))
setB.SNS.core.plot <- COSMIC.exposure.SNS.core.setB %>%
  ggplot(aes(x = dist, y = value)) + ylim(0, 8000) + xlim(-2500,2500) +
  geom_col(aes(fill = variable), width = 100, alpha = 0.8) +
  scale_fill_manual(values = c("#000000","#08519C", "#2171B5", "#4292C6", "#6BAED6", "#9ECAE1")) +
  ggtitle("COSMIC signatures set B around SNS core origins") + xlab("Distance from origins (bp)") + ylab("Absolute contribution\n(number of mutations") + labs(fill = "COSMIC\nsignatures") +
  theme_classic() + theme(aspect.ratio=1.5)
COSMIC.exposure.SNS.stochastic.setB <- COSMIC.exposure.SNS.stochastic %>% filter(variable %in% c(setB, "SBS87"))
setB.SNS.stochastic.plot <- COSMIC.exposure.SNS.stochastic.setB %>%
  ggplot(aes(x = dist, y = value)) + ylim(0, 8000) + xlim(-2500,2500) +
  geom_col(aes(fill = variable), width = 100, alpha = 0.8) +
  scale_fill_manual(values = c("#000000","#08519C", "#6BAED6", "#9ECAE1")) +
  ggtitle("COSMIC signatures set B around SNS stochastic origins") + xlab("Distance from origins (bp)") + ylab("Absolute contribution\n(number of mutations") + labs(fill = "COSMIC\nsignatures") +
  theme_classic() + theme(aspect.ratio=1.5)

# For COSMIC signatures Set C
COSMIC.exposure.SNS.core.setC <- COSMIC.exposure.SNS.core %>% filter(variable %in% c(setC, "SBS6")) 
setC.SNS.core.plot <- COSMIC.exposure.SNS.core.setC %>%
  ggplot(aes(x = dist, y = value)) + ylim(0,11000) + xlim(-2500,2500) +
  geom_col(aes(fill = variable), width = 100, alpha = 0.8) +
  scale_fill_manual(values = c("#999999", "#D4D2D2", "#800026", "#FC4E2A", "#FD8D3C", "#FED976", "#FFEDA0")) +
  ggtitle("COSMIC signatures set C around SNS core origins") + xlab("Distance from origins (bp)") + ylab("Absolute contribution\n(number of mutations") + labs(fill = "COSMIC\nsignatures") +
  theme_classic() + theme(aspect.ratio=1.5)
COSMIC.exposure.SNS.stochastic.setC <- COSMIC.exposure.SNS.stochastic %>% filter(variable %in% c(setC, "SBS6")) 
setC.SNS.stochastic.plot <- COSMIC.exposure.SNS.stochastic.setC %>%
  ggplot(aes(x = dist, y = value)) + ylim(0,11000) + xlim(-2500,2500) +
  geom_col(aes(fill = variable), width = 100, alpha = 0.8) +
  scale_fill_manual(values = c("#999999", "#D4D2D2", "#FD8D3C", "#FED976", "#FFEDA0")) +
  ggtitle("COSMIC signatures set C around SNS stochastic origins") + xlab("Distance from origins (bp)") + ylab("Absolute contribution\n(number of mutations") + labs(fill = "COSMIC\nsignatures") +
  theme_classic() + theme(aspect.ratio=1.5)

# Replot IniSeq origins with different scales for comparison

setA.plot <- COSMIC.exposure.IniSeq.setA %>% mutate(signature = fct_relevel(variable, "SBS1")) %>%
  ggplot(aes(x = dist, y = value)) +
  geom_col(aes(fill = signature), width = 100, alpha = 0.8) + ylim(0, 2000) + xlim(-2500,2500) +
  scale_fill_manual(values = c("#D4D2D2")) +
  ggtitle("set A") + xlab("Distance from origins (bp)") + ylab("Absolute contribution\n(number of mutations") + labs(fill = "COSMIC\nsignatures") +
  theme_classic() + theme(aspect.ratio=1.5)
setB.plot <- COSMIC.exposure.IniSeq.setB %>%
  ggplot(aes(x = dist, y = value)) + ylim(0, 8000) + xlim(-2500,2500) +
  geom_col(aes(fill = variable), width = 100, alpha = 0.8) +
  scale_fill_manual(values = c("#000000", "#08306B","#08519C", "#2171B5", "#4292C6", "#6BAED6", "#9ECAE1")) +
  ggtitle("set B") + xlab("Distance from origins (bp)") + ylab("Absolute contribution\n(number of mutations") + labs(fill = "COSMIC\nsignatures") +
  theme_classic() + theme(aspect.ratio=1.5)
setC.plot <- COSMIC.exposure.IniSeq.setC %>%
  ggplot(aes(x = dist, y = value)) +
  geom_col(aes(fill = variable), width = 100, alpha = 0.8) + ylim(0, 11000) + xlim(-2500,2500) +
  scale_fill_manual(values = c("#999999", "#D4D2D2", "#800026", "#FC4E2A", "#FD8D3C", "#FED976")) +
  ggtitle("set C") + xlab("Distance from origins (bp)") + ylab("Absolute contribution\n(number of mutations") + labs(fill = "COSMIC\nsignatures") +
  theme_classic() + theme(aspect.ratio=1.5)

pdf("./Rplot/Fig_S3C.pdf", width=12, height=12, useDingbats=FALSE)
ggarrange(setA.plot, setB.plot, setC.plot, 
          setA.SNS.core.plot, setB.SNS.core.plot, setC.SNS.core.plot, 
          setA.SNS.stochastic.plot, setB.SNS.stochastic.plot, setC.SNS.stochastic.plot, 
          ncol = 3, nrow = 3)
dev.off()

############################################################################
# Assess the contribution of polE and polD around constitutive origins

# Download individual signatures SBS10a,b (POLE) and SBS10c,d(POLD) from the COSMIC signature websites:
# https://cancer.sanger.ac.uk/signatures/sbs/
# and fit their contribution to constitutive origins

SBS10a <- read.csv("./Figure_2/v3.2_SBS10a_PROFILE.txt", sep = "\t") %>% dplyr::select(tri = X, SBS10a = SBS10a_GRCh38)
SBS10b <- read.csv("./Figure_2/v3.2_SBS10b_PROFILE.txt", sep = "\t") %>% dplyr::select(tri = X, SBS10b = SBS10b_GRCh38)
SBS10c <- read.csv("./Figure_2/v3.2_SBS10c_PROFILE_cxbZBIy.txt", sep = "\t") %>% dplyr::select(tri = X, SBS10c = SBS10c_GRCh38)
SBS10d <- read.csv("./Figure_2/v3.2_SBS10d_PROFILE_IJAo3WV.txt", sep = "\t") %>% dplyr::select(tri = X, SBS10d = SBS10d_GRCh38)
SBS.pol <- SBS10a %>% full_join(SBS10b) %>% full_join(SBS10c) %>% full_join(SBS10d) %>% dplyr::select(-tri)

Pol.fit.IniSeq <- fit_to_signatures(mut.mat.2, as.matrix(SBS.pol))
Pol.fit.IniSeq <- as.data.frame(Pol.fit.IniSeq$contribution)
Pol.fit.IniSeq <- as.data.frame(t(Pol.fit.IniSeq)) %>% dplyr::select_if(colSums(.) != 0) %>% add_rownames(var = "dist")
Pol.fit.IniSeq$dist <- as.numeric(as.character(Pol.fit.IniSeq$dist))
Pol.fit.IniSeq <- Pol.fit.IniSeq %>% arrange(dist)

# Compute mean of exposure outside origins and normalise data set

SBS10b.mean <- mean(Pol.fit.IniSeq$SBS10b[c(1:20,180:200)])
SBS10c.mean <- mean(Pol.fit.IniSeq$SBS10c[c(1:20,180:200)])
Pol.fit.IniSeq <- Pol.fit.IniSeq %>% mutate(SBS10b.norm = SBS10b/SBS10b.mean, SBS10c.norm = SBS10c/SBS10c.mean)

# Fit exposures

# SBS10b
# Fit using a sigmoidal on half of the profile
Sig.SBS10b.pos <- Pol.fit.IniSeq %>% filter(dist > 0)
model.Sig.SBS10b.pos <- drm(SBS10b.norm ~ log(dist), fct = G.4(), data = Sig.SBS10b.pos)  # G.4 : Gompertz curve
# Reconstruct profile
model.Sig.SBS10b <- c(rev(predict(model.Sig.SBS10b.pos)), predict(model.Sig.SBS10b.pos))
# Check fit
plot(Pol.fit.IniSeq$dist, Pol.fit.IniSeq$SBS10b.norm)
points(Pol.fit.IniSeq$dist, model.Sig.SBS10b, col = "red")

# SBS10c
# Fit using a sigmoidal on half of the profile
Sig.SBS10c.pos <- Pol.fit.IniSeq %>% filter(dist > 0)
model.Sig.SBS10c.pos <- drm(SBS10c.norm ~ log(dist), fct = G.4(), data = Sig.SBS10c.pos)  # G.4 : Gompertz curve
# Reconstruct profile
model.Sig.SBS10c <- c(rev(predict(model.Sig.SBS10c.pos)), predict(model.Sig.SBS10c.pos))
# Check fit
plot(Pol.fit.IniSeq$dist, Pol.fit.IniSeq$SBS10c.norm)
points(Pol.fit.IniSeq$dist, model.Sig.SBS10c, col = "red")

# Add fits to data

Pol.fit.IniSeq <- Pol.fit.IniSeq %>% mutate(SBS10b.fit = model.Sig.SBS10b, SBS10c.fit = model.Sig.SBS10c)
# Set negative values to 0
Pol.fit.IniSeq$SBS10b.fit[Pol.fit.IniSeq$SBS10b.fit<0] <- 0

# Plot

pdf("./Rplot/Fig_3D.pdf", width=8, height=4, useDingbats=FALSE)
coeff <- SBS10c.mean/SBS10b.mean
Pol.fit.IniSeq %>% ggplot(aes(x=dist)) + xlim(-5000,5000) +
  geom_point(aes(y=SBS10b.norm*SBS10b.mean), col = "#009E73", size = 1) + 
  geom_point(aes(y=SBS10c.norm*SBS10b.mean), col = "#4E84C4", size = 1) +
  geom_line(aes(y=model.Sig.SBS10b*SBS10b.mean), col = "#009E73", size = 1 ) + 
  geom_line(aes(y=model.Sig.SBS10c*SBS10b.mean), col = "#4E84C4", size = 1) + 
  scale_y_continuous(
    # Features of the first axis
    name = "SBS10b - PolE",
    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*coeff, name="SBS10c - PolD"),
    limits = c(-1,100)
  ) +
  ggtitle("Absolute contribution (number of mutations)") + xlab("Distance from origins (nt)") +
  theme_bw() + theme(aspect.ratio=1)
dev.off()

############################################################################
############################################################################
###             Distribution of DNA damage and repair factors            ###
###                      at constitutive origins                         ###
############################################################################
############################################################################

############################################################################
# Format dataset
# IniSeq origins positions are computed on hg38, hg19 and hg18 genome build
# for compatibility with available genomic data

IniSeq.ori <- read.table("./Dataset/Final_HL_enrchiment_of_2_for_3h_All_042019.gDNA_rm_dup_single_nt_position_by_100_bp_windows.bed", header = F)
colnames(IniSeq.ori) <- c("chr", "start", "end", "EFF")
# Compute positions of origins
IniSeq.peak.pos.hg38 <- IniSeq.ori %>% mutate(chr = paste0("chr", chr)) %>% mutate(center = round((start+end)/2)) %>% 
  mutate(start=as.integer(center), end=as.integer(center+1)) %>% dplyr::select(chr, start, end, EFF)
# Save hg38 coordinate
write.table(IniSeq.peak.pos.hg38, "./Figure_2/Repair/IniSeq.center.hg38.bed", sep="\t", col.names = F, row.names = F, quote = F)
# Create Grange object
IniSeq.hg38.gr <- makeGRangesFromDataFrame(IniSeq.peak.pos.hg38, keep.extra.columns=TRUE)
# liftover with rtracklayer
# download Liftover chain from UCSC
# http://hgdownload.cse.ucsc.edu/goldenpath/hg38/liftOver/
# https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/
# liftover to hg19 with rtracklayer
chain.hg38tohg19 <- import.chain("./Figure_2/Repair/hg38ToHg19.over.chain")
IniSeq.hg19.gr <- liftOver(IniSeq.hg38.gr, chain.hg38tohg19)
IniSeq.peak.pos.hg19 <- as.data.frame(IniSeq.hg19.gr@unlistData) %>% dplyr::select(chr = seqnames, start, end, EFF)
# Save hg19 coordinate
write.table(IniSeq.peak.pos.hg19, "./Figure_2/Repair/IniSeq.center.hg19.bed", sep="\t", col.names = F, row.names = F, quote = F)
# liftover to hg18 with rtracklayer
chain.hg19tohg18 <- import.chain("./Figure_2/Repair/hg19ToHg18.over.chain")
IniSeq.hg18.gr <- liftOver(IniSeq.hg19.gr, chain.hg19tohg18)
IniSeq.peak.pos.hg18 <- as.data.frame(IniSeq.hg18.gr@unlistData) %>% dplyr::select(chr = seqnames, start, end, EFF)
# Save hg18 coordinate
write.table(IniSeq.peak.pos.hg18, "./Figure_2/Repair/IniSeq.center.hg18.bed", sep="\t", col.names = F, row.names = F, quote = F)

# Prepare bed files reporting 10kb domains around origins

IniSeq.10kb.hg38 <- IniSeq.peak.pos.hg38  %>% mutate(end = as.integer(start + 10000), start = as.integer(start - 10000))
IniSeq.10kb.hg19 <- IniSeq.peak.pos.hg19  %>% mutate(end = as.integer(start + 10000), start = as.integer(start - 10000))
IniSeq.10kb.hg18 <- IniSeq.peak.pos.hg18  %>% mutate(end = as.integer(start + 10000), start = as.integer(start - 10000))
write.table(IniSeq.10kb.hg38, "./Figure_2/Repair/IniSeq.10kb.hg38.bed", sep="\t", col.names = F, row.names = F, quote = F)
write.table(IniSeq.10kb.hg19, "./Figure_2/Repair/IniSeq.10kb.hg19.bed", sep="\t", col.names = F, row.names = F, quote = F)
write.table(IniSeq.10kb.hg18, "./Figure_2/Repair/IniSeq.10kb.hg18.bed", sep="\t", col.names = F, row.names = F, quote = F)

# Split origins 10 kb domains into windows of 100 bp
# Using bedtools

# sort -k1,1 -k2,2n IniSeq.10kb.hg38.bed > IniSeq.10kb.hg38.sorted.bed
# bedtools makewindows -b IniSeq.10kb.hg38.sorted.bed -w 100 -i srcwinnum > IniSeq.10kb.hg38.100bp.windows.bed
# sort -k1,1 -k2,2n IniSeq.10kb.hg38.100bp.windows.bed > IniSeq.10kb.hg38.100bp.windows.sorted.bed
# 
# sort -k1,1 -k2,2n IniSeq.10kb.hg19.bed > IniSeq.10kb.hg19.sorted.bed
# bedtools makewindows -b IniSeq.10kb.hg19.sorted.bed -w 100 -i srcwinnum > IniSeq.10kb.hg19.100bp.windows.bed
# sort -k1,1 -k2,2n IniSeq.10kb.hg19.100bp.windows.bed > IniSeq.10kb.hg19.100bp.windows.sorted.bed
# 
# sort -k1,1 -k2,2n IniSeq.10kb.hg18.bed > IniSeq.10kb.hg18.sorted.bed
# bedtools makewindows -b IniSeq.10kb.hg18.sorted.bed -w 100 -i srcwinnum > IniSeq.10kb.hg18.100bp.windows.bed
# sort -k1,1 -k2,2n IniSeq.10kb.hg18.100bp.windows.bed > IniSeq.10kb.hg18.100bp.windows.sorted.bed

# Compute efficiency bins

quantile(IniSeq.ori$EFF, c(0.333,0.666))
#   33.3% 66.6% 
#   0.8418 0.9065 

############################################################################
# Analysis of TOP2 distribution

# Datasets #1
# GSM4205700 - MCF7 - hg19 from Martínez-García et al. PLoS Comput. Biol. 2021 17(1):e1007814 

# Convert bigwig to bedgraph
# bigWigToBedGraph ./TOP2/GSM4205700_MCF7_R1_TOP2B.bigWig ./TOP2/GSM4205700.TOP2B.MCF7.hg19.bedgraph

# Compute mean coverage for each IniSeq windows using bedGraph inputs
# bedtools map -a IniSeq.10kb.hg19.100bp.windows.sorted.bed -b ./TOP2/GSM4205700.TOP2B.MCF7.hg19.bedgraph -c 4 -o mean -null 0 > ./TOP2/GSM4205700.TOP2B.MCF7.hg19.IniSeq.cov.bedgraph

# Open and format results

TOP2B.MCF7.cov <- read.table("./Figure_2/Repair/TOP2/GSM4205700.TOP2B.MCF7.hg19.IniSeq.cov.bedgraph", header = F)
colnames(TOP2B.MCF7.cov) <- c("chr", "start", "end", "INFO", "COV")
TOP2B.MCF7.cov <- TOP2B.MCF7.cov %>% separate(INFO, c("EFF", "bin"), sep = "_") %>% mutate(dist = (as.numeric(bin)-100)*100)

# Compute fold enrichment in coverage
TOP2B.MCF7.cov.plot <- TOP2B.MCF7.cov %>% group_by(dist) %>% summarise(COV.mean = mean(COV)) %>% mutate(COV.fold = log2(COV.mean/mean(COV.mean[c(1:20,180:200)]))) %>% mutate(cell = "TOP2B in MCF7 cells")

# Bin by origins efficiency
TOP2B.MCF7.low.cov <- TOP2B.MCF7.cov %>% filter(EFF < 0.8418) %>% group_by(dist) %>% summarise(COV.mean = mean(COV)) %>% mutate(COV.fold = log2(COV.mean/mean(COV.mean[c(1:20,180:200)]))) %>% mutate(efficiency = "low")  
TOP2B.MCF7.medium.cov <- TOP2B.MCF7.cov %>% filter(EFF >= 0.8418 & EFF < 0.9065) %>% group_by(dist) %>% summarise(COV.mean = mean(COV)) %>% mutate(COV.fold = log2(COV.mean/mean(COV.mean[c(1:20,180:200)]))) %>% mutate(efficiency = "medium")  
TOP2B.MCF7.high.cov <- TOP2B.MCF7.cov %>% filter(EFF >= 0.9065) %>% group_by(dist) %>% summarise(COV.mean = mean(COV)) %>% mutate(COV.fold = log2(COV.mean/mean(COV.mean[c(1:20,180:200)]))) %>% mutate(efficiency = "high")  
TOP2B.MCF7.eff.cov.plot <- rbind(TOP2B.MCF7.low.cov, TOP2B.MCF7.medium.cov, TOP2B.MCF7.high.cov)

# Plot coverage

TOP2B.EFF.plot <- TOP2B.MCF7.eff.cov.plot  %>% ggplot(aes(x = dist, y = COV.fold, color = efficiency)) +
  geom_line(size = 0.5) + xlim(-5000,5000) +
  scale_color_manual(values=c("#FC4E07", "#0072B2", "#CC79A7")) +
  xlab("Distance from origin (bp)") + ylab("TOP2B enrichment (log2 FC)") + ggtitle("TOP2B enrichment in MCF7 cells around IniSeq origins") +
  theme_bw() + theme(aspect.ratio=1)

# Datasets #2
# GSM2442946 - MCF10 - hg18 from hg18 from Dellino et al. Nature Genetics 2019

# Convert bigwig to bedgraph
# bigWigToBedGraph ./TOP2/GSM2442946_MCF10AsiSIER_NT_TOP2B_chip.bw ./TOP2/GSM2442946.TOP2B.MCF10.hg18.bedgraph

# Compute mean coverage for each IniSeq windows using bedGraph inputs
# bedtools map -a IniSeq.10kb.hg18.100bp.windows.sorted.bed -b ./TOP2/GSM2442946.TOP2B.MCF10.hg18.bedgraph -c 4 -o mean -null 0 > ./TOP2/GSM2442946.TOP2B.MCF10.hg18.IniSeq.cov.bedgraph

# Open and format results

TOP2B.MCF10A.cov <- read.table("./Figure_2/Repair/TOP2/GSM2442946.TOP2B.MCF10.hg18.IniSeq.cov.bedgraph", header = F)
colnames(TOP2B.MCF10A.cov) <- c("chr", "start", "end", "INFO", "COV")
TOP2B.MCF10A.cov <- TOP2B.MCF10A.cov %>% separate(INFO, c("EFF", "bin"), sep = "_") %>% mutate(dist = (as.numeric(bin)-100)*100)

# Compute fold enrichment in coverage
TOP2B.MCF10A.cov.plot <- TOP2B.MCF10A.cov %>% group_by(dist) %>% summarise(COV.mean = mean(COV)) %>% mutate(COV.fold = log2(COV.mean/mean(COV.mean[c(1:20,180:200)]))) %>% mutate(cell = "TOP2B in MCF10A cells")

# Average TOP2B signal over all cell lines

TOP2B.MCF7 <- TOP2B.MCF7.cov.plot %>% dplyr::select(dist, COV.MCF7 = COV.fold)
TOP2B.MCF10A <- TOP2B.MCF10A.cov.plot %>% dplyr::select(dist, COV.MCF10A = COV.fold)
TOP2B <- TOP2B.MCF7 %>% right_join(TOP2B.MCF10A) %>% mutate(COV.fold = (COV.MCF7+COV.MCF10A)/2, Factor = "TOP2B") %>% dplyr::select(dist, COV.fold, Factor)

############################################################################
# Analysis of TOP1 distribution

# Datasets #1
# GSM3182665 - MCF10A - hg18 from Dellino et al. Nature Genetics 2019

# Convert bigwig to bedgraph
# bigWigToBedGraph ./TOP1/GSM3182665_MCF10AsiSIER_NT_TOP1_chip.bw ./TOP1/GSM3182665.TOP1.MCF10A.hg18.bedgraph

# Compute mean coverage for each IniSeq windows using bedGraph inputs
# bedtools map -a IniSeq.10kb.hg18.100bp.windows.sorted.bed -b ./TOP1/GSM3182665.TOP1.MCF10A.hg18.bedgraph -c 4 -o mean -null 0 > ./TOP1/GSM3182665.TOP1.MCF10A.hg18.IniSeq.cov.bedgraph

# Open and format results

TOP1.MCF10A.cov <- read.table("./Figure_2/Repair/TOP1/GSM3182665.TOP1.MCF10A.hg18.IniSeq.cov.bedgraph", header = F)
colnames(TOP1.MCF10A.cov) <- c("chr", "start", "end", "INFO", "COV")
TOP1.MCF10A.cov <- TOP1.MCF10A.cov %>% separate(INFO, c("EFF", "bin"), sep = "_") %>% mutate(dist = (as.numeric(bin)-100)*100)

# Compute fold enrichment in coverage
TOP1.MCF10A.cov.plot <- TOP1.MCF10A.cov %>% group_by(dist) %>% summarise(COV.mean = mean(COV)) %>% mutate(COV.fold = log2(COV.mean/mean(COV.mean[c(1:20,180:200)]))) %>% mutate(cell = "TOP1 in MCF10A cells")

# Bin by origins efficiency
TOP1.MCF10A.low.cov <- TOP1.MCF10A.cov %>% filter(EFF < 0.8418) %>% group_by(dist) %>% summarise(COV.mean = mean(COV)) %>% mutate(COV.fold = log2(COV.mean/mean(COV.mean[c(1:20,180:200)]))) %>% mutate(efficiency = "low")  
TOP1.MCF10A.medium.cov <- TOP1.MCF10A.cov %>% filter(EFF >= 0.8418 & EFF < 0.9065) %>% group_by(dist) %>% summarise(COV.mean = mean(COV)) %>% mutate(COV.fold = log2(COV.mean/mean(COV.mean[c(1:20,180:200)]))) %>% mutate(efficiency = "medium")  
TOP1.MCF10A.high.cov <- TOP1.MCF10A.cov %>% filter(EFF >= 0.9065) %>% group_by(dist) %>% summarise(COV.mean = mean(COV)) %>% mutate(COV.fold = log2(COV.mean/mean(COV.mean[c(1:20,180:200)]))) %>% mutate(efficiency = "high")  
TOP1.MCF10A.eff.cov.plot <- rbind(TOP1.MCF10A.low.cov, TOP1.MCF10A.medium.cov, TOP1.MCF10A.high.cov)

# Plot coverage

TOP1.EFF.plot <- TOP1.MCF10A.eff.cov.plot  %>% ggplot(aes(x = dist, y = COV.fold, color = efficiency)) +
  geom_line(size = 0.5) + xlim(-5000,5000) +
  scale_color_manual(values=c("#FC4E07", "#0072B2", "#CC79A7")) +
  xlab("Distance from origin (bp)") + ylab("TOP1 enrichment (log2 FC)") + ggtitle("TOP1 enrichment in MCF10A cells around IniSeq origins") +
  theme_bw() + theme(aspect.ratio=1)

# Datasets #2
# GSM2058666 - HCT116 - hg19 from Baranello et al. Cell 2016 Apr 7;165(2):357-71

# Convert Wig to bigWig
# wigToBigWig ./TOP1/GSM2058666_HCT116_wt_Top1_ChIPSeq_ctrl_normalized.profile.wig hg19.chrom.sizes ./TOP1/GSM2058666_TOP1_HCT116_hg19.bw
# Delete Wig file

# Convert bigwig to bedgraph
# bigWigToBedGraph ./TOP1/GSM2058666_TOP1_HCT116_hg19.bw ./TOP1/GSM2058666_TOP1_HCT116_hg19.bedgraph

# Compute mean coverage for each IniSeq windows using bedGraph inputs
# bedtools map -a IniSeq.10kb.hg19.100bp.windows.sorted.bed -b ./TOP1/GSM2058666_TOP1_HCT116_hg19.bedgraph -c 4 -o mean -null 0 > ./TOP1/GSM2058666_TOP1_HCT116_hg19.IniSeq.cov.bedgraph

# Open and format results

TOP1.HCT116.cov <- read.table("./Figure_2/Repair/TOP1/GSM2058666_TOP1_HCT116_hg19.IniSeq.cov.bedgraph", header = F)
colnames(TOP1.HCT116.cov) <- c("chr", "start", "end", "INFO", "COV")
TOP1.HCT116.cov <- TOP1.HCT116.cov %>% separate(INFO, c("EFF", "bin"), sep = "_") %>% mutate(dist = (as.numeric(bin)-100)*100)

# Compute fold enrichment in coverage
TOP1.HCT116.cov.plot <- TOP1.HCT116.cov %>% group_by(dist) %>% summarise(COV.mean = mean(COV)) %>% mutate(COV.fold = log2(COV.mean/mean(COV.mean[c(1:20,180:200)]))) %>% mutate(cell = "TOP1 in HCT116 cells")

# Datasets #3
# GSM1385717 - HCT116 - hg19 from Baranello et al. Cell 2016 Apr 7;165(2):357-71
# Top-seq in the presence of CPT (Campthotecin) for mapping catalytically engaged TOP1

# Convert Wig to bigWig
# wigToBigWig ./TOP1/GSM1385718_HCT116_wt_+CPT_Top1_Seq_rep2_allreads.wig hg19.chrom.sizes ./TOP1/GSM1385717_TOP1_seq_2_HCT116_hg19.bw
# Delete Wig file

# Convert bigwig to bedgraph
# bigWigToBedGraph ./TOP1/GSM1385717_TOP1_seq_2_HCT116_hg19.bw ./TOP1/GSM1385717_TOP1_seq_2_HCT116_hg19.bedgraph

# Compute mean coverage for each IniSeq windows using bedGraph inputs
# bedtools map -a IniSeq.10kb.hg19.100bp.windows.sorted.bed -b ./TOP1/GSM1385717_TOP1_seq_2_HCT116_hg19.bedgraph -c 4 -o mean -null 0 > ./TOP1/GSM1385717_TOP1_seq_2_HCT116_hg19.IniSeq.cov.bedgraph

# Open and format results

TOP1.seq.2.HCT116.cov <- read.table("./Figure_2/Repair/TOP1/GSM1385717_TOP1_seq_2_HCT116_hg19.IniSeq.cov.bedgraph", header = F)
colnames(TOP1.seq.2.HCT116.cov) <- c("chr", "start", "end", "INFO", "COV")
TOP1.seq.2.HCT116.cov <- TOP1.seq.2.HCT116.cov %>% separate(INFO, c("EFF", "bin"), sep = "_") %>% mutate(dist = (as.numeric(bin)-100)*100)

# Compute fold enrichment in coverage
TOP1.seq.2.HCT116.cov.plot <- TOP1.seq.2.HCT116.cov %>% group_by(dist) %>% summarise(COV.mean = mean(COV)) %>% mutate(COV.fold = log2(COV.mean/mean(COV.mean[c(1:20,180:200)]))) %>% mutate(cell = "TOP1 seq in HCT116 cells")

# Average TOP1 signal over all cell lines

TOP1.MCF10A <- TOP1.MCF10A.cov.plot %>% dplyr::select(dist, COV.MCF10A = COV.fold)
TOP1.HCT116 <- TOP1.HCT116.cov.plot %>% dplyr::select(dist, COV.HCT116 = COV.fold)
TOP1.seq.HCT116 <- TOP1.seq.2.HCT116.cov.plot %>% dplyr::select(dist, COV.HCT116.2 = COV.fold)
TOP1 <- TOP1.MCF10A %>% right_join(TOP1.HCT116) %>% right_join(TOP1.seq.HCT116) %>% mutate(COV.fold = (COV.MCF10A+COV.HCT116+COV.HCT116.2)/3, Factor = "TOP1") %>%
  dplyr::select(dist, COV.fold, Factor)

############################################################################
# Analysis of dsDNA breaks by BLISS-seq

# Datasets #1 
# GSM3444984 - TK6 - hg19 from Gothe et al. Mol Cell 2019 Jul 25;75(2):267-283.e12

# Convert bigwig to bedgraph
# bigWigToBedGraph ./BLISSseq/GSM3444984_TK6_DMSO_rep1.bw ./BLISSseq/GSM3444984.BLISSseq.TK6.hg19.bedgraph

# Compute mean coverage for each IniSeq windows using bedGraph inputs
# bedtools map -a IniSeq.10kb.hg19.100bp.windows.sorted.bed -b ./BLISSseq/GSM3444984.BLISSseq.TK6.hg19.bedgraph -c 4 -o mean -null 0 > ./BLISSseq/GSM3444984.BLISSseq.TK6.hg19.IniSeq.cov.bedgraph

# Open and format results

BLISSseq.TK6.cov <- read.table("./Figure_2/Repair/BLISSseq/GSM3444984.BLISSseq.TK6.hg19.IniSeq.cov.bedgraph", header = F)
colnames(BLISSseq.TK6.cov) <- c("chr", "start", "end", "INFO", "COV")
BLISSseq.TK6.cov <- BLISSseq.TK6.cov %>% separate(INFO, c("EFF", "bin"), sep = "_") %>% mutate(dist = (as.numeric(bin)-100)*100)

# Compute fold enrichment in coverage
BLISSseq.TK6.cov.plot <- BLISSseq.TK6.cov %>% group_by(dist) %>% summarise(COV.mean = mean(COV)) %>% mutate(COV.fold = log2(COV.mean/mean(COV.mean[c(1:20,180:200)]))) %>% mutate(cell = "TOP1 in TK6 cells")

# Datasets #2 
# GSM3444988 - K562 - hg19 from Gothe et al. Mol Cell 2019 Jul 25;75(2):267-283.e12

# Convert bigwig to bedgraph
# bigWigToBedGraph ./BLISSseq/GSM3444988_K562_DMSO.bw ./BLISSseq/GSM3444988.BLISSseq.K562.hg19.bedgraph

# Compute mean coverage for each IniSeq windows using bedGraph inputs
# bedtools map -a IniSeq.10kb.hg19.100bp.windows.sorted.bed -b ./BLISSseq/GSM3444988.BLISSseq.K562.hg19.bedgraph -c 4 -o mean -null 0 > ./BLISSseq/GSM3444988.BLISSseq.K562.hg19.IniSeq.cov.bedgraph

# Open and format results

BLISSseq.K562.cov <- read.table("./Figure_2/Repair/BLISSseq/GSM3444988.BLISSseq.K562.hg19.IniSeq.cov.bedgraph", header = F)
colnames(BLISSseq.K562.cov) <- c("chr", "start", "end", "INFO", "COV")
BLISSseq.K562.cov <- BLISSseq.K562.cov %>% separate(INFO, c("EFF", "bin"), sep = "_") %>% mutate(dist = (as.numeric(bin)-100)*100)

# Compute fold enrichment in coverage
BLISSseq.K562.cov.plot <- BLISSseq.K562.cov %>% group_by(dist) %>% summarise(COV.mean = mean(COV)) %>% mutate(COV.fold = log2(COV.mean/mean(COV.mean[c(1:20,180:200)]))) %>% mutate(cell = "TOP1 in K562 cells")

# Datasets #3 
# GSM3687236 - CD34 - hg19 from Gothe et al. Mol Cell 2019 Jul 25;75(2):267-283.e12

# Convert bigwig to bedgraph
# bigWigToBedGraph ./BLISSseq/GSM3687236_CD34_DMSO_rep1.bw ./BLISSseq/GSM3687236.BLISSseq.K562.hg19.bedgraph

# Compute mean coverage for each IniSeq windows using bedGraph inputs
# bedtools map -a IniSeq.10kb.hg19.100bp.windows.sorted.bed -b ./BLISSseq/GSM3687236.BLISSseq.K562.hg19.bedgraph -c 4 -o mean -null 0 > ./BLISSseq/GSM3687236.BLISSseq.K562.hg19.IniSeq.cov.bedgraph

# Open and format results

BLISSseq.CD34.cov <- read.table("./Figure_2/Repair/BLISSseq/GSM3687236.BLISSseq.K562.hg19.IniSeq.cov.bedgraph", header = F)
colnames(BLISSseq.CD34.cov) <- c("chr", "start", "end", "INFO", "COV")
BLISSseq.CD34.cov <- BLISSseq.CD34.cov %>% separate(INFO, c("EFF", "bin"), sep = "_") %>% mutate(dist = (as.numeric(bin)-100)*100)

# Compute fold enrichment in coverage
BLISSseq.CD34.cov.plot <- BLISSseq.CD34.cov %>% group_by(dist) %>% summarise(COV.mean = mean(COV)) %>% mutate(COV.fold = log2(COV.mean/mean(COV.mean[c(1:20,180:200)]))) %>% mutate(cell = "TOP1 in CD34 cells")

# Prepare plot averaging signal from all cell lines

BLISSseq.TK6 <- BLISSseq.TK6.cov.plot %>% dplyr::select(dist, COV.TK6 = COV.fold)
BLISSseq.K562 <- BLISSseq.K562.cov.plot %>% dplyr::select(dist, COV.K562 = COV.fold)
BLISSseq.CD34 <- BLISSseq.CD34.cov.plot %>% dplyr::select(dist, COV.CD34 = COV.fold)
BLISSseq <- BLISSseq.TK6 %>% right_join(BLISSseq.K562) %>% right_join(BLISSseq.CD34) %>% mutate(COV.fold = (COV.TK6+COV.K562+COV.CD34)/3, Factor = "BLISS-seq") %>% dplyr::select(dist, COV.fold, Factor)

############################################################################
# Anaysis of BRCA1 distribution

# Datasets #1
# GSM1340576 - MCF10A - hg18 from Gardini et al. EMBO J 2014 16 33 890-905

# Convert bigwig to bedgraph
# bigWigToBedGraph ./BRCA1/GSM1340576_D9.bw ./BRCA1/GSM1340576.BRCA1.MCF10A.hg18.bedgraph

# Compute mean coverage for each IniSeq windows using bedGraph inputs
# bedtools map -a IniSeq.10kb.hg18.100bp.windows.sorted.bed -b ./BRCA1/GSM1340576.BRCA1.MCF10A.hg18.bedgraph -c 4 -o mean -null 0 > ./BRCA1/GSM1340576.BRCA1.MCF10A.hg18.IniSeq.cov.bedgraph

# Open and format results

BRCA1.MCF10.cov <- read.table("./Figure_2/Repair/BRCA1/GSM1340576.BRCA1.MCF10A.hg18.IniSeq.cov.bedgraph", header = F)
colnames(BRCA1.MCF10.cov) <- c("chr", "start", "end", "INFO", "COV")
BRCA1.MCF10.cov <- BRCA1.MCF10.cov %>% separate(INFO, c("EFF", "bin"), sep = "_") %>% mutate(dist = (as.numeric(bin)-100)*100)

# Compute fold enrichment in coverage
BRCA1.MCF10.cov.plot <- BRCA1.MCF10.cov %>% group_by(dist) %>% summarise(COV.mean = mean(COV)) %>% mutate(COV.fold = log2(COV.mean/mean(COV.mean[c(1:20,180:200)]))) %>% mutate(cell = "MCF10A")

# Bin by origins efficiency
BRCA1.MCF10.low.cov <- BRCA1.MCF10.cov %>% filter(EFF < 0.8418) %>% group_by(dist) %>% summarise(COV.mean = mean(COV)) %>% mutate(COV.fold = log2(COV.mean/mean(COV.mean[c(1:20,180:200)]))) %>% mutate(efficiency = "low")  
BRCA1.MCF10.medium.cov <- BRCA1.MCF10.cov %>% filter(EFF >= 0.8418 & EFF < 0.9065) %>% group_by(dist) %>% summarise(COV.mean = mean(COV)) %>% mutate(COV.fold = log2(COV.mean/mean(COV.mean[c(1:20,180:200)]))) %>% mutate(efficiency = "medium")  
BRCA1.MCF10.high.cov <- BRCA1.MCF10.cov %>% filter(EFF >= 0.9065) %>% group_by(dist) %>% summarise(COV.mean = mean(COV)) %>% mutate(COV.fold = log2(COV.mean/mean(COV.mean[c(1:20,180:200)]))) %>% mutate(efficiency = "high")  
BRCA1.MCF10.eff.cov.plot <- rbind(BRCA1.MCF10.low.cov, BRCA1.MCF10.medium.cov, BRCA1.MCF10.high.cov)

# Plot coverage

BRCA1.EFF.plot <- BRCA1.MCF10.eff.cov.plot  %>% ggplot(aes(x = dist, y = COV.fold, color = efficiency)) +
  geom_line(size = 0.5) +
  scale_color_manual(values=c("#FC4E07", "#0072B2", "#CC79A7")) + xlim(-5000,5000) +
  xlab("Distance from origin (bp)") + ylab("BRCA1 coverage (log2 FC)") + ggtitle("BRCA1 coverage in MCF10 cells around IniSeq origins") +
  theme_bw() + theme(aspect.ratio=1)

# Datasets #2
# GSM935517 - H1-hESC - hg19 from ENCODE Project Consortium Nature 2012  489 57-74

# Convert bigwig to bedgraph
# bigWigToBedGraph ./BRCA1/GSM935517_hg19_wgEncodeSydhTfbsH1hescBrca1IggrabSig.bigWig ./BRCA1/GSM935517.BRCA1.H1.hg19.bedgraph

# Compute mean coverage for each IniSeq windows using bedGraph inputs
# bedtools map -a IniSeq.10kb.hg19.100bp.windows.sorted.bed -b ./BRCA1/GSM935517.BRCA1.H1.hg19.bedgraph -c 4 -o mean -null 0 > ./BRCA1/GSM935517.BRCA1.H1.hg19.IniSeq.cov.bedgraph

# Open and format results

BRCA1.H1.cov <- read.table("./Figure_2/Repair/BRCA1/GSM935517.BRCA1.H1.hg19.IniSeq.cov.bedgraph", header = F)
colnames(BRCA1.H1.cov) <- c("chr", "start", "end", "INFO", "COV")
BRCA1.H1.cov <- BRCA1.H1.cov %>% separate(INFO, c("EFF", "bin"), sep = "_") %>% mutate(dist = (as.numeric(bin)-100)*100)

# Compute fold enrichment in coverage
BRCA1.H1.cov.plot <- BRCA1.H1.cov %>% group_by(dist) %>% summarise(COV.mean = mean(COV)) %>% mutate(COV.fold = log2(COV.mean/mean(COV.mean[c(1:20,180:200)]))) %>% mutate(cell = "H1")

# Datasets #3
# GSM935552 - HeLa - hg19 from ENCODE Project Consortium Nature 2012  489 57-74

# Convert bigwig to bedgraph
# bigWigToBedGraph ./BRCA1/GSM935552_hg19_wgEncodeSydhTfbsHelas3Brca1a300IggrabSig.bigWig ./BRCA1/GSM935552.BRCA1.HeLa.hg19.bedgraph

# Compute mean coverage for each IniSeq windows using bedGraph inputs
# bedtools map -a IniSeq.10kb.hg19.100bp.windows.sorted.bed -b ./BRCA1/GSM935552.BRCA1.HeLa.hg19.bedgraph -c 4 -o mean -null 0 > ./BRCA1/GSM935552.BRCA1.HeLa.hg19.IniSeq.cov.bedgraph

# Open and format results

BRCA1.HeLa.cov <- read.table("./Figure_2/Repair/BRCA1/GSM935552.BRCA1.HeLa.hg19.IniSeq.cov.bedgraph", header = F)
colnames(BRCA1.HeLa.cov) <- c("chr", "start", "end", "INFO", "COV")
BRCA1.HeLa.cov <- BRCA1.HeLa.cov %>% separate(INFO, c("EFF", "bin"), sep = "_") %>% mutate(dist = (as.numeric(bin)-100)*100)

# Compute fold enrichment in coverage
BRCA1.HeLa.cov.plot <- BRCA1.HeLa.cov %>% group_by(dist) %>% summarise(COV.mean = mean(COV)) %>% mutate(COV.fold = log2(COV.mean/mean(COV.mean[c(1:20,180:200)]))) %>% mutate(cell = "HeLa")

# Datasets #4
# GSM935377 - GM12878 - hg19 from ENCODE Project Consortium Nature 2012  489 57-74

# Convert bigwig to bedgraph
# bigWigToBedGraph ./BRCA1/GSM935377_hg19_wgEncodeSydhTfbsGm12878Brca1a300IggmusSig.bigWig ./BRCA1/GSM935377.BRCA1.GM12878.hg19.bedgraph

# Compute mean coverage for each IniSeq windows using bedGraph inputs
# bedtools map -a IniSeq.10kb.hg19.100bp.windows.sorted.bed -b ./BRCA1/GSM935377.BRCA1.GM12878.hg19.bedgraph -c 4 -o mean -null 0 > ./BRCA1/GSM935377.BRCA1.GM12878.hg19.IniSeq.cov.bedgraph

# Open and format results

BRCA1.GM12878.cov <- read.table("./Figure_2/Repair/BRCA1/GSM935377.BRCA1.GM12878.hg19.IniSeq.cov.bedgraph", header = F)
colnames(BRCA1.GM12878.cov) <- c("chr", "start", "end", "INFO", "COV")
BRCA1.GM12878.cov <- BRCA1.GM12878.cov %>% separate(INFO, c("EFF", "bin"), sep = "_") %>% mutate(dist = (as.numeric(bin)-100)*100)

# Compute fold enrichment in coverage
BRCA1.GM12878.cov.plot <- BRCA1.GM12878.cov %>% group_by(dist) %>% summarise(COV.mean = mean(COV)) %>% mutate(COV.fold = log2(COV.mean/mean(COV.mean[c(1:20,180:200)]))) %>% mutate(cell = "GM12878")

# Average BRCA1 signal over cell lines

BRCA1.MCF10 <- BRCA1.MCF10.cov.plot %>% dplyr::select(dist, COV.MCF10 = COV.fold)
BRCA1.H1 <- BRCA1.H1.cov.plot %>% dplyr::select(dist, COV.H1 = COV.fold)
BRCA1.HeLa <- BRCA1.HeLa.cov.plot %>% dplyr::select(dist, COV.HeLa = COV.fold)
BRCA1 <- BRCA1.GM12878.cov.plot %>% dplyr::select(dist, COV.GM12878 = COV.fold ) %>% right_join(BRCA1.MCF10) %>%
  right_join(BRCA1.H1) %>% right_join(BRCA1.HeLa) %>% mutate(COV.fold = (COV.MCF10+COV.H1+COV.HeLa+COV.GM12878)/4) %>% 
  mutate(Factor = "BRCA1") %>% dplyr::select(dist, COV.fold, Factor)

############################################################################
# Anaysis of XRCC proteins distribution

# Datasets #1 - XRCC5
# GSM3393608 - K562 - hg19 from ENCODE Project Consortium Nature 2012  489 57-74

# Convert bigwig to bedgraph
# bigWigToBedGraph ./XRCCs/GSM3393608_XRCC5_ChIP_K562_B7.pval.signal.bw ./XRCCs/GSM3393608.XRCC5.K562.hg19.bedgraph

# Compute mean coverage for each IniSeq windows using bedGraph inputs
# bedtools map -a IniSeq.10kb.hg19.100bp.windows.sorted.bed -b ./XRCCs/GSM3393608.XRCC5.K562.hg19.bedgraph -c 4 -o mean -null 0 > ./XRCCs/GSM3393608.XRCC5.K562.hg19.IniSeq.cov.bedgraph

# Open and format results

XRCC5.K562.cov <- read.table("./Figure_2/Repair/XRCCs/GSM3393608.XRCC5.K562.hg19.IniSeq.cov.bedgraph", header = F)
colnames(XRCC5.K562.cov) <- c("chr", "start", "end", "INFO", "COV")
XRCC5.K562.cov <- XRCC5.K562.cov %>% separate(INFO, c("EFF", "bin"), sep = "_") %>% mutate(dist = (as.numeric(bin)-100)*100)

# Compute fold enrichment in coverage
XRCC5.K562.cov.plot <- XRCC5.K562.cov %>% group_by(dist) %>% summarise(COV.mean = mean(COV)) %>% mutate(COV.fold = log2(COV.mean/mean(COV.mean[c(1:20,180:200)]))) %>% mutate(cell = "K562")

# Datasets #2 - XRCC5
# GSM3393560 - HepG2 - hg19 from Xiao et al Cell 2019 178 107-121

# Convert bigwig to bedgraph
# bigWigToBedGraph ./XRCCs/GSM3393560_XRCC5_ChIP_HepG2_B16.pval.signal.bw ./XRCCs/GSM3393560.XRCC5.HepG2.hg19.bedgraph

# Compute mean coverage for each IniSeq windows using bedGraph inputs
# bedtools map -a IniSeq.10kb.hg19.100bp.windows.sorted.bed -b ./XRCCs/GSM3393560.XRCC5.HepG2.hg19.bedgraph -c 4 -o mean -null 0 > ./XRCCs/GSM3393560.XRCC5.HepG2.hg19.IniSeq.cov.bedgraph

# Open and format results

XRRC5.HepG2.cov <- read.table("./Figure_2/Repair/XRCCs/GSM3393560.XRCC5.HepG2.hg19.IniSeq.cov.bedgraph", header = F)
colnames(XRRC5.HepG2.cov) <- c("chr", "start", "end", "INFO", "COV")
XRRC5.HepG2.cov <- XRRC5.HepG2.cov %>% separate(INFO, c("EFF", "bin"), sep = "_") %>% mutate(dist = (as.numeric(bin)-100)*100)

# Compute fold enrichment in coverage
XRRC5.HepG2.cov.plot <- XRRC5.HepG2.cov %>% group_by(dist) %>% summarise(COV.mean = mean(COV)) %>% mutate(COV.fold = log2(COV.mean/mean(COV.mean[c(1:20,180:200)]))) %>% mutate(cell = "HepG2")

# Average XRRC5 signal over cell lines

XRCC5.K562 <- XRCC5.K562.cov.plot %>% dplyr::select(dist, COV.K562 = COV.fold)
XRCC5 <- XRRC5.HepG2.cov.plot %>% dplyr::select(dist, COV.HepG2 = COV.fold ) %>% right_join(XRCC5.K562) %>%
  mutate(COV.fold = (COV.K562+COV.HepG2)/2, Factor = "XRCC5") %>% dplyr::select(dist, COV.fold, Factor)

############################################################################
# Anaysis of PARP1 distribution

# Datasets #1
# GSM3182666 - MCF10 - hg18 from Dellino et al. Nature Genetics 2019

# Convert bigwig to bedgraph
# bigWigToBedGraph ./PARP1/GSM3182666_MCF10AsiSIER_NT_PARP1_chip.bw ./PARP1/GSM3182666.PARP1.MCF10.hg18.bedgraph

# Compute mean coverage for each IniSeq windows using bedGraph inputs
# bedtools map -a IniSeq.10kb.hg18.100bp.windows.sorted.bed -b ./PARP1/GSM3182666.PARP1.MCF10.hg18.bedgraph -c 4 -o mean -null 0 > ./PARP1/GSM3182666.PARP1.MCF10.hg18.IniSeq.cov.bedgraph

# Open and format results

PARP1.MCF10.cov <- read.table("./Figure_2/Repair/PARP1/GSM3182666.PARP1.MCF10.hg18.IniSeq.cov.bedgraph", header = F)
colnames(PARP1.MCF10.cov) <- c("chr", "start", "end", "INFO", "COV")
PARP1.MCF10.cov <- PARP1.MCF10.cov %>% separate(INFO, c("EFF", "bin"), sep = "_") %>% mutate(dist = (as.numeric(bin)-100)*100)

# Compute fold enrichment in coverage
PARP1.MCF10.cov.plot <- PARP1.MCF10.cov %>% group_by(dist) %>% summarise(COV.mean = mean(COV)) %>% mutate(COV.fold = log2(COV.mean/mean(COV.mean[c(1:20,180:200)]))) %>% mutate(cell = "MCF10")

# Datasets #2 & #3
# GSE61916 contains PARP1 ChIP-seq in two different cell lines - hg19 from Nalabothula et al.PLoS One 2015;10(8):e0135410.
# GSM1517305 - MCF7
# GSM1517306 - MDA-MB231

# Merge wig files and generate bigwig file

# gzcat ./PARP1/GSE61916/GSE61916_RAW/GSM1517305/*.wig.gz > ./PARP1/GSM1517305.PARP1.MCF7.hg19.wig
# gzcat ./PARP1/GSE61916/GSE61916_RAW/GSM1517306/*.wig.gz > ./PARP1/GSM1517306.PARP1.MDA.MB231.hg19.wig

# Remove headers
# grep -v ^track ./PARP1/GSM1517305.PARP1.MCF7.hg19.wig > ./PARP1/GSM1517305.PARP1.MCF7.hg19.2.wig
# grep -v ^track ./PARP1/GSM1517306.PARP1.MDA.MB231.hg19.wig > ./PARP1/GSM1517306.PARP1.MDA.MB231.hg19.2.wig

# Convert wig to bigwig
# wigToBigWig -clip ./PARP1/GSM1517305.PARP1.MCF7.hg19.2.wig hg19.chrom.sizes ./PARP1/GSM1517305.PARP1.MCF7.hg19.bw
# wigToBigWig -clip ./PARP1/GSM1517306.PARP1.MDA.MB231.hg19.2.wig hg19.chrom.sizes ./PARP1/GSM1517306.PARP1.MDA.MB231.hg19.bw

# Convert bigwig to bedgraph
# bigWigToBedGraph ./PARP1/GSM1517305.PARP1.MCF7.hg19.bw ./PARP1/GSM1517305.PARP1.MCF7.hg19.bedgraph
# bigWigToBedGraph ./PARP1/GSM1517306.PARP1.MDA.MB231.hg19.bw ./PARP1/GSM1517306.PARP1.MDA.MB231.hg19.bedgraph

# Compute mean coverage for each IniSeq windows using bedGraph inputs
# bedtools map -a IniSeq.10kb.hg19.100bp.windows.sorted.bed -b ./PARP1/GSM1517305.PARP1.MCF7.hg19.bedgraph -c 4 -o mean -null 0 > ./PARP1/GSM1517305.PARP1.MCF7.hg19.IniSeq.cov.bedgraph
# bedtools map -a IniSeq.10kb.hg19.100bp.windows.sorted.bed -b ./PARP1/GSM1517306.PARP1.MDA.MB231.hg19.bedgraph -c 4 -o mean -null 0 > ./PARP1/GSM1517306.PARP1.MDA.MB231.hg19.IniSeq.cov.bedgraph

# Clean and remove intermediates files

# Open and format results

PARP1.MCF7.cov <- read.table("./Figure_2/Repair/PARP1/GSM1517305.PARP1.MCF7.hg19.IniSeq.cov.bedgraph", header = F)
colnames(PARP1.MCF7.cov) <- c("chr", "start", "end", "INFO", "COV")
PARP1.MCF7.cov <- PARP1.MCF7.cov %>% separate(INFO, c("EFF", "bin"), sep = "_") %>% mutate(dist = (as.numeric(bin)-100)*100)

# Compute fold enrichment in coverage
PARP1.MCF7.cov.plot <- PARP1.MCF7.cov %>% group_by(dist) %>% summarise(COV.mean = mean(COV)) %>% mutate(COV.fold = log2(COV.mean/mean(COV.mean[c(1:20,180:200)]))) %>% mutate(cell = "MCF7")

# Open and format results

PARP1.MDA.MB231.cov <- read.table("./Figure_2/Repair/PARP1/GSM1517306.PARP1.MDA.MB231.hg19.IniSeq.cov.bedgraph", header = F)
colnames(PARP1.MDA.MB231.cov) <- c("chr", "start", "end", "INFO", "COV")
PARP1.MDA.MB231.cov <- PARP1.MDA.MB231.cov %>% separate(INFO, c("EFF", "bin"), sep = "_") %>% mutate(dist = (as.numeric(bin)-100)*100)

# Compute fold enrichment in coverage
PARP1.MDA.MB231.cov.plot <- PARP1.MDA.MB231.cov %>% group_by(dist) %>% summarise(COV.mean = mean(COV)) %>% mutate(COV.fold = log2(COV.mean/mean(COV.mean[c(1:20,180:200)]))) %>% mutate(cell = "MDA MB231") %>% mutate(COV.fold = replace(COV.fold, COV.fold > 0.70, NA))

# Bin by origins efficiency
PARP1.MDA.MB231.low.cov <- PARP1.MDA.MB231.cov %>% filter(EFF < 0.8418) %>% group_by(dist) %>% summarise(COV.mean = mean(COV)) %>% mutate(COV.fold = log2(COV.mean/mean(COV.mean[c(1:20,180:200)]))) %>% mutate(efficiency = "low")  
PARP1.MDA.MB231.medium.cov <- PARP1.MDA.MB231.cov %>% filter(EFF >= 0.8418 & EFF < 0.9065) %>% group_by(dist) %>% summarise(COV.mean = mean(COV)) %>% mutate(COV.fold = log2(COV.mean/mean(COV.mean[c(1:20,180:200)]))) %>% mutate(efficiency = "medium")  
PARP1.MDA.MB231.high.cov <- PARP1.MDA.MB231.cov %>% filter(EFF >= 0.9065) %>% group_by(dist) %>% summarise(COV.mean = mean(COV)) %>% mutate(COV.fold = log2(COV.mean/mean(COV.mean[c(1:20,180:200)]))) %>% mutate(efficiency = "high")  %>% mutate(COV.fold = replace(COV.fold, COV.fold > 1, NA))  
PARP1.MDA.MB231.eff.cov.plot <- rbind(PARP1.MDA.MB231.low.cov, PARP1.MDA.MB231.medium.cov, PARP1.MDA.MB231.high.cov)

# Plot coverage

PARP1.EFF.plot <- PARP1.MDA.MB231.eff.cov.plot  %>% ggplot(aes(x = dist, y = COV.fold, color = efficiency)) +
  geom_line(size = 0.5) +
  scale_color_manual(values=c("#FC4E07", "#0072B2", "#CC79A7")) + xlim(-5000,5000) +
  xlab("Distance from origin (bp)") + ylab("PARP1 coverage (log2 FC)") + ggtitle("PARP1 coverage in MDA MB231 around IniSeq origins") +
  theme_bw() + theme(aspect.ratio=1)

# Average PARP signal over cell lines
PARP.MCF10 <- PARP1.MCF10.cov.plot %>% dplyr::select(dist, COV.MCF10 = COV.fold )
PARP1.MDA.MB231 <- PARP1.MDA.MB231.cov.plot  %>% dplyr::select(dist, COV.MDA = COV.fold )
PARP <- PARP1.MCF7.cov.plot %>% dplyr::select(dist, COV.MCF7 = COV.fold ) %>% right_join(PARP.MCF10)  %>% right_join(PARP1.MDA.MB231) %>% mutate(COV.fold = (COV.MCF10+COV.MCF7+COV.MDA)/3, Factor = "PARP1") %>% replace(is.na(.), 0) %>% 
  dplyr::select(dist, COV.fold, Factor)

############################################################################
# Anaysis of RPA distribution

# Datasets #1
# GSM2033104 - HeLa - hg19 from Zhang et al. Mol Cell 2017  65 272-284

# Convert bigwig to bedgraph
# bigWigToBedGraph ./RPA/GSM2033104_h32.bw ./RPA/GSM2033104.RPA2.HeLa.hg19.bedgraph

# Compute mean coverage for each IniSeq windows using bedGraph inputs
# bedtools map -a IniSeq.10kb.hg19.100bp.windows.sorted.bed -b ./RPA/GSM2033104.RPA2.HeLa.hg19.bedgraph -c 4 -o mean -null 0 > ./RPA/GSM2033104.RPA2.HeLa.hg19.IniSeq.cov.bedgraph

# Open and format results

RPA2.HeLa.cov <- read.table("./Figure_2/Repair/RPA/GSM2033104.RPA2.HeLa.hg19.IniSeq.cov.bedgraph", header = F)
colnames(RPA2.HeLa.cov) <- c("chr", "start", "end", "INFO", "COV")
RPA2.HeLa.cov <- RPA2.HeLa.cov %>% separate(INFO, c("EFF", "bin"), sep = "_") %>% mutate(dist = (as.numeric(bin)-100)*100)

# Compute fold enrichment in coverage
RPA2.HeLa.cov.plot <- RPA2.HeLa.cov %>% group_by(dist) %>% summarise(COV.mean = mean(COV)) %>% mutate(COV.fold = log2(COV.mean/mean(COV.mean[c(1:20,180:200)]))) %>% mutate(cell = "RPA2")

# Datasets #2
# GSM2033106 - HeLa - hg19 from Zhang et al. Mol Cell 2017  65 272-284

# Convert bigwig to bedgraph
# bigWigToBedGraph ./RPA/GSM2033106_h34.bw ./RPA/GSM2033106.RPA1.HeLa.hg19.bedgraph

# Compute mean coverage for each IniSeq windows using bedGraph inputs
# bedtools map -a IniSeq.10kb.hg19.100bp.windows.sorted.bed -b ./RPA/GSM2033106.RPA1.HeLa.hg19.bedgraph -c 4 -o mean -null 0 > ./RPA/GSM2033106.RPA1.HeLa.hg19.IniSeq.cov.bedgraph

# Open and format results

RPA1.HeLa.cov <- read.table("./Figure_2/Repair/RPA/GSM2033106.RPA1.HeLa.hg19.IniSeq.cov.bedgraph", header = F)
colnames(RPA1.HeLa.cov) <- c("chr", "start", "end", "INFO", "COV")
RPA1.HeLa.cov <- RPA1.HeLa.cov %>% separate(INFO, c("EFF", "bin"), sep = "_") %>% mutate(dist = (as.numeric(bin)-100)*100)

# Compute fold enrichment in coverage
RPA1.HeLa.cov.plot <- RPA1.HeLa.cov %>% group_by(dist) %>% summarise(COV.mean = mean(COV)) %>% mutate(COV.fold = log2(COV.mean/mean(COV.mean[c(1:20,180:200)]))) %>% mutate(cell = "RPA1")

# Average RPA signals

RPA1 <- RPA1.HeLa.cov.plot %>% dplyr::select(dist, COV.RPA1 = COV.fold)
RPA <- RPA2.HeLa.cov.plot %>% dplyr::select(dist, COV.RPA2 = COV.fold ) %>% right_join(RPA1) %>%
  mutate(COV.fold = (COV.RPA1+COV.RPA2)/2, Factor = "RPA") %>% dplyr::select(dist, COV.fold, Factor)

############################################################################
# Anaysis of phosphorylated RPA distribution

# Datasets #1 - phospho-RPA2-S33
# GSM2891676 - HeLa - hg19 from Promonet et al Nat Commun 2020 11 3940

# Convert bigwig to bedgraph
# bigWigToBedGraph ./RPA/GSM2891676_prpa_HeLa_IP_2_1.bw ./RPA/GSM2891676.PRPA.HeLa.hg19.bedgraph

# Compute mean coverage for each IniSeq windows using bedGraph inputs
# bedtools map -a IniSeq.10kb.hg19.100bp.windows.sorted.bed -b ./RPA/GSM2891676.PRPA.HeLa.hg19.bedgraph -c 4 -o mean -null 0 > ./RPA/GSM2891676.PRPA.HeLa.hg19.IniSeq.cov.bedgraph

# Open and format results

PRPA.HeLa.cov <- read.table("./Figure_2/Repair/RPA/GSM2891676.PRPA.HeLa.hg19.IniSeq.cov.bedgraph", header = F)
colnames(PRPA.HeLa.cov) <- c("chr", "start", "end", "INFO", "COV")
PRPA.HeLa.cov <- PRPA.HeLa.cov %>% separate(INFO, c("EFF", "bin"), sep = "_") %>% mutate(dist = (as.numeric(bin)-100)*100)

# Compute fold enrichment in coverage
PRPA.HeLa.cov.plot <- PRPA.HeLa.cov %>% group_by(dist) %>% summarise(COV.mean = mean(COV)) %>% mutate(COV.fold = log2(COV.mean/mean(COV.mean[c(1:20,180:200)]))) %>% mutate(cell = "p-RPA2-S33")
PRPA <- PRPA.HeLa.cov.plot %>% mutate(Factor = "p-RPA") %>% dplyr::select(dist, COV.fold, Factor)

############################################################################
# Analysis of ORC proteins distribution

# Datasets #1 - ORC1
# GSM922790 - HeLa S3 - hg18 from Dellino et al. Genome Res. Biol. 2013 23 1-11 

# Compute mean coverage for each IniSeq windows using bedGraph inputs
# bedtools map -a IniSeq.10kb.hg18.100bp.windows.sorted.bed -b ./ORC1/GSM922790_ChIPseq_Orc1_GradientHela.bedGraph -c 4 -o mean -null 0 > ./ORC1/GSM922790.ORC1.HELA.S3.hg18.IniSeq.cov.bedgraph

# Open and format results

ORC1.HELA.cov <- read.table("./Figure_2/Repair/ORC1/GSM922790.ORC1.HELA.S3.hg18.IniSeq.cov.bedgraph", header = F)
colnames(ORC1.HELA.cov) <- c("chr", "start", "end", "INFO", "COV")
ORC1.HELA.cov <- ORC1.HELA.cov %>% separate(INFO, c("EFF", "bin"), sep = "_") %>% mutate(dist = (as.numeric(bin)-100)*100)

# Compute fold enrichment in coverage
ORC1.HELA.cov.plot <- ORC1.HELA.cov %>% group_by(dist) %>% summarise(COV.mean = mean(COV)) %>% mutate(COV.fold = log2(COV.mean/mean(COV.mean[c(1:20,180:200)]))) %>% mutate(cell = "ORC1 - HeLa")

# Datasets #2 - ORC2
# GSM1717888 - K562 - hg19 from Miotto et al. PNAS 2016 113 E4810-E4819 

# Convert bigwig to bedgraph
# bigWigToBedGraph ./ORC2/GSM1717888_Orc2_K562_rep1.bw ./ORC2/GSM1717888.ORC2.K562.hg19.bedgraph

# Compute mean coverage for each IniSeq windows using bedGraph inputs
# bedtools map -a IniSeq.10kb.hg19.100bp.windows.sorted.bed -b ./ORC2/GSM1717888.ORC2.K562.hg19.bedgraph -c 4 -o mean -null 0 > ./ORC2/GSM1717888.ORC2.K562.hg19.IniSeq.cov.bedgraph

# Open and format results

ORC2.K562.cov <- read.table("./Figure_2/Repair/ORC2/GSM1717888.ORC2.K562.hg19.IniSeq.cov.bedgraph", header = F)
colnames(ORC2.K562.cov) <- c("chr", "start", "end", "INFO", "COV")
ORC2.K562.cov <- ORC2.K562.cov %>% separate(INFO, c("EFF", "bin"), sep = "_") %>% mutate(dist = (as.numeric(bin)-100)*100)

# Compute fold enrichment in coverage
ORC2.K562.cov.plot <- ORC2.K562.cov %>% group_by(dist) %>% summarise(COV.mean = mean(COV)) %>% mutate(COV.fold = log2(COV.mean/mean(COV.mean[c(1:20,180:200)]))) %>% mutate(cell = "ORC1 - K562")

# Average ORC1 and ORC2 signals

ORC1 <- ORC1.HELA.cov.plot %>% dplyr::select(dist, COV.ORC1 = COV.fold)
ORC <- ORC2.K562.cov.plot %>% dplyr::select(dist, COV.ORC2 = COV.fold ) %>% right_join(ORC1) %>%
  mutate(COV.fold = (COV.ORC1+COV.ORC2)/2, Factor = "ORCs") %>% dplyr::select(dist, COV.fold, Factor)

############################################################################
# Anaysis of RAD51 distribution

# Datasets #1
# GSM2442936 - MCF10A - hg18 from Dellino et al Nat Genet 2019 51 1011-1023

# Convert bigwig to bedgraph
# bigWigToBedGraph ./RAD51/GSM2442936_MCF10AsiSIER_4OHT2h_RAD51_chip.bw ./RAD51/GSM2442936.RAD51.MCF10A.hg18.bedgraph

# Compute mean coverage for each IniSeq windows using bedGraph inputs
# bedtools map -a IniSeq.10kb.hg18.100bp.windows.sorted.bed -b ./RAD51/GSM2442936.RAD51.MCF10A.hg18.bedgraph -c 4 -o mean -null 0 > ./RAD51/GSM2442936.RAD51.MCF10A.hg18.IniSeq.cov.bedgraph

# Open and format results

RAD51.MCF10A.cov <- read.table("./Figure_2/Repair/RAD51/GSM2442936.RAD51.MCF10A.hg18.IniSeq.cov.bedgraph", header = F)
colnames(RAD51.MCF10A.cov) <- c("chr", "start", "end", "INFO", "COV")
RAD51.MCF10A.cov <- RAD51.MCF10A.cov %>% separate(INFO, c("EFF", "bin"), sep = "_") %>% mutate(dist = (as.numeric(bin)-100)*100)

# Compute fold enrichment in coverage
RAD51.MCF10A.cov.plot <- RAD51.MCF10A.cov %>% group_by(dist) %>% summarise(COV.mean = mean(COV)) %>% mutate(COV.fold = log2(COV.mean/mean(COV.mean[c(1:20,180:200)]))) %>% mutate(cell = "RAD51")

RAD51 <- RAD51.MCF10A.cov.plot %>% mutate(Factor = "RAD51") %>% dplyr::select(dist, COV.fold, Factor)

############################################################################
# Prepare and save plots

# Plot coverage of multiple DNA damage sensing / repair factors

Set.A <- rbind(TOP2B, BRCA1, XRCC5, BLISSseq, RAD51, ORC)
Set.A.plot <- Set.A  %>% ggplot(aes(x = dist, y = COV.fold, color = Factor)) +
  geom_line(size = 0.5) + xlim(-5000,5000) +
  scale_color_manual(values=c("#2EBAED", "#DE1C14", "#E69F00", "#999999", "#ADCC54", "#CC79A7")) +
  xlab("Distance from origin (bp)") + ylab("Enrichment over \nbackground (log2 FC)") + ggtitle("Set A") +
  theme_bw() + theme(aspect.ratio=1)

Set.B <- rbind(TOP1, PARP, PRPA, RPA)
Set.B.plot <- Set.B  %>% ggplot(aes(x = dist, y = COV.fold, color = Factor)) +
  geom_line(size = 0.5) + xlim(-5000,5000) +
  scale_color_manual(values=c("#2EBAED", "#DE1C14", "#E69F00", "#ADCC54", "#CC79A7")) +
  xlab("Distance from origin (bp)") + ylab("Enrichment over \nbackground (log2 FC)") + ggtitle("Set B") +
  theme_bw() + theme(aspect.ratio=1)

pdf("./Rplot/Fig_3B_3C.pdf", width=8, height=5, useDingbats=FALSE)
ggarrange(Set.A.plot, Set.B.plot, ncol = 2, nrow = 1)
dev.off()

# Origin efficiency dependency for selected factors

pdf("./Rplot/Fig_3E.pdf", width=15, height=4, useDingbats=FALSE)
ggarrange(TOP2B.EFF.plot, BRCA1.EFF.plot, TOP1.EFF.plot, PARP1.EFF.plot, ncol = 4, nrow = 1)
dev.off()








