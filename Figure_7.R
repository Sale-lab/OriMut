
# P. Murat, MRC Laboratory of Molecular Biology, July 2022
# Code for generating Figure 7 with associated supporting Figures
# in Murat et al. DNA replication initiation shapes the mutational landscape and expression of the human genome

library(dplyr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(stringr)
library(MutationalPatterns)
library(ggpubr)
library(reshape2)
library(lsa)
library(drc)
library(tidyr)
library(seqinr)
library(Biostrings)
library(vegan)
library(trackViewer)
library(zoo)
library(forcats)
library(Biostrings)

setwd("/Volumes/pmurat/OriMut")

############################################################################
############################################################################
###               Define probability density functions                   ###
############################################################################
############################################################################

############################################################################
# PDF for mutation rates at constitutive origins

population.SNP.IniSeq.closest <- read.csv("./Figure_1/population.SNP.IniSeq.closest.10kb.bed", sep="\t", header = F)
colnames(population.SNP.IniSeq.closest) <- c("chr", "pos", "REF", "ALT", "SAO", "CAF", "ori.pos", "EFF", "dist")
population.SNP.IniSeq.closest.bin <- population.SNP.IniSeq.closest %>% 
  mutate(dist = (as.numeric(cut(dist, breaks = seq(-10000, 10000, 100)))-100)*100) %>% 
  mutate(mut.type = ifelse((REF == "G" | REF == "C"), "GC", "AT")) %>% group_by(dist) %>% 
  summarise(GC.mut = sum(mut.type == "GC"), AT.mut = sum(mut.type == "AT")) %>%
  mutate(mut.rate = (GC.mut+AT.mut)/9351, mut.rate.fold = log2(mut.rate/mean(mut.rate[c(1:20,180:200)]))) %>% drop_na()
# Fit experimental data using using a multi-peak gaussian profile (nls function)

model.mutation.rate <- nls(mut.rate.fold ~ (C1*exp(-(dist-mean1)^2/(2*sigma1^2)) + C2*exp(-(dist-mean2)^2/(2*sigma2^2)) + C3*exp(-(dist-mean3)^2/(2*sigma3^2))),
                           data = population.SNP.IniSeq.closest.bin,
                           start=list(C1=0.15, mean1=0, sigma1=50, C2=-0.06, mean2=-1500, sigma2=50, C3=-0.06, mean3=1500, sigma3=50),
                           algorithm="port") 
# Build a symmetrical profile
model.mutation.rate.2 <- c(predict(model.mutation.rate)[1:100], rev(predict(model.mutation.rate)[1:100]))
# Prepare and save probability density function
MutFreq.pdf <- cbind.data.frame(dist = population.SNP.IniSeq.closest.bin$dist,
                        rate = model.mutation.rate.2) %>% mutate(prob = (2^rate)-0.5)
write.table(MutFreq.pdf, "./Figure_7/PDF_rate.tsv", sep="\t", col.names = T, row.names = F, quote = F)
MutFreq.pdf <- read.csv("./Figure_7/PDF_rate.tsv", sep = "\t")
# Save pdf plot
pdf("./Rplot/Fig_S7A.pdf", width=3, height=3, useDingbats=FALSE)
MutFreq.pdf %>% ggplot(aes(x=dist, y=prob)) +
  geom_line(size = 1) +
  labs(x = "Distance from origin (bp)", y = "Mutation probability") +
  theme_bw() + theme(aspect.ratio=0.9)
dev.off()

############################################################################
# PDF for COSMIC signatures exposure at constitutive origins

# In order to capture the mutational processes active at pre-origins
# the mutation matrix is corrected for base composition composition

# Load original 96 trinucleotide mutation count matrix

# TEST
mut.mat <- readRDS("./Figure_2/IniSeq_SNP_mutation_matrix.rds")

# Load base composition of IniSeq origins

IniSeq.inter.summary <- read.csv("./Figure_1/IniSeq.inter.10kb.base.composition.csv") %>% mutate(C.G = G+C, T.A = A+T) %>% dplyr::select(dist, C = C.G, T = T.A) %>% arrange(dist) %>% tibble::column_to_rownames("dist")
IniSeq.inter.summary <- t(IniSeq.inter.summary)

# Correct mutation matrix by the base composition at IniSeq origins

# Reorder columns
mut.mat <- mut.mat[ , order(as.numeric(as.character(colnames(mut.mat))))]
# Extract REF
mut.mat.cor <- as.data.frame(mut.mat) %>% add_rownames(var = "Tri.mut") %>% tidyr::separate(Tri.mut, c("D", "R", "A", "U"), sep = "\\[|>|]", remove = F) %>% dplyr::select(-D, -A, -U, -Tri.mut)
mut.mat.cor <- as.data.frame(mut.mat.cor)
row.names(mut.mat.cor) <- row.names(mut.mat)

# Create an empty matrix
mut.mat.cor.2 <- matrix(nrow = 0, ncol = 200)
# Fill it with corrected mutation counts
for (i in 1:nrow(mut.mat.cor)) {
  ref <- mut.mat.cor$R[i]
  base.comp.ref <- IniSeq.inter.summary[(which(row.names(IniSeq.inter.summary) == ref)),]
  mut.mat.cor.ref <- mut.mat.cor[i,] %>% dplyr::select(-R)
  p <- (mut.mat.cor.ref/base.comp.ref)*1000000   # mutations per Mb
  mut.mat.cor.2 <- rbind(mut.mat.cor.2, p)
}
# Reorder rows by original name
mut.mat.cor.3 <- mut.mat.cor.2[order(match(rownames(mut.mat.cor.2), row.names(mut.mat))), , drop = FALSE]

# Correct mutation matrix by the trinucleotide composition at IniSeq origins

# Extract  IniSeq origins trinucleotide composition
IniSeq.10kb.gr <- import("./Figure_7/IniSeq.ori.10kb.100nt.split.bed")
IniSeq.views <- Views(Hsapiens, IniSeq.10kb.gr)
IniSeq.tri <- trinucleotideFrequency(IniSeq.views, step=1, as.prob = FALSE)
IniSeq.tri.df <- cbind.data.frame(bin = IniSeq.10kb.gr$name, IniSeq.tri)
IniSeq.tri.summary <- IniSeq.tri.df %>% group_by(bin) %>% summarise_at(colnames(IniSeq.tri), sum) %>% mutate(dist = (as.numeric(bin)-100)*100)
row.names(IniSeq.tri.summary) <- IniSeq.tri.summary$dist

# Convert table to pyrimidine trinucleotide count
IniSeq.tri.summary.pyr <- IniSeq.tri.summary %>% dplyr::select(-bin, -dist)
IniSeq.tri.summary.pyr <- as.data.frame(t(IniSeq.tri.summary.pyr))
colnames(IniSeq.tri.summary.pyr) <- row.names(IniSeq.tri.summary)
IniSeq.tri.summary.pyr <- as.data.frame(IniSeq.tri.summary.pyr) %>% add_rownames(var = "tri")
# Need a loop to extract pyrimidine context
tri2 <- vector()
for (i in IniSeq.tri.summary.pyr$tri) {
  trinuc <- str_split(i, "")[[1]]
  ref <- as.character(trinuc[2])
  if (ref == "A" | ref == "G") {
    p <- toupper(c2s(rev(comp(trinuc))))} else {
      p <- paste(trinuc, collapse = "")
    }
  tri2 <- append(tri2, p)
}
IniSeq.tri.summary.pyr$tri2 <- tri2
IniSeq.tri.summary.pyr.2 <- IniSeq.tri.summary.pyr %>% dplyr::select(-tri, -tri2)
# Compute counts
IniSeq.tri.summary.pyr.3 <- IniSeq.tri.summary.pyr %>% dplyr::select(-tri) %>% group_by(tri2) %>% summarise_at(colnames(IniSeq.tri.summary.pyr.2), sum)
IniSeq.tri.summary.pyr.3 <- as.data.frame(IniSeq.tri.summary.pyr.3)

# Correct mutation matrix by the observed frequencies of trinucleotides at IniSeq origins
mut.mat.cor.5 <- as.data.frame(mut.mat.cor.3) %>% add_rownames(var = "Tri.mut") %>% tidyr::separate(Tri.mut, c("D", "R", "A", "U"), sep = "\\[|>|]", remove = F) %>% 
  mutate(tri = paste0(D, R, U, sep = "")) %>% dplyr::select(-D, -R, -A, -U, -Tri.mut)
mut.mat.cor.5 <- as.data.frame(mut.mat.cor.5)
row.names(mut.mat.cor.5) <- row.names(mut.mat.cor.3)
# Create an empty matrix
mut.mat.cor.6 <- matrix(nrow = 0, ncol = 200)
for (i in 1:nrow(mut.mat.cor.5)) {
  trinuc <- mut.mat.cor.5$tri[i]
  trinuc.count <- IniSeq.tri.summary.pyr.3 %>% filter(tri2 == trinuc) %>% dplyr::select(-tri2)
  trinuc.count <- trinuc.count[,order(names(trinuc.count))]
  count.mut <- mut.mat.cor.5 %>% filter(tri == trinuc) %>% dplyr::select(-tri)
  count.mut <- count.mut[,order(names(count.mut))]
  p <- mapply('/', count.mut, trinuc.count)*1000000  # mutation per Mb
  rownames(p) <- rownames(count.mut)
  mut.mat.cor.6 <- rbind(mut.mat.cor.6, p)
}
# Remove duplicated rows
mut.mat.cor.7 <-mut.mat.cor.6[!duplicated(mut.mat.cor.6), ]
# Reorder rows by original name
mut.mat.cor.7 <- mut.mat.cor.7[order(match(rownames(mut.mat.cor.7), row.names(mut.mat.cor.3))), drop = FALSE]

# Fit COSMIC signatures to the 96 trinucleotide mutation count matrix corrected for base composition

# Find optimal contribution of COSMIC signatures to mutational landscape of IniSeq origins
# Use strict refitting 

COSMIC.signatures <- get_known_signatures()
COSMIC.fit.IniSeq.cor <- fit_to_signatures_strict(mut.mat.cor.7, COSMIC.signatures, max_delta = 0.004)

# Assess contribution of relevant COSMIC signatures to the mutational landscape of IniSeq origins corrected for base composition
COSMIC.fit.IniSeq.cor.2 <- as.data.frame(COSMIC.fit.IniSeq.cor$fit_res$contribution)
COSMIC.fit.IniSeq.cor.2 <- apply(COSMIC.fit.IniSeq.cor.2,2,function(x){x/sum(x)})
COSMIC.fit.IniSeq.cor.2 <- as.data.frame(t(COSMIC.fit.IniSeq.cor.2)) %>% select_if(colSums(.) != 0) %>% add_rownames(var = "dist")
COSMIC.fit.IniSeq.cor.2$dist <- as.numeric(as.character(COSMIC.fit.IniSeq.cor.2$dist))

# Remove signatures with significant exposure in less than five bins and renormalise dataframe

COSMIC.fit.IniSeq.cor.3 <- COSMIC.fit.IniSeq.cor.2[, which(as.numeric(colSums(COSMIC.fit.IniSeq.cor.2 != 0)) > 5)]
COSMIC.fit.IniSeq.cor.4 <- COSMIC.fit.IniSeq.cor.3[,2:8]/rowSums(COSMIC.fit.IniSeq.cor.3[,2:8])
COSMIC.fit.IniSeq.cor.4 <- cbind.data.frame(dist = COSMIC.fit.IniSeq.cor.3$dist, COSMIC.fit.IniSeq.cor.4)

# Signature identified
# SBS1, SBS5, SBS10b, SBS12, SBS26, SBS37, SBS87

# Fit exposures to prepare PDFs 

Sig.pdf.plot <- melt(COSMIC.fit.IniSeq.cor.4, id.vars = "dist")
colnames(Sig.pdf.plot) <- c("dist", "Signature", "Probability")

# SBS1
Sig.pdf.SBS1 <- Sig.pdf.plot %>% filter(Signature == "SBS1")
Sig.pdf.SBS1 <- Sig.pdf.SBS1[order(Sig.pdf.SBS1$dist),]
# Fit using a sigmoidal on half of the profile
Sig.pdf.SBS1.pos <- Sig.pdf.SBS1 %>% filter(dist > 0)
model.Sig.pdf.SBS1.pos <- drm(Probability ~ log(dist), fct = G.4(), data = Sig.pdf.SBS1.pos)  # G.4 : Gompertz curve
# Reconstruct profile
model.Sig.pdf.SBS1 <- c(rev(predict(model.Sig.pdf.SBS1.pos)), predict(model.Sig.pdf.SBS1.pos))
# Check fit
plot(Sig.pdf.SBS1$dist, Sig.pdf.SBS1$Probability)
points(Sig.pdf.SBS1$dist, model.Sig.pdf.SBS1, col = "red")

# SBS5
Sig.pdf.SBS5 <- Sig.pdf.plot %>% filter(Signature == "SBS5")
Sig.pdf.SBS5 <- Sig.pdf.SBS5[order(Sig.pdf.SBS5$dist),]
# Fit using a multi gaussian profile
model.Sig.pdf.SBS5 <- nls(Probability-0.44 ~ (C1*exp(-(dist-mean1)^2/(2*sigma1^2)) + C2*exp(-(dist-mean2)^2/(2*sigma2^2)) + C3*exp(-(dist-mean3)^2/(2*sigma3^2))),
                          data = Sig.pdf.SBS5,
                          start=list(C1=0.2, mean1=-1000, sigma1=50, C2=-0.4, mean2=0, sigma2=50, C3=0.2, mean3=1000, sigma3=50),
                          algorithm="port")
# Check fit
plot(Sig.pdf.SBS5$dist, Sig.pdf.SBS5$Probability-0.44)
points(Sig.pdf.SBS5$dist, predict(model.Sig.pdf.SBS5), col = "red")
# Reconstruct profile
model.Sig.pdf.SBS5.2 <- predict(model.Sig.pdf.SBS5)+0.44
# Make the model symetrical around 0
model.Sig.pdf.SBS5.3 <- c(rev(model.Sig.pdf.SBS5.2[101:200]), model.Sig.pdf.SBS5.2[101:200])
# Assign 0 to distance values bin 100 and 101 to fit observed values
model.Sig.pdf.SBS5.3[100] <- 0
model.Sig.pdf.SBS5.3[101] <- 0
# Check final fit
plot(Sig.pdf.SBS5$dist, Sig.pdf.SBS5$Probability)
points(Sig.pdf.SBS5$dist, model.Sig.pdf.SBS5.3, col = "red")

# SBS10b
Sig.pdf.SBS10b <- Sig.pdf.plot %>% filter(Signature == "SBS10b")
Sig.pdf.SBS10b <- Sig.pdf.SBS10b[order(Sig.pdf.SBS10b$dist),]
# Fit using a sigmoidal on half of the profile using a log scale
Sig.pdf.SBS10b.neg <- Sig.pdf.SBS10b %>% filter(dist < 0)
model.Sig.pdf.SBS10b.neg <- drm(Probability ~ log(-dist), fct = G.4(), data = Sig.pdf.SBS10b.neg)   # G.4 : Gompertz curve
# Check fit
plot(log(-Sig.pdf.SBS10b.neg$dist), Sig.pdf.SBS10b.neg$Probability)
points(log(-Sig.pdf.SBS10b.neg$dist),predict(model.Sig.pdf.SBS10b.neg), col = "red")
# Reconstruct profile
model.Sig.pdf.SBS10b <- c(predict(model.Sig.pdf.SBS10b.neg), 0, 0, rev(predict(model.Sig.pdf.SBS10b.neg)))
# Check fit
plot(Sig.pdf.SBS10b$dist, Sig.pdf.SBS10b$Probability)
points(Sig.pdf.SBS10b$dist, model.Sig.pdf.SBS10b, col = "red")

# SBS12
Sig.pdf.SBS12 <- Sig.pdf.plot %>% filter(Signature == "SBS12")
Sig.pdf.SBS12 <- Sig.pdf.SBS12[order(Sig.pdf.SBS12$dist),]
# Fit a gaussian
model.Sig.pdf.SBS12 <- nls(Probability ~ C1*exp(-(dist-mean1)^2/(2*sigma1^2)),
                           data = Sig.pdf.SBS12,
                           start=list(C1=0.1, mean1=0, sigma1=50),
                           algorithm="port") 
# Make the model symetrical around 0
model.Sig.pdf.SBS12.3 <- c(predict(model.Sig.pdf.SBS12)[0:100], rev(predict(model.Sig.pdf.SBS12)[0:100]))
# Check final fit
plot(Sig.pdf.SBS12$dist, Sig.pdf.SBS12$Probability)
points(Sig.pdf.SBS12$dist, model.Sig.pdf.SBS12.3, col = "red")

# SBS26
Sig.pdf.SBS26 <- Sig.pdf.plot %>% filter(Signature == "SBS26")
Sig.pdf.SBS26 <- Sig.pdf.SBS26[order(Sig.pdf.SBS26$dist),]
# Fit using a sigmoidal on half of the profile using a log scale
Sig.pdf.SBS26.pos <- Sig.pdf.SBS26 %>% filter(dist > 0)
model.Sig.pdf.SBS26.pos <- drm(Probability ~ log(dist), fct = G.4(), data = Sig.pdf.SBS26.pos)   # G.4 : Gompertz curve
# Reconstruct profile
model.Sig.pdf.SBS26 <- c(rev(predict(model.Sig.pdf.SBS26.pos)), predict(model.Sig.pdf.SBS26.pos))
# Check fit
plot(Sig.pdf.SBS26$dist, Sig.pdf.SBS26$Probability)
points(Sig.pdf.SBS26$dist, model.Sig.pdf.SBS26, col = "red")

# SBS37
Sig.pdf.SBS37 <- Sig.pdf.plot %>% filter(Signature == "SBS37")
Sig.pdf.SBS37 <- Sig.pdf.SBS37[order(Sig.pdf.SBS37$dist),]
# Fit using a sigmoidal on half of the profile using a log scale
Sig.pdf.SBS37.pos <- Sig.pdf.SBS37 %>% filter(dist > 0)
model.Sig.pdf.SBS37.pos <- drm(Probability ~ log(dist), fct = G.4(), data = Sig.pdf.SBS37.pos)   # G.4 : Gompertz curve
# Reconstruct profile
model.Sig.pdf.SBS37 <- c(rev(predict(model.Sig.pdf.SBS37.pos)), predict(model.Sig.pdf.SBS37.pos))
# Check fit
plot(Sig.pdf.SBS37$dist, Sig.pdf.SBS37$Probability)
points(Sig.pdf.SBS37$dist, model.Sig.pdf.SBS37, col = "red")

# SBS87
Sig.pdf.SBS87 <- Sig.pdf.plot %>% filter(Signature == "SBS87")
Sig.pdf.SBS87 <- Sig.pdf.SBS87[order(Sig.pdf.SBS87$dist),]
# Fit using a sigmoidal on half of the profile using a log scale
Sig.pdf.SBS87.neg <- Sig.pdf.SBS87 %>% filter(dist < 0)
model.Sig.pdf.SBS87.neg <- drm(Probability ~ log(-dist), fct = G.4(), data = Sig.pdf.SBS87.neg)   # G.4 : Gompertz curve
# Check fit
plot(log(-Sig.pdf.SBS87.neg$dist), Sig.pdf.SBS87.neg$Probability)
points(log(-Sig.pdf.SBS87.neg$dist),predict(model.Sig.pdf.SBS87.neg), col = "red")
# Reconstruct profile
model.Sig.pdf.SBS87 <- c(predict(model.Sig.pdf.SBS87.neg), 0.11, 0.11, rev(predict(model.Sig.pdf.SBS87.neg)))
# Check fit
plot(Sig.pdf.SBS87$dist, Sig.pdf.SBS87$Probability)
points(Sig.pdf.SBS87$dist, model.Sig.pdf.SBS87, col = "red")

# Reconstruct PDF with fitted profiles

Sig.pdf.fit <- cbind.data.frame(dist = Sig.pdf.SBS1$dist,
                                SBS1 = model.Sig.pdf.SBS1,
                                SBS5 = model.Sig.pdf.SBS5.3,
                                SBS10b = model.Sig.pdf.SBS10b,
                                SBS12 = model.Sig.pdf.SBS12.3,
                                SBS26 = model.Sig.pdf.SBS26,
                                SBS37 = model.Sig.pdf.SBS37.3,
                                SBS87 = model.Sig.pdf.SBS87)

# Normalise probability to 1
Sig.pdf.fit.norm <- Sig.pdf.fit[,-1]/rowSums(Sig.pdf.fit[,-1])
Sig.pdf.fit.norm$dist <- Sig.pdf.fit$dist

# Prepare plot
Sig.pdf.fit.norm.plot <- melt(Sig.pdf.fit.norm, "dist", variable.name = "Signature", value.name = "Probability")

pdf("./Rplot/Fig_S7B.pdf", width=5, height=5, useDingbats=FALSE)
Sig.pdf.fit.norm.plot %>% ggplot(aes(x=dist, y=Probability, color = Signature)) +
  geom_line(size = 1) +
  scale_color_manual(values = c("#999999", "#D4D2D2", "#ADCC54","#FFEDA0", "#FD8D3C", "#08306B", "#D6604D", "#6BAED6")) +
  labs(x = "Distance from origin (bp)", y = "Mutational signature probability") +
  theme_bw() + theme(aspect.ratio=0.9)
dev.off()

# Assign 0 to negative probability
Sig.pdf.fit.norm[Sig.pdf.fit.norm$SBS1 < 0,]$SBS1 <- 0
Sig.pdf.fit.norm[Sig.pdf.fit.norm$SBS26 < 0,]$SBS26 <- 0
Sig.pdf.fit.norm[Sig.pdf.fit.norm$SBS87 < 0,]$SBS87 <- 0

# Save PDFs
write.table(Sig.pdf.fit.norm, "./Figure_7/PDF_mutation_fit.tsv", sep="\t", col.names = T, row.names = F, quote = F)

############################################################################
# PDF for replication strand bias

# Open result of prior analysis

CG.AT.IniSeq.assymetry <- read.csv("./Figure_1/CG.AT.IniSeq.assymetry.csv")
CG.GC.IniSeq.assymetry <- read.csv("./Figure_1/CG.GC.IniSeq.assymetry.csv")
TA.CG.IniSeq.assymetry <- read.csv("./Figure_1/TA.CG.IniSeq.assymetry.csv")
TA.GC.IniSeq.assymetry <- read.csv("./Figure_1/TA.GC.IniSeq.assymetry.csv")

# C>A and T>C replication strand bias profiles display sigmoidal shapes
# Fitting using drm function

model.IniSeq.CG.AT <- drm(bias ~ Distance, fct = L.4(), data = CG.AT.IniSeq.assymetry)
model.IniSeq.TA.CG <- drm(bias ~ Distance, fct = L.4(), data = TA.CG.IniSeq.assymetry)

# C>G and T>A replication strand bias profiles display sigmoidal shapes
# Fitting using drm function on a loess fit to minimise noise

lo <- loess(bias~Distance, data = CG.GC.IniSeq.assymetry, span = 0.1)
lo.fit <- predict(lo)
CG.GC.IniSeq.assymetry <- CG.GC.IniSeq.assymetry %>% mutate(lo = lo.fit)
model.IniSeq.CG.GC <- drm(lo.fit ~ Distance, fct = L.4(), data = CG.GC.IniSeq.assymetry)

lo <- loess(bias~Distance, data = TA.GC.IniSeq.assymetry, span = 0.1)
lo.fit <- predict(lo)
TA.GC.IniSeq.assymetry <- TA.GC.IniSeq.assymetry %>% mutate(lo = lo.fit)
model.IniSeq.TA.GC <- drm(lo.fit ~ Distance, fct = L.4(), data = TA.GC.IniSeq.assymetry)

# C>T replication strand bias profile displays a multi-peak gaussian profile
# Fitting using nls function

model.IniSeq.CG.TA <- nls(bias ~ (C1*exp(-(Distance-mean1)^2/(2*sigma1^2)) + C2*exp(-(Distance-mean2)^2/(2*sigma2^2))),
                          data = CG.TA.IniSeq.assymetry,
                          start=list(C1=-0.3, mean1=-100, sigma1=50, C2=0.3, mean2=100, sigma2=50),
                          algorithm="port") 
# Make fit symmetrical
model.IniSeq.CG.TA.2 <- c(-rev(predict(model.IniSeq.CG.TA)[101:200]), predict(model.IniSeq.CG.TA)[101:200])

# T>A replication strand bias profile is a flat profile

model.IniSeq.TA.AT <- rep(0, 200)

# Combine profiles

Assymetry.profiles <- cbind.data.frame(dist = CG.AT.IniSeq.assymetry$Distance,
                        'C>A' = predict(model.IniSeq.CG.AT),
                        'G>T' = rev(predict(model.IniSeq.CG.AT)),
                        'C>G' = predict(model.IniSeq.CG.GC),
                        'G>C' = rev(predict(model.IniSeq.CG.GC)),
                        'C>T' = model.IniSeq.CG.TA.2,
                        'G>A' = rev(model.IniSeq.CG.TA.2),
                        'T>A' = model.IniSeq.TA.AT,
                        'A>T' = rev(model.IniSeq.TA.AT),
                        'T>C' = predict(model.IniSeq.TA.CG),
                        'A>G' = rev(predict(model.IniSeq.TA.CG)),
                        'T>G' = predict(model.IniSeq.TA.GC),
                        'A>C' = rev(predict(model.IniSeq.TA.GC)))

# Prepare PDFs
# Convert log2 ratio to observed frequency and then probability
# the probability of fixing a given mutation is given as P = fct/(1+fct) where fct is the observed frequency
Assymetry.pdf <- cbind.data.frame(dist = Assymetry.profiles$dist, (2^(Assymetry.profiles[,-1]))/(1+(2^(Assymetry.profiles[,-1]))))

Assymetry.pdf.C.mut.plot <- melt(Assymetry.pdf, "dist", variable.name = "Mutation", value.name = "Probability") %>% filter(Mutation == "C>A" | Mutation == "C>G" | Mutation == "C>T")
Assymetry.pdf.T.mut.plot <- melt(Assymetry.pdf, "dist", variable.name = "Mutation", value.name = "Probability") %>% filter(Mutation == "T>A" | Mutation == "T>C" | Mutation == "T>G")
Assymetry.pdf.G.mut.plot <- melt(Assymetry.pdf, "dist", variable.name = "Mutation", value.name = "Probability") %>% filter(Mutation == "G>A" | Mutation == "G>C" | Mutation == "G>T")
Assymetry.pdf.A.mut.plot <- melt(Assymetry.pdf, "dist", variable.name = "Mutation", value.name = "Probability") %>% filter(Mutation == "A>C" | Mutation == "A>G" | Mutation == "A>T")

a <- Assymetry.pdf.C.mut.plot %>% ggplot(aes(x=dist, y=Probability, color = Mutation)) +
  geom_line(size = 1) +
  scale_color_manual(values = c("#2EBAED", "#999999", "#DE1C14")) +
  labs(x = "Distance from origin (bp)", y = "Replication bias probability", title = "C mutations") +
  theme_bw() + theme(aspect.ratio=0.9)

b <- Assymetry.pdf.T.mut.plot %>% ggplot(aes(x=dist, y=Probability, color = Mutation)) +
  geom_line(size = 1) +
  scale_color_manual(values = c("#0072B2", "#ADCC54", "#D55E00")) +
  labs(x = "Distance from origin (bp)", y = "Replication bias probability", title = "T mutations") +
  theme_bw() + theme(aspect.ratio=0.9)

pdf("./Rplot/Fig_S7C.pdf", width=8, height=4, useDingbats=FALSE)
ggarrange(a, b, ncol = 2, nrow = 1)
dev.off()

# Save PDF
write.table(Assymetry.pdf, "./Figure_7/PDF_asymmetry.tsv", sep="\t", col.names = T, row.names = F, quote = F)

############################################################################
############################################################################
###                   In Silico model of evolution                       ###
############################################################################
############################################################################

# The model is based on 5 independent steps/functions:                                                                                      
# (i) Select the number of mutations per 100 bp bins according to the observed mutation rate around origins                                 
# (ii) Select the nature of the mutation events (in their trinucleotide context) according to the observed exposures to COSMIC signatures   
# (iii) Select the strand in which the mutations are fixed according to the observed replication strand bias                                
# (iv) Replace mutated bases in original sequence                                                                                           
# (v) Perform several rounds of in silico evolution to recapitulate origin sequences

# Open PDFs

MutFreq.pdf <- read.csv("./Figure_7/PDF_rate.tsv", sep = "\t")
cosmix.prob <- read.csv("./Figure_7/PDF_mutation_fit.tsv", sep = "\t")
Assymetry.pdf <- read.csv("./Figure_7/PDF_asymmetry.tsv", sep = "\t")
colnames(Assymetry.pdf) <- c("dist","C>A","G>T","C>G","G>C","C>T","G>A","T>A","A>T","T>C","A>G","T>G","A>C")

############################################################################
# (i) Select the number of mutations per 100 bp bins according to the observed mutation rate around origins
# Take a number of total number of bases to mutate (MutNumber) and output their distribution around the origin

BinSelect <- function(MutNumber) {
  sample(MutFreq.pdf$dist, MutNumber, replace = T, prob = MutFreq.pdf$prob)
}

############################################################################
# (ii) Select the nature of the mutation events (in their trinucleotide context) according to the observed exposures to COSMIC signatures

# Recover COSMIC signatures

# Recover mutational signatures
cosmix <- get_known_signatures()
Ori.sig <- c("SBS1", "SBS5", "SBS6", "SBS10b", "SBS12", "SBS26", "SBS37", "SBS87")
cosmix.ori <- as.data.frame(cosmix[, which(colnames(cosmix) %in% Ori.sig)])

# Recover mutation types
nmf_res <- readRDS("./Figure_2/Signatures_IniSeq_SNP.rds")
cosmix.ori.seq <- rownames(as.data.frame(nmf_res$signatures))

# prepare mutation table

Mut <- vector()
Mut.context <- vector()
tri.ref.watson <- vector()
tri.alt.watson <- vector()
tri.ref.crick <- vector()
tri.alt.crick <- vector()
for (i in 1:length(cosmix.ori.seq)) {
  p <- cosmix.ori.seq[i]
  D <- strsplit(p, split = "")[[1]][1]   # down
  R <- strsplit(p, split = "")[[1]][3]   # REF
  U <- strsplit(p, split = "")[[1]][7]   # up
  A <- strsplit(p, split = "")[[1]][5]   # ALT
  q <- paste0(D, R, U)
  r <- paste0(D, A, U)
  s <- toupper(c2s(rev(comp(c(D, R, U)))))
  t <- toupper(c2s(rev(comp(c(D, A, U)))))
  u <- paste0(R, ">", A)
  Mut.context <- append(Mut.context, p)
  tri.ref.watson <- append(tri.ref.watson, q)
  tri.alt.watson <- append(tri.alt.watson, r)
  tri.ref.crick <- append(tri.ref.crick, s)
  tri.alt.crick <- append(tri.alt.crick, t)
  Mut <- append(Mut, u)
}
cosmix.df.freq <- cbind.data.frame(Mut, Mut.context, tri.ref.watson, tri.alt.watson, tri.ref.crick, tri.alt.crick, cosmix.ori)
# replace 0 values by 0.000001
cosmix.df.freq[cosmix.df.freq == 0] <- 0.000001

Mut.table <- cosmix.df.freq %>% dplyr::select(Mut, Mut.context, tri.ref.watson, tri.alt.watson, tri.ref.crick, tri.alt.crick)

# Define function that accept a list of bins to mutate as input (MutBins) and ouput a mutation event as output

MutEvent <- function(x) {   # x is a row from MutBins generated by BaseSelect
  bin <- as.numeric(x[2])
  prob.sig <- cosmix.prob %>% filter(dist == bin) %>% dplyr::select(-dist)
  # Select a signature based on the probability
  selected.sig <- sample(colnames(prob.sig), size = 1, prob = prob.sig[1,])
  # Select mutational signature probabilities
  prob.mut <- cosmix.df.freq %>% dplyr::select(selected.sig)
  # Select mutation events based on the selected COSMIC signature
  mut <- sample(cosmix.df.freq$Mut.context, size = 1, prob = prob.mut[,1])
  return(mut)
}

############################################################################
# (iii) Select the strand in which the mutations are fixed according to the observed replication strand bias

# Define function that accept 2 parameters as input (bin, Mut) and output a mutational event as an output
# based on the observed replication strand bias

MutSelect <- function(x) { # x is a row from the MutEvent dataframe generated by binding the results from MutEvent to MutBins
  bin <- as.numeric(x[1])
  Mut <- as.character(x[3])
  tri.ref.watson <- x[4]
  tri.alt.watson <- x[5]
  tri.ref.crick <- x[6]
  tri.alt.crick <- x[7]
  # Select bin in asymmetry PDF
  Asy.bin <- Assymetry.pdf %>% filter(dist == bin) %>% dplyr::select(Mut)
  tri.ref.select <- sample(c(tri.ref.watson, tri.ref.crick), size = 1, prob = c(Asy.bin[,1], 1-Asy.bin[,1]))
  return(tri.ref.select)
}

############################################################################
# (iv) Replace mutated bases in original sequence 

# Write function that input a sequence, a bin number, a ref trinucleotide and a alt trinucleotide to output the modified sequence

MutFix <- function(seq, x){ # x is a row from the MutSelect dataframe generated by binding the results from MutSelect to MutEvent
  bin <- as.numeric(x[1])+9900
  ref <- as.character(x[2])
  alt <-as.character(x[3])
  seq.ref <- toupper(c2s(seq))
  # identify all ref patterns and randomly select one to mutate
  ref.pos <- str_locate_all(seq.ref, ref)[[1]][,1]
  ref.pos.bin <- ref.pos[which(ref.pos >= bin & ref.pos < bin + 100)]
  if (length(ref.pos.bin) == 0) {
    seq.mut <- seq} else {
      ref.pos.select <- sample(ref.pos.bin, size = 1)
      # mutate
      sub.seq.mut <- seq.ref %>% substr(ref.pos.select-1, ref.pos.select+3) %>% str_replace(ref, alt)
      seq.ref.mut <- paste0(substr(seq.ref, 1, ref.pos.select-2), sub.seq.mut, substr(seq.ref, ref.pos.select+4,nchar(seq.ref)))
      # convert sequenc to vector
      seq.mut <- s2c(tolower(seq.ref.mut))
    }
  return(seq.mut)
}

# Write function that input a sequence (seq), the number of bases to mutate (n) and output a mutated seq (seq.mut)

MutSeq <- function(seq, MutNumber) {
  # select the bins to mutate
  MutBins.df <- cbind.data.frame(Exp = c(1:MutNumber), bin = BinSelect(MutNumber))
  # compute a mutation events
  MutEvent.df <- cbind.data.frame(MutBins.df, Mut.context = apply(MutBins.df, 1, MutEvent)) %>% left_join(Mut.table, by = "Mut.context") %>% dplyr::select(-Exp)
  # select bases to be mutated and resulting mutations (in their trinucleotide context)
  MutSelect.df <- cbind.data.frame(MutEvent.df, tri.ref.select = apply(MutEvent.df, 1, MutSelect)) %>% mutate(tri.alt.select = case_when(tri.ref.select == tri.ref.watson ~ tri.alt.watson, TRUE ~ tri.alt.crick)) %>% 
    dplyr::select(bin, tri.ref.select, tri.alt.select)
  n <- nrow(MutSelect.df)
  out.seq <- vector("list", n)
  out.seq[[1]] <- MutFix(seq, MutSelect.df[1,])
  for (i in 2:n) out.seq[[i]] <- MutFix(out.seq[[i-1]], MutSelect.df[i-1,])
  seq.mut <- out.seq[[n]]
  return(seq.mut)
}

############################################################################
# (v) Perform several rounds of in silico evolution to recapitulate origin sequences

# Generate a multi fasta file of random 20kb DNA sequences

# Random sequences are generated using the RSAT - random sequence generator
# http://rsat.sb-roscoff.fr/random-seq_form.cgi
# using a Homo_sapiens_GRCh38 specific Markov model with an oligonucleotide size of 6
# from DNA sequences calibrated on non-coding upstream sequences

# Options
# Sequence length 20000
# Number 100
# Format fasta
# Line width 0
# Mask none
# Organism Homo_sapiens_GRCh38
# DNA sequences calibrated on non-coding upstream sequences
# oligo size 6

###########################################################################################
# Apply one round of evolution to a multi fasta file

# Write a function to evolve a multifasta file

Evolve <- function(x, n) { # input is a DNAStringSet object
  line_header <- vector()
  DNA_sequence <- vector()
  for(i in 1:length(x)) {
    p <- paste("Seq_", i, sep="")
    q <- c2s(toupper(MutSeq(s2c(c2s(x[[i]])), n)))
    line_header <- append(line_header, p)
    DNA_sequence <- append(DNA_sequence, q)
  }
  x.mut <- DNAStringSet(DNA_sequence)
  names(x.mut) <- line_header
  return(x.mut) # output is a DNAStringSet object
}

############################################################################
# (vi) In Silico Evolution

# Apply rounds of evolution 
# Starting from Seqs_c_0.fasta

n <- 100  # bases to mutate
for (i in 0:5000) {
  print(i)
  print(Sys.time())
  infile <- paste("./Figure_7/Exp1/", "Seqs_c_", i, ".fasta", sep = "")
  outfile <- paste("./Figure_7/Exp1/", "Seqs_c_", (i+1), ".fasta", sep = "")
  seqs <- readDNAStringSet(infile)
  seqs.mut <- Evolve(seqs, n)
  writeXStringSet(seqs.mut, outfile)
}

# Experiments were performed in batches of 100 sequences
# 10 experiments for 1,000 sequences

############################################################################
# Combine batch experiments

# Combine sequences for cycles of interest

set <- c(0,2,4,8,16,32,64,128,256,512,1024,2048,4096)

for (i in (set)) {
  print(i)
  Exp1.file.name <- paste("./Figure_7/Exp1/Seqs_c_", i, ".fasta", sep ="") 
  Exp2.file.name <- paste("./Figure_7/Exp2/Seqs_c_", i, ".fasta", sep ="") 
  Exp3.file.name <- paste("./Figure_7/Exp3/Seqs_c_", i, ".fasta", sep ="") 
  Exp4.file.name <- paste("./Figure_7/Exp4/Seqs_c_", i, ".fasta", sep ="") 
  Exp5.file.name <- paste("./Figure_7/Exp5/Seqs_c_", i, ".fasta", sep ="") 
  Exp6.file.name <- paste("./Figure_7/Exp6/Seqs_c_", i, ".fasta", sep ="") 
  Exp7.file.name <- paste("./Figure_7/Exp7/Seqs_c_", i, ".fasta", sep ="") 
  Exp8.file.name <- paste("./Figure_7/Exp8/Seqs_c_", i, ".fasta", sep ="") 
  Exp9.file.name <- paste("./Figure_7/Exp9/Seqs_c_", i, ".fasta", sep ="") 
  Exp10.file.name <- paste("./Figure_7/Exp10/Seqs_c_", i, ".fasta", sep ="")
  Exp1.seqs <- readDNAStringSet(Exp1.file.name)
  Exp2.seqs <- readDNAStringSet(Exp2.file.name)
  Exp3.seqs <- readDNAStringSet(Exp3.file.name)
  Exp4.seqs <- readDNAStringSet(Exp4.file.name)
  Exp5.seqs <- readDNAStringSet(Exp5.file.name)
  Exp6.seqs <- readDNAStringSet(Exp6.file.name)
  Exp7.seqs <- readDNAStringSet(Exp7.file.name)
  Exp8.seqs <- readDNAStringSet(Exp8.file.name)
  Exp9.seqs <- readDNAStringSet(Exp9.file.name)
  Exp10.seqs <- readDNAStringSet(Exp10.file.name)
  all.seqs <- DNAStringSet(c(Exp1.seqs, Exp2.seqs, Exp3.seqs, Exp4.seqs, Exp5.seqs, Exp6.seqs, Exp7.seqs, Exp8.seqs, Exp9.seqs, Exp10.seqs))
  seq_id <- c(1:1000)
  seq_names <- paste("Seq_", seq_id, sep ="")
  names(all.seqs) <- seq_names
  outfile <- paste("./Figure_7/Exp1000/", "Seqs_c_", i, ".fasta", sep = "")
  writeXStringSet(all.seqs, outfile)
}

############################################################################
############################################################################
###           Analysis of in silico evolution experiment                 ###
###                     Sequence features                                ###
############################################################################
############################################################################

############################################################################
# Create functions for computing base composition and other sequence features in rolling windows

# Define regex and functions

# G-quadruplex regex

gquad <- "G{3,}.{1,7}G{3,}.{1,7}G{3,}.{1,7}G{3,}"
cquad <- "C{3,}.{1,7}C{3,}.{1,7}C{3,}.{1,7}C{3,}"

# Count average base composition in bins

SeqFeat <- function(x) { # input is a multi fasta file
  A.mat <- matrix(nrow = 0, ncol = 200)
  C.mat <- matrix(nrow = 0, ncol = 200)
  G.mat <- matrix(nrow = 0, ncol = 200)
  T.mat <- matrix(nrow = 0, ncol = 200)
  GA.mat <- matrix(nrow = 0, ncol = 200)
  GC.mat <- matrix(nrow = 0, ncol = 200)
  CpG.mat <- matrix(nrow = 0, ncol = 200)
  g4.mat <- matrix(nrow = 0, ncol = 200)
  c4.mat <- matrix(nrow = 0, ncol = 200)
  quad.mat <- matrix(nrow = 0, ncol = 200)
  for (j in 1:length(x))  {
    #    print(j/length(x))
    A.count <- vector()
    C.count <- vector()
    G.count <- vector()
    T.count <- vector()
    GA <- vector()
    GC <- vector()
    CpG <- vector()
    g4.count <- vector()
    c4.count <- vector()
    quad.count <- vector()
    seq.split <- splitseq(x[[j]], frame = 0, 100)  # split sequences in 200 sequences of 100 bp
    for (i in 1:length(seq.split)) {
      mat <- as.matrix(count(strsplit(seq.split[i], split = "")[[1]], 1))  # compute count matrix
      p <- (mat[1] + mat[3]) / (mat[1] + mat[2] + mat[3] + mat[4])  # GA content
      q <- (mat[2] + mat[3]) / (mat[1] + mat[2] + mat[3] + mat[4])  # GC content
      A.count <- append(A.count, mat[1]) # count A
      C.count <- append(C.count, mat[2]) # count C
      G.count <- append(G.count, mat[3]) # count G
      T.count <- append(T.count, mat[4]) # count T
      GA <- append(GA, p)
      GC <- append(GC, q)
      CpG <- append(CpG, as.matrix(count(strsplit(seq.split[i], split = "")[[1]], 2))[7])
      g4.count <- append(g4.count, stringr::str_count(toupper(seq.split[i]), gquad))
      c4.count <- append(c4.count, stringr::str_count(toupper(seq.split[i]), cquad))
      quad.count <- append(quad.count, stringr::str_count(toupper(seq.split[i]), gquad) + stringr::str_count(toupper(seq.split[i]), cquad))
    }
    A.mat <- rbind(A.mat, A.count)
    C.mat <- rbind(C.mat, C.count)
    G.mat <- rbind(G.mat, G.count)
    T.mat <- rbind(T.mat, T.count)
    GA.mat <- rbind(GA.mat, GA)
    GC.mat <- rbind(GC.mat, GC)
    CpG.mat <- rbind(CpG.mat, CpG)
    g4.mat <- rbind(g4.mat, g4.count)
    c4.mat <- rbind(c4.mat, c4.count)
    quad.mat <- rbind(quad.mat, quad.count)
  }
  A.mean <- colMeans(A.mat, na.rm = TRUE)
  C.mean <- colMeans(C.mat, na.rm = TRUE)
  G.mean <- colMeans(G.mat, na.rm = TRUE)
  T.mean <- colMeans(T.mat, na.rm = TRUE)
  A.sum <- colSums(A.mat, na.rm = TRUE)
  C.sum <- colSums(C.mat, na.rm = TRUE)
  G.sum <- colSums(G.mat, na.rm = TRUE)
  T.sum <- colSums(T.mat, na.rm = TRUE)
  GA.mean <- colMeans(GA.mat, na.rm = TRUE)
  GC.mean <- colMeans(GC.mat, na.rm = TRUE)
  CpG.mean <- colMeans(CpG.mat, na.rm = TRUE)
  g4.mean <- colMeans(g4.mat, na.rm = TRUE)
  c4.mean <- colMeans(c4.mat, na.rm = TRUE)
  quad.mean <- colMeans(quad.mat, na.rm = TRUE)
  bin <- seq(-10000, 9900, 100)
  Ori.nuc.count <- cbind.data.frame(bin, A.mean, C.mean, G.mean, T.mean, A.sum, C.sum, G.sum, T.sum, GA.mean, GC.mean, CpG.mean, g4.mean, c4.mean, quad.mean)
  # Add base composition skews
  Ori.nuc.count <- Ori.nuc.count %>% mutate(GC.skew = (G.mean-C.mean)/(G.mean+C.mean), AT.skew = (A.mean-T.mean)/(A.mean+T.mean), skew = GC.skew-AT.skew)
  return(Ori.nuc.count)    # output is a dataframe returning sequence features
}

# Function to assess the goodness of the model
# RMSE (Root Mean Square Error) is the standard deviation of the residuals (prediction errors) and is expressed in the same unit of the prediction

RMSE <- function (x,y) {
  res <- (x-y)[50:150]   # consider only origins +/- 5kb
  RSS <- c(crossprod(res))
  MSE <- RSS/length(res)
  RMSE <- sqrt(MSE)
  return(RMSE)
}

############################################################################
# Compute sequence features observed at origins 

IniSeq.ori <- read.table("./Figure_1/IniSeq.center.bed")
colnames(IniSeq.ori) <- c("chr", "start", "end", "EFF")
IniSeq.ori <- IniSeq.ori %>% mutate(start = as.integer(start-10000), end = as.integer(start+20000))
# Save
write.table(IniSeq.ori, "./Figure_7/IniSeq.ori.10kb.bed", sep="\t", col.names = F, row.names = F, quote = F)
# Create grange object
IniSeq.ori.gr <- import("./Figure_7/IniSeq.ori.10kb.bed")
# Recover sequences
IniSeq.ori.seq <- getSeq(Hsapiens, IniSeq.ori.gr)
# Save multifasta
writeXStringSet(IniSeq.ori.seq, file="./Figure_7/IniSeq.ori.10kb.fasta")
# Compute features
ori.seqs <- read.fasta("./Figure_7/IniSeq.ori.10kb.fasta")
ori.features <- SeqFeat(ori.seqs)
# Save results
write.table(ori.features, "./Figure_7/IniSeq.ori.features.tsv", sep="\t", col.names = T, row.names = F, quote = F)
# Load
ori.features <- read.csv("./Figure_7/IniSeq.ori.features.tsv", sep="\t") %>% mutate(cycle = "observed")

############################################################################
# Compute sequence features of evolved sequences (in silico evolution experiment)   

# Initialise dataframe

SeqFeat.df <- tibble(bin = integer(), A.mean = integer(), C.mean = integer(), G.mean = integer(), T.mean = integer(),
                         A.sum = integer(), C.sum = integer(), G.sum = integer(), T.sum = integer(),
                         GA.mean = integer(), GC.mean = integer(), CpG.mean = integer(),  
                         g4.mean = integer(), c4.mean = integer(), quad.mean = integer(),
                         GC.skew = integer(), AT.skew = integer(), skew = integer(),
                         cycle = character())

# Select set of sequences

set <- c(0,2,4,8,16,32,64,128,256,512,1024,2048,4096)

# Compute sequence features

for (i in set) {
  print(i)
  print(Sys.time())
  infile <- paste("./Figure_7/Exp1000/", "Seqs_c_", i, ".fasta", sep = "")
  seqs <- read.fasta(infile)
  p <- SeqFeat(seqs)
  q <- p %>% mutate(cycle = i)
  SeqFeat.df <- rbind(SeqFeat.df, q)
}

write.table(SeqFeat.df, "./Figure_7/Exp1000_results.tsv", sep="\t", col.names = T, row.names = F, quote = F)

# Prepare plots

SeqFeat.df <- read.csv("./Figure_7/Exp1000_results.tsv", sep="\t")

# Prepare data for plotting (ordering cycle number)

SeqFeat.df$cycle <- as.character(SeqFeat.df$cycle)
SeqFeat.df <- SeqFeat.df %>% mutate(Cycle = fct_relevel(cycle, str_sort(unique(SeqFeat.df$cycle), numeric = T)))

# Select color palette
library(wesanderson)
palette <- wes_palette("Zissou1", length(unique(SeqFeat.df$Cycle)), type = "continuous")
# add observed values
ori.features.plot <- ori.features %>% mutate(Cycle = "observed")
SeqFeat.df <- rbind(SeqFeat.df, ori.features.plot)

# Plots reporting the feature profiles over several cycle
# together with RMSE values associated to each features (comparing observed and cycle 4096 values)

RMSE(SeqFeat.df[which(SeqFeat.df$Cycle == 4096),]$GC.mean, SeqFeat.df[which(SeqFeat.df$Cycle == "observed"),]$GC.mean) # 0.04061432
GC.plot <- SeqFeat.df %>% ggplot(aes(x=bin, y=GC.mean, color=Cycle)) +
  geom_line() +
  geom_point(cex = 0.1) +
  labs(x = "Distance from origin (bp)", y = "%GC", title = "GC content around origins | RMSE = 0.04061432") +
  scale_colour_manual(values=c(palette, "black")) +
  theme_bw() + theme(aspect.ratio=0.5)

RMSE(SeqFeat.df[which(SeqFeat.df$Cycle == 4096),]$AT.skew, SeqFeat.df[which(SeqFeat.df$Cycle == "observed"),]$AT.skew) # 0.01022134
AT.skew.plot <- SeqFeat.df %>% ggplot(aes(x=bin, y=AT.skew, color=Cycle)) +
  geom_line() +
  geom_point(cex = 0.1) + 
  labs(x = "Distance from origin (bp)", y = "AT skew", title = "RMSE = 0.01022134") +
  scale_colour_manual(values=c(palette, "black")) +
  theme_bw() + theme(aspect.ratio=0.5)

RMSE(SeqFeat.df[which(SeqFeat.df$Cycle == 4096),]$GC.skew, SeqFeat.df[which(SeqFeat.df$Cycle == "observed"),]$GC.skew) # 0.02779178
GC.skew.plot <- SeqFeat.df %>% ggplot(aes(x=bin, y=GC.skew, color=Cycle)) +
  geom_line() +
  geom_point(cex = 0.1) + xlim(-5000,5000) +
  labs(x = "Distance from origin (bp)", y = "GC skew", title = "RMSE = 0.02779178") +
  scale_colour_manual(values=c(palette, "black")) +
  theme_bw() + theme(aspect.ratio=1.1)

RMSE(SeqFeat.df[which(SeqFeat.df$Cycle == 4096),]$CpG.mean, SeqFeat.df[which(SeqFeat.df$Cycle == "observed"),]$CpG.mean) # 1.250286
CpG.plot <- SeqFeat.df %>% ggplot(aes(x=bin, y=CpG.mean, color=Cycle)) +
  geom_line() +
  geom_point(cex = 0.1) + xlim(-5000,5000) +
  labs(x = "Distance from origin (bp)", y = "CpG density", title = "RMSE = 1.250286") +
  scale_colour_manual(values=c(palette, "black")) +
  theme_bw() + theme(aspect.ratio=1.1)

RMSE(SeqFeat.df[which(SeqFeat.df$Cycle == 4096),]$quad.mean, SeqFeat.df[which(SeqFeat.df$Cycle == "observed"),]$quad.mean) # 0.04039779
quad.plot <- SeqFeat.df %>% ggplot(aes(x=bin, y=quad.mean, color=Cycle)) +
  geom_line() +
  geom_point(cex = 0.1) + xlim(-5000,5000) +
  labs(x = "Distance from origin (bp)", y = "G4 density", title = "RMSE = 0.04039779") +
  scale_colour_manual(values=c(palette, "black")) +
  theme_bw() + theme(aspect.ratio=1.1)

# Save plots

pdf("./Rplot/Fig_7B_7C.pdf", width=8, height=8, useDingbats=FALSE)
ggarrange(GC.plot, AT.skew.plot, ncol = 1, nrow = 2)
dev.off()

pdf("./Rplot/Fig_S7D_S7E_S7F.pdf", width=18, height=4, useDingbats=FALSE)
ggarrange(GC.skew.plot, CpG.plot, quad.plot, ncol = 3, nrow = 1)
dev.off()

############################################################################
############################################################################
###           Analysis of in silico evolution experiment                 ###
###                     Sequence diversity                               ###
############################################################################
############################################################################

############################################################################
# Sequence diverisity is quantified as Shannon diversity index
# analysing kmer (k = 6) frequency within origin sequences

# Define function for computing kmer occurrences

kmer.freq.calculation <- function(x) {  # x is a DNAStringSet object
  # initialise tibble on first 100 nucleotides
  kmer.freq.ori <- tibble()
  n <- length(x)
  seq.split.1 <- subseq(x, start= 1, end= 100)
  kmer.freq.1 <- oligonucleotideFrequency(seq.split.1, 6, as.prob= TRUE) # to have a frequency as probabilities
  kmer.freq.ori <- as.data.frame(colMeans(kmer.freq.1)) %>% add_rownames(var = "k-mer")
  colnames(kmer.freq.ori) <- c("k-mer", "1")
  # Compute values for all bins
  for (i in 1:199) {
    #    print(i)
    seq.split <- subseq(x, start= (i*100)+1, end= (i*100)+100)
    kmer.freq <- oligonucleotideFrequency(seq.split, 6, as.prob= TRUE)
    p <- as.data.frame(colMeans(kmer.freq)) %>% add_rownames(var = "k-mer")
    colnames(p) <- c("k-mer", i*100)
    kmer.freq.ori <- kmer.freq.ori %>% full_join(p, by = "k-mer")
  }
  return(kmer.freq.ori)
}

# Define function for computing sequence diversity (Shannon diversity index)
# The Shannon (or Shannon–Wiener) index is defined as H = -sum p_i log(b) p_i, where p_i is the proportional abundance of species i and b is the base of the logarithm.

Seq.diversity <- function(x) {  # x is a dataframe reporting the occurences of all k-mers (k=6) in a given sequences
  # in our case x is the result from the kmer.freq.calculation function
  dist <- colnames(x[,-1])
  diversity <- vector()
  for (i in dist) {
    kmer.dist <- x[,which(colnames(x) == i)]
    p <- diversity(kmer.dist, "shannon")
    diversity <- append(diversity, p)
  }
  kmer.diversity <- cbind.data.frame(dist = as.numeric(dist), diversity)
  return(kmer.diversity)
}

############################################################################
# Compute Shannon diversity indexes for IniSeq origins

ori.seqs <- readDNAStringSet("./Figure_7/IniSeq.ori.10kb.fasta")  # 9,351 origins
Ref.kmer <- kmer.freq.calculation(ori.seqs)
Ref.kmer.diversity <- Seq.diversity(Ref.kmer) %>% mutate(cycle = "observed")

############################################################################
# Compute Shannon diversity indexes for evolved sequences

SeqDiv.df <- data_frame(dist = integer(), diversity = integer(), cycle = character()) 
for (i in set) {
  print(i)
  print(Sys.time())
  infile <- paste("./Figure_7/Exp1000_first_round/", "Seqs_c_", i, ".fasta", sep = "")
  seqs <- readDNAStringSet(infile)
  p <- kmer.freq.calculation(seqs)
  q <- Seq.diversity(p) %>% mutate(cycle = i)
  SeqDiv.df <- rbind(SeqDiv.df, q)
}

SeqDiv.df <- rbind(SeqDiv.df, Ref.kmer.diversity)

write.table(SeqDiv.df, "./Figure_7/Exp1000_seq_diversity.tsv", sep="\t", col.names = T, row.names = F, quote = F)

# Load results

SeqDiv.df.plot <- read.csv("./Figure_7/Exp1000_seq_diversity.tsv", sep="\t")

# Prepare data for plotting (ordering cycle number)

SeqDiv.df.plot$cycle <- as.character(SeqDiv.df.plot$cycle)
SeqDiv.df.plot <- SeqDiv.df.plot %>% mutate(dist = dist - 10000) %>% mutate(Cycle = fct_relevel(cycle, str_sort(unique(SeqDiv.df.plot$cycle), numeric = T)))

# Correct observed values to account for bigger sample size
# There is 9351/1000=9.351 more sequences
# Need to add values/(9.351*ln(4096)) to the computed values

SeqDiv.df.plot[which(SeqDiv.df.plot$cycle == "observed"),]$diversity <- SeqDiv.df.plot[which(SeqDiv.df.plot$cycle == "observed"),]$diversity+(SeqDiv.df.plot[which(SeqDiv.df.plot$cycle == "observed"),]$diversity/(9.351*log(4096)))

palette.4 <- wes_palette("Zissou1", length(unique(as.character(SeqDiv.df.plot$Cycle)))-1, type = "continuous")

RMSE(SeqDiv.df.plot[which(SeqDiv.df.plot$Cycle == 4096),]$diversity, SeqDiv.df.plot[which(SeqDiv.df.plot$Cycle == "observed"),]$diversity) # 0.07124831
div.plot <- SeqDiv.df.plot %>% ggplot(aes(x=dist, y=diversity, color=Cycle)) +
  geom_line() + 
  geom_point(cex = 0.1) +
  labs(x = "Distance from origin (bp)", y = "Shannon diversity index", title = "k-mers (k = 6) diversity around origins | RMSE = 0.07124831") +
  scale_colour_manual(values=c(palette.4, "black")) +
  theme_bw() + theme(aspect.ratio=0.5)

pdf("./Rplot/Fig_7D.pdf", width=8, height=4, useDingbats=FALSE)
div.plot
dev.off()

############################################################################
# Compute Shannon diversity indexes for IniSeq origins when binned by efficiency

# Compute kmer occurrences for origins when binned by efficiency

quantile(IniSeq.ori$EFF, c(0.333,0.666))
#    33.3%  66.6% 
#   0.8775 0.9316 

IniSeq.ori.low <- IniSeq.ori %>% filter(EFF < 0.8775) # 3,111 origins
IniSeq.ori.medium <- IniSeq.ori %>% filter(EFF >= 0.8775 & EFF < 0.9316) # 3,115 origins
IniSeq.ori.high <- IniSeq.ori %>% filter(EFF >= 0.9316) # 3,125 origins
# Save
write.table(IniSeq.ori.low, "./Figure_7/IniSeq.ori.low.10kb.bed", sep="\t", col.names = F, row.names = F, quote = F)
write.table(IniSeq.ori.medium, "./Figure_7/IniSeq.ori.medium.10kb.bed", sep="\t", col.names = F, row.names = F, quote = F)
write.table(IniSeq.ori.high, "./Figure_7/IniSeq.ori.high.10kb.bed", sep="\t", col.names = F, row.names = F, quote = F)
# Create grange object
IniSeq.ori.low.gr <- import("./Figure_7/IniSeq.ori.low.10kb.bed")
IniSeq.ori.medium.gr <- import("./Figure_7/IniSeq.ori.medium.10kb.bed")
IniSeq.ori.high.gr <- import("./Figure_7/IniSeq.ori.high.10kb.bed")
# Recover sequences
IniSeq.ori.low.seq <- getSeq(Hsapiens, IniSeq.ori.low.gr)
IniSeq.ori.medium.seq <- getSeq(Hsapiens, IniSeq.ori.medium.gr)
IniSeq.ori.high.seq <- getSeq(Hsapiens, IniSeq.ori.high.gr)
# Save multifasta
writeXStringSet(IniSeq.ori.low.seq, file="./Figure_7/IniSeq.ori.low.10kb.fasta")
writeXStringSet(IniSeq.ori.medium.seq, file="./Figure_7/IniSeq.ori.medium.10kb.fasta")
writeXStringSet(IniSeq.ori.high.seq, file="./Figure_7/IniSeq.ori.high.10kb.fasta")

# Compute kmer occurrences for origins when binned by efficiency

ori.low.seqs <- readDNAStringSet("./Figure_7/IniSeq.ori.low.10kb.fasta")
ori.medium.seqs <- readDNAStringSet("./Figure_7/IniSeq.ori.medium.10kb.fasta")
ori.high.seqs <- readDNAStringSet("./Figure_7/IniSeq.ori.high.10kb.fasta")

Ref.kmer.low <- kmer.freq.calculation(ori.low.seqs)
Ref.kmer.medium <- kmer.freq.calculation(ori.medium.seqs)
Ref.kmer.high <- kmer.freq.calculation(ori.high.seqs)

# Compute Shannon diversity index

Ref.kmer.low.diversity <- Seq.diversity(Ref.kmer.low) %>% mutate(dist = dist - 10000) %>% mutate(class = "low")
Ref.kmer.medium.diversity <- Seq.diversity(Ref.kmer.medium) %>% mutate(dist = dist - 10000) %>% mutate(class = "medium")
Ref.kmer.high.diversity <- Seq.diversity(Ref.kmer.high) %>% mutate(dist = dist - 10000) %>% mutate(class = "high")

# Prepare plot

Ref.kmer.diversity.plot <- rbind(Ref.kmer.low.diversity, Ref.kmer.medium.diversity, Ref.kmer.high.diversity)

pdf("./Rplot/Fige_7E.pdf", width=8, height=4, useDingbats=FALSE)
Ref.kmer.diversity.plot %>% ggplot(aes(x = dist, y = diversity, color = class)) +
  geom_line(size = 0.5) + xlim(-10000,10000) +
  scale_color_manual(values=c("#FC4E07", "#0072B2", "#CC79A7")) +
  xlab("Distance from origin (bp)") + ylab("Shannon diversity index") + ggtitle("Sequence diversity (kmers, k = 6) around origins") +
  theme_bw() + theme(aspect.ratio=0.5)
dev.off()

# Compute statistics

Ref.kmer.diversity.stat <- Ref.kmer.diversity.plot %>% filter((dist <= -500 & dist >= -2500) | (dist >= 500 & dist <= 2500)) %>% mutate(round = round(diversity, digits = 2))
test.table <- t(table(Ref.kmer.diversity.stat$class, Ref.kmer.diversity.stat$round))
test <- chisq.test(test.table)
test  # p-value = 0.0008263

############################################################################
############################################################################
###              Nucleotide substitution rate estimation                 ###
###                            GERP score                                ###
############################################################################
############################################################################

# Genomic evolutionary rate profiling (GERP) scores 
# were used as a measure of nucleotide diversity between species.
# The single-nucleotide resolution bigWig file (gerp_conservation_scores.homo_sapiens.GRCh38.bw)
# reporting substitution rates computed from the multiple alignment of 111 mammal genomes
# was obtained from ftp://ftp.ensembl.org/pub/release101/compara/conservation_scores/111_mammals.gerp_conservation_score/.

####################################################################################
# Load genomic coordinates of origins

SNS.core <- read.table("./Figure_1/SNS.core.inter.center.10kb.bed")
SNS.stochastic <- read.table("./Figure_1/SNS.stochastic.inter.center.10kb.bed")
IniSeq <- read.table("./Figure_1/IniSeq.inter.center.10kb.bed")
colnames(SNS.core) <- c("chr", "start", "end", "INFO")
colnames(SNS.stochastic) <- c("chr", "start", "end", "INFO")
colnames(IniSeq) <- c("chr", "start", "end", "INFO")

# Remove chr from chromosome name

SNS.core <- SNS.core %>% mutate(chr = str_remove(chr, "chr"))
SNS.stochastic <- SNS.stochastic %>% mutate(chr = str_remove(chr, "chr"))
IniSeq <- IniSeq %>% mutate(chr = str_remove(chr, "chr"))

# Subset Iniseq origins by efficiency

# Define low, medium and high efficiency values as before
#  33.3% 66.6% 
#  0.8775 0.9316 

IniSeq.low <- IniSeq %>% filter(INFO < 0.8775)
IniSeq.medium <- IniSeq %>% filter(INFO >= 0.8775 & INFO < 0.9316)
IniSeq.high <- IniSeq %>% filter(INFO >= 0.9316)

####################################################################################
# Compute mean of GERP scores in rolling windows around origins

# all IniSeq origins

IniSeq.GERP.score <- rep(0, 19902)  # number of values computed in rolling windows
for (i in 1:nrow(IniSeq)) {
  print(i/nrow(IniSeq))
  ori <- IniSeq[i,]
  # Recover GERP scores
  GERP.ori <- importScore("./Dataset/gerp_conservation_scores.homo_sapiens.GRCh38.bw", format= "BigWig", ranges=GRanges(ori$chr, IRanges(ori$start, ori$end)))
  # Create matrix
  mat <- matrix(nrow = 20001, ncol = 0) %>% cbind.data.frame(pos = c(ori$start:ori$end)) %>%
    left_join(cbind.data.frame(pos = start(GERP.ori$dat), score = GERP.ori$dat$score), by = "pos")
  # Add 0 as first and last values
  mat$score[1] <- 0
  mat$score[20001] <- 0
  # Fill empty values
  mat.fill <- na.locf(na.locf(mat), fromLast = TRUE)
  # Compute values in rolling windows
  p <- rollmean(mat.fill$score, 100)  # rolling window of 100 bp
  IniSeq.GERP.score <- IniSeq.GERP.score + p
}
# Prepare final table
IniSeq.GERP.score.result <- cbind.data.frame(dist = c(-9950:9951), GERP.score = IniSeq.GERP.score/nrow(IniSeq), class = "IniSeq")
plot(IniSeq.GERP.score.result$dist, IniSeq.GERP.score.result$GERP.score)
# Correct signal from background level
IniSeq.GERP.score.baseline <- IniSeq.GERP.score.result %>% filter(dist <= -8000 | dist >= 8000)
IniSeq.GERP.score.result.corrected <- IniSeq.GERP.score.result %>% mutate(GERP.score.corr = -log2(GERP.score/mean(IniSeq.GERP.score.baseline$GERP.score)))
# Save result
write.table(IniSeq.GERP.score.result.corrected, "./Figure_7/GERP_score_IniSeq_all.tsv", quote=FALSE, sep='\t', row.names = F, col.names = T)

# SNS core origins

SNS.core.GERP.score <- rep(0, 19902)  # number of values computed in rolling windows
for (i in 1:nrow(SNS.core)) {
  print(i/nrow(SNS.core))
  ori <- SNS.core[i,]
  # Recover GERP scores
  GERP.ori <- importScore("./Dataset/gerp_conservation_scores.homo_sapiens.GRCh38.bw", format= "BigWig", ranges=GRanges(ori$chr, IRanges(ori$start, ori$end)))
  # Create matrix
  mat <- matrix(nrow = 20001, ncol = 0) %>% cbind.data.frame(pos = c(ori$start:ori$end)) %>%
    left_join(cbind.data.frame(pos = start(GERP.ori$dat), score = GERP.ori$dat$score), by = "pos")
  # Add 0 as first and last values
  mat$score[1] <- 0
  mat$score[20001] <- 0
  # Fill empty values
  mat.fill <- na.locf(na.locf(mat), fromLast = TRUE)
  # Compute values in rolling windows
  p <- rollmean(mat.fill$score, 100)  # rolling window of 100 bp
  SNS.core.GERP.score <- SNS.core.GERP.score + p
}
# Prepare final table
SNS.core.GERP.score.result <- cbind.data.frame(dist = c(-9950:9951), GERP.score = SNS.core.GERP.score/nrow(SNS.core), class = "SNS core")
plot(SNS.core.GERP.score.result$dist, SNS.core.GERP.score.result$GERP.score)
# Correct signal from background level
SNS.core.GERP.score.baseline <- SNS.core.GERP.score.result %>% filter(dist <= -8000 | dist >= 8000)
SNS.core.GERP.score.result.corrected <- SNS.core.GERP.score.result %>% mutate(GERP.score.corr = -log2(GERP.score/mean(SNS.core.GERP.score.baseline$GERP.score)))
# Save results
write.table(SNS.core.GERP.score.result.corrected, "./Figure_7/GERP_score_SNS.core.tsv", quote=FALSE, sep='\t', row.names = F, col.names = T)

# SNS stochastic origins

SNS.stochastic.GERP.score <- rep(0, 19902)  # number of values computed in rolling windows
for (i in 1:nrow(SNS.stochastic)) {
  print(i/nrow(SNS.stochastic))
  ori <- SNS.stochastic[i,]
  # Recover GERP scores
  GERP.ori <- importScore("./Dataset/gerp_conservation_scores.homo_sapiens.GRCh38.bw", format= "BigWig", ranges=GRanges(ori$chr, IRanges(ori$start, ori$end)))
  # Create matrix
  mat <- matrix(nrow = 20001, ncol = 0) %>% cbind.data.frame(pos = c(ori$start:ori$end)) %>%
    left_join(cbind.data.frame(pos = start(GERP.ori$dat), score = GERP.ori$dat$score), by = "pos")
  # Add 0 as first and last values
  mat$score[1] <- 0
  mat$score[20001] <- 0
  # Fill empty values
  mat.fill <- na.locf(na.locf(mat), fromLast = TRUE)
  # Compute values in rolling windows
  p <- rollmean(mat.fill$score, 100)  # rolling window of 100 bp
  SNS.stochastic.GERP.score <- SNS.stochastic.GERP.score + p
}
# Prepare final table
SNS.stochastic.GERP.score.result <- cbind.data.frame(dist = c(-9950:9951), GERP.score = SNS.stochastic.GERP.score/nrow(SNS.stochastic), class = "SNS stochastic")
plot(SNS.stochastic.GERP.score.result$dist, SNS.stochastic.GERP.score.result$GERP.score)
# Correct signal from background level
SNS.stochastic.GERP.score.baseline <- SNS.stochastic.GERP.score.result %>% filter(dist <= -8000 | dist >= 8000)
SNS.stochastic.GERP.score.result.corrected <- SNS.stochastic.GERP.score.result %>% mutate(GERP.score.corr = -log2(GERP.score/mean(SNS.stochastic.GERP.score.baseline$GERP.score)))
# Save results
write.table(SNS.stochastic.GERP.score.result.corrected, "./Figure_7/GERP_score_SNS.stochastic.tsv", quote=FALSE, sep='\t', row.names = F, col.names = T)

# IniSeq origins of low efficiency

IniSeq.low.GERP.score <- rep(0, 19902)  # number of values computed in rolling windows
for (i in 1:nrow(IniSeq.low)) {
  print(i/nrow(IniSeq.low))
  ori <- IniSeq.low[i,]
  # Recover GERP scores
  GERP.ori <- importScore("./Dataset/gerp_conservation_scores.homo_sapiens.GRCh38.bw", format= "BigWig", ranges=GRanges(ori$chr, IRanges(ori$start, ori$end)))
  # Create matrix
  mat <- matrix(nrow = 20001, ncol = 0) %>% cbind.data.frame(pos = c(ori$start:ori$end)) %>%
    left_join(cbind.data.frame(pos = start(GERP.ori$dat), score = GERP.ori$dat$score), by = "pos")
  # Add 0 as first and last values
  mat$score[1] <- 0
  mat$score[20001] <- 0
  # Fill empty values
  mat.fill <- na.locf(na.locf(mat), fromLast = TRUE)
  # Compute values in rolling windows
  p <- rollmean(mat.fill$score, 100)  # rolling window of 100 bp
  IniSeq.low.GERP.score <- IniSeq.low.GERP.score + p
}
# Prepare final table
IniSeq.low.GERP.score.result <- cbind.data.frame(dist = c(-9950:9951), GERP.score = IniSeq.low.GERP.score/nrow(IniSeq.low), class = "low")
plot(IniSeq.low.GERP.score.result$dist, IniSeq.low.GERP.score.result$GERP.score)
# Correct signal from background level
IniSeq.low.GERP.score.baseline <- IniSeq.low.GERP.score.result %>% filter(dist <= -8000 | dist >= 8000)
IniSeq.low.GERP.score.result.corrected <- IniSeq.low.GERP.score.result %>% mutate(GERP.score.corr = -log2(GERP.score/mean(IniSeq.low.GERP.score.baseline$GERP.score)))
# Save result
write.table(IniSeq.low.GERP.score.result.corrected, "./Figure_7/GERP_score_IniSeq_low.tsv", quote=FALSE, sep='\t', row.names = F, col.names = T)

# IniSeq origins of medium efficiency

IniSeq.medium.GERP.score <- rep(0, 19902)  # number of values computed in rolling windows
for (i in 1:nrow(IniSeq.medium)) {
  print(i/nrow(IniSeq.medium))
  ori <- IniSeq.medium[i,]
  # Recover GERP scores
  GERP.ori <- importScore("./Dataset/gerp_conservation_scores.homo_sapiens.GRCh38.bw", format= "BigWig", ranges=GRanges(ori$chr, IRanges(ori$start, ori$end)))
  # Create matrix
  mat <- matrix(nrow = 20001, ncol = 0) %>% cbind.data.frame(pos = c(ori$start:ori$end)) %>%
    left_join(cbind.data.frame(pos = start(GERP.ori$dat), score = GERP.ori$dat$score), by = "pos")
  # Add 0 as first and last values
  mat$score[1] <- 0
  mat$score[20001] <- 0
  # Fill empty values
  mat.fill <- na.locf(na.locf(mat), fromLast = TRUE)
  # Compute values in rolling windows
  p <- rollmean(mat.fill$score, 100)  # rolling window of 100 bp
  IniSeq.medium.GERP.score <- IniSeq.medium.GERP.score + p
}
# Prepare final table
IniSeq.medium.GERP.score.result <- cbind.data.frame(dist = c(-9950:9951), GERP.score = IniSeq.medium.GERP.score/nrow(IniSeq.medium), class = "medium")
plot(IniSeq.medium.GERP.score.result$dist, IniSeq.medium.GERP.score.result$GERP.score)
# Correct signal from background level
IniSeq.medium.GERP.score.baseline <- IniSeq.medium.GERP.score.result %>% filter(dist <= -8000 | dist >= 8000)
IniSeq.medium.GERP.score.result.corrected <- IniSeq.medium.GERP.score.result %>% mutate(GERP.score.corr = -log2(GERP.score/mean(IniSeq.medium.GERP.score.baseline$GERP.score)))
# Save result
write.table(IniSeq.medium.GERP.score.result.corrected, "./Figure_7/GERP_score_IniSeq_medium.tsv", quote=FALSE, sep='\t', row.names = F, col.names = T)

# IniSeq origins of high efficiency

IniSeq.high.GERP.score <- rep(0, 19902)  # number of values computed in rolling windows
for (i in 1:nrow(IniSeq.high)) {
  print(i/nrow(IniSeq.high))
  ori <- IniSeq.high[i,]
  # Recover GERP scores
  GERP.ori <- importScore("./Dataset/gerp_conservation_scores.homo_sapiens.GRCh38.bw", format= "BigWig", ranges=GRanges(ori$chr, IRanges(ori$start, ori$end)))
  # Create matrix
  mat <- matrix(nrow = 20001, ncol = 0) %>% cbind.data.frame(pos = c(ori$start:ori$end)) %>%
    left_join(cbind.data.frame(pos = start(GERP.ori$dat), score = GERP.ori$dat$score), by = "pos")
  # Add 0 as first and last values
  mat$score[1] <- 0
  mat$score[20001] <- 0
  # Fill empty values
  mat.fill <- na.locf(na.locf(mat), fromLast = TRUE)
  # Compute values in rolling windows
  p <- rollmean(mat.fill$score, 100)  # rolling window of 100 bp
  IniSeq.high.GERP.score <- IniSeq.high.GERP.score + p
}
# Prepare final table
IniSeq.high.GERP.score.result <- cbind.data.frame(dist = c(-9950:9951), GERP.score = IniSeq.high.GERP.score/nrow(IniSeq.high), class = "high")
plot(IniSeq.high.GERP.score.result$dist, IniSeq.high.GERP.score.result$GERP.score)
# Correct signal from background level
IniSeq.high.GERP.score.baseline <- IniSeq.high.GERP.score.result %>% filter(dist <= -8000 | dist >= 8000)
IniSeq.high.GERP.score.result.corrected <- IniSeq.high.GERP.score.result %>% mutate(GERP.score.corr = -log2(GERP.score/mean(IniSeq.high.GERP.score.baseline$GERP.score)))
# Save result
write.table(IniSeq.high.GERP.score.result.corrected, "./Figure_7/GERP_score_IniSeq_high.tsv", quote=FALSE, sep='\t', row.names = F, col.names = T)

####################################################################################
# Prepare plots

# Load data
IniSeq.GERP.score.result.corrected <- read.csv("./Figure_7/GERP_score_IniSeq_all.tsv", sep='\t')
SNS.core.GERP.score.result.corrected <- read.csv("./Figure_7/GERP_score_SNS.core.tsv", sep='\t')
SNS.stochastic.GERP.score.result.corrected <- read.csv("./Figure_7/GERP_score_SNS.stochastic.tsv", sep='\t')
IniSeq.low.GERP.score.result.corrected <- read.csv("./Figure_7/GERP_score_IniSeq_low.tsv", sep='\t')
IniSeq.medium.GERP.score.result.corrected <- read.csv("./Figure_7/GERP_score_IniSeq_medium.tsv", sep='\t')
IniSeq.high.GERP.score.result.corrected <- read.csv("./Figure_7/GERP_score_IniSeq_high.tsv", sep='\t')

# Plot
IniSeq.GERP <- rbind(IniSeq.low.GERP.score.result.corrected, IniSeq.medium.GERP.score.result.corrected, IniSeq.high.GERP.score.result.corrected)
Ori.GERP <- rbind(IniSeq.GERP.score.result.corrected, SNS.core.GERP.score.result.corrected, SNS.stochastic.GERP.score.result.corrected)

IniSeq.GERP.plot <- IniSeq.GERP %>% ggplot(aes(x = dist, y = GERP.score.corr, color = class)) +
  geom_line(size = 0.5) + 
  scale_color_manual(values=c("#FC4E07", "#0072B2", "#CC79A7")) +
  xlab("Distance from origin (bp)") + ylab("GERP score (log2 FC)") + ggtitle("GERP score variation around origins") +
  theme_bw() + theme(aspect.ratio=0.5) + scale_y_reverse()

Ori.GERP.plot <- Ori.GERP %>% ggplot(aes(x = dist, y = GERP.score.corr, color = class)) +
  geom_line(size = 0.5) +
  scale_color_manual(values=c("#FC4E07", "#009E73", "#0072B2")) +
  xlab("Distance from origin (bp)") + ylab("GERP score (log2 FC)") + ggtitle("GERP score variation around origins") +
  theme_bw() + theme(aspect.ratio=0.5) + scale_y_reverse()

pdf("./Rplot/Fig_7E.pdf", width=8, height=4, useDingbats=FALSE)
IniSeq.GERP.plot
dev.off()

pdf("./Rplot/Fig_S7G.pdf", width=8, height=4, useDingbats=FALSE)
Ori.GERP.plot
dev.off()

# Compute statistics

IniSeq.low.flank <- IniSeq.low.GERP.score.result.corrected %>% filter((dist <= -500 & dist >= -2500) | (dist >= 500 & dist <= 2500))
IniSeq.medium.flank <- IniSeq.medium.GERP.score.result.corrected %>% filter((dist <= -500 & dist >= -2500) | (dist >= 500 & dist <= 2500))
IniSeq.high.flank <- IniSeq.high.GERP.score.result.corrected %>% filter((dist <= -500 & dist >= -2500) | (dist >= 500 & dist <= 2500))
test.df <- rbind(IniSeq.low.flank, IniSeq.medium.flank, IniSeq.high.flank) %>% mutate(round = round(GERP.score.corr, digits = 2))
boxplot(GERP.score.corr ~ class,data = test.df)
test <- chisq.test(t(table(test.df$class, test.df$round)))
test  # p-value < 2.2e-16

IniSeq.flank <- IniSeq.GERP.score.result.corrected %>% filter((dist <= -500 & dist >= -2500) | (dist >= 500 & dist <= 2500))
SNS.core.flank <- SNS.core.GERP.score.result.corrected %>% filter((dist <= -500 & dist >= -2500) | (dist >= 500 & dist <= 2500))
SNS.stochastic.flank <- SNS.stochastic.GERP.score.result.corrected %>% filter((dist <= -500 & dist >= -2500) | (dist >= 500 & dist <= 2500))
test.df <- rbind(IniSeq.flank, SNS.core.flank, SNS.stochastic.flank) %>% mutate(round = round(GERP.score.corr, digits = 2))
boxplot(GERP.score.corr ~ class,data = test.df)
test <- chisq.test(t(table(test.df$class, test.df$round)))
test  # p-value < 2.2e-16









