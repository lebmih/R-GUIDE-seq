# Working folder ----------------------------------------------------------
# Input your working folder path here
setwd("__")


# Libraries ---------------------------------------------------------------
library(stringr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggridges)
library(seqinr)
library(stringi)
library(DescTools)
library(Biostrings)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
library(BSgenome.Hsapiens.UCSC.hg38)


# Alignment ---------------------------------------------------------------
clones.symmetry.simple <- read.table("__SymmetryReadsSimplified.txt",
                                     header = TRUE, sep = "\t", dec = ".",
                                     stringsAsFactors = FALSE)
# Specify the sgRNA sequences here
sgrna1 <- "TCAGGACCTGCGCTGCTGAG*GG"
sgrna2 <- "TCCACTCCTGATGCTGTTAT*GG"
# Specify the number of rows in the plot
number.of.seqs.to.visualize <- 30
# Specify the parameter you would like to use to sort the plot by
clones.symmetry.simple.tmp <- clones.symmetry.simple %>%
  arrange(min.ed) %>%
  filter(sgrna == "sgrna1")

clones.symmetry.simple.subset <- head(clones.symmetry.simple.tmp, number.of.seqs.to.visualize)
clones.symmetry.simple.subset <- clones.symmetry.simple.subset %>%
  mutate(seq.id = seq(number.of.seqs.to.visualize))

seqs <- clones.symmetry.simple.subset$target.seq

seqs <- c(sgrna1,seqs)

seq.df <- data.frame(seq.id = seq(0,length(seqs)-1),
                     seq = seqs)
seq.df <- seq.df %>%
  separate(seq, into = as.character(seq(0,str_length(seqs[1]))), sep = "") %>%
  select(-"0") %>%
  pivot_longer("1":"23", names_to = "pos", values_to = "nt")

seq.df <- seq.df %>%
  left_join(seq.df %>%
              filter(seq.id == 0) %>%
              select(-seq.id), by = "pos") %>%
  mutate(nt = ifelse(seq.id == 0, nt.x,
                     ifelse(nt.x == nt.y, ".",
                            nt.x))) %>%
  select(-nt.x,-nt.y) %>%
  mutate(seq.id = as.character(seq.id))

seq.df <- seq.df %>%
  left_join(clones.symmetry.simple.subset %>%
              select(chrom, cut.pos, ctrl, trt.wo.inh, trt.w.inh, gene, seq.id) %>%
              mutate(seq.id = as.character(seq.id)),
            by = "seq.id")


clones.symmetry.simple.subset <- clones.symmetry.simple.subset %>%
  separate(gene, into = "gene", sep = "\\|")

ggplot(seq.df)+
  geom_tile(aes(x = pos, y = seq.id, fill = nt,
                color = nt, alpha = nt, linewidth = nt))+
  geom_text(aes(x = pos, y = seq.id, label = nt),
            size = 6)+
  geom_text(data = clones.symmetry.simple.subset,
            aes(y = rev(seq.id), label = chrom), x = 24, hjust = 0)+
  geom_text(data = clones.symmetry.simple.subset,
            aes(y = rev(seq.id), label = cut.pos), x = 26, hjust = 0)+
  geom_text(data = clones.symmetry.simple.subset,
            aes(y = rev(seq.id), label = gene), x = 30, hjust = 0)+
  geom_text(data = clones.symmetry.simple.subset,
              aes(y = rev(seq.id), label = trt.wo.inh), x = 33, hjust = 0)+
  geom_text(data = clones.symmetry.simple.subset,
            aes(y = rev(seq.id), label = trt.w.inh), x = 35, hjust = 0)+
  geom_text(data = clones.symmetry.simple.subset,
            aes(y = rev(seq.id), label = ctrl), x = 37, hjust = 0)+
  scale_color_manual(values = c("A" = "black",
                                "C" = "black",
                                "T" = "black",
                                "G" = "black",
                                "N" = "black",
                                "*" = "black",
                                "." = "white"))+
  scale_linewidth_manual(values = c("A" = 0.3,
                                    "C" = 0.3,
                                    "T" = 0.3,
                                    "G" = 0.3,
                                    "N" = 0.3,
                                    "*" = 0.3,
                                    "." = 0))+
  scale_alpha_manual(values = c("A" = 1,
                                "C" = 1,
                                "T" = 1,
                                "G" = 1,
                                "N" = 1,
                                "*" = 1,
                                "." = 0))+
  scale_fill_manual(values = c("A" = "#ff6500",
                               "C" = "#00a1ff",
                               "T" = "#00c000",
                               "G" = "#faff00",
                               "N" = "#c2c2c2",
                               "*" = "#c2c2c2",
                               "." = "white"))+
  scale_y_discrete(name = "Target",
                   limits = as.character(rev(seq(0,number.of.seqs.to.visualize))),
                   labels = as.character(rev(c("", seq(number.of.seqs.to.visualize)))))+
  annotate("text", x = 24, y = number.of.seqs.to.visualize+1, label = "chrom", hjust = 0)+
  annotate("text", x = 26, y = number.of.seqs.to.visualize+1, label = "cut, nt", hjust = 0)+
  annotate("text", x = 30, y = number.of.seqs.to.visualize+1, label = "gene", hjust = 0)+
  annotate("text", x = 33, y = number.of.seqs.to.visualize+1, label = "reads,\ntreated", hjust = 0)+
  annotate("text", x = 35, y = number.of.seqs.to.visualize+1, label = "reads,\n+inh", hjust = 0)+
  annotate("text", x = 37, y = number.of.seqs.to.visualize+1, label = "reads,\ncontrol", hjust = 0)+
  scale_x_discrete(name = "Position, nt",
                   limits = as.character(seq(1,39)),
                   labels = as.character(c(seq(1,23), rep("",16))))+
  theme_classic()+
  theme(legend.position = "none",
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_text(hjust=0.3))

ggsave("___SortedBy___.svg",
       width = 11, height = 0.3*length(seqs), device = "svg")
