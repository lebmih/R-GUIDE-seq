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

# Functions ---------------------------------------------------------------
# This function aids the search for the particular sequence via text search,
# creating a library of regular expressions matching to the pattern with
# specified number of mutations.
create.mutants.lib <- function(sequence, mut.n = 1, if.collapse = TRUE){
  mut.lib <- c()
  for (i in seq(str_length(sequence))) {
    sequence.tmp <- sequence
    substr(sequence.tmp, i, i) <- "."
    mut.lib <- c(mut.lib, sequence.tmp)
  }
  if(if.collapse == TRUE) {
    mut.lib <- str_c(mut.lib, collapse = "|")
  }
  return(mut.lib)
}
# This function is similar to the previous one but allows for indels.
create.mutants.lib.indel <- function(sequence, mut.n = 1, if.collapse = TRUE){
  mut.lib <- c()
  for (i in seq(str_length(sequence))) {
    sequence.tmp <- paste0(substr(sequence,0,i-1), ".{0,2}",
                           substr(sequence,i+1, str_length(sequence)))
    mut.lib <- c(mut.lib, sequence.tmp)
  }
  if(if.collapse == TRUE) {
    mut.lib <- str_c(mut.lib, collapse = "|")
  }
  return(mut.lib)
}
# A simple reverse-complement transformation.
rc <- function(string) {
  string <- stri_reverse(string)
  string <- unlist(base::strsplit(string, ""))
  a <- which(string == "A")
  t <- which(string == "T")
  c <- which(string == "C")
  g <- which(string == "G")
  string[a] <- "T"
  string[t] <- "A"
  string[c] <- "G"
  string[g] <- "C"
  string <- stri_flatten(string)
  return(string)
}

# Import files ------------------------------------------------------------
# The script uses PEAR-assembled reads translated into FASTA format as input.
files <- list.files(path = "reads/assembled/", pattern = ".fasta", full.names = TRUE)

reads <- data.frame()
for (i in files) {
  reads.fa <- read.fasta(i, as.string = TRUE)
  filename <- i %>%
    str_split("\\/") %>% unlist() %>%
    last() %>% str_remove_all("\\.\\w*")
  sample <- filename %>% str_split("_") %>% unlist()
  sample <- sample[4]
  reads.tmp <- data.frame(read.id = names(unlist(reads.fa)),
                          seq = str_to_upper(unlist(reads.fa)),
                          sample = sample,
                          row.names = NULL,
                          stringsAsFactors = FALSE)
  reads <- reads %>% bind_rows(reads.tmp)
  print(i)
}
rm(reads.fa, reads.tmp, i, filename, sample)

# Specify the dsODN sequence here
dsodn.full <- "GTTTAATTGAGTTGTCATATGTTAATAACGGTAT"
names(dsodn.full) <- "dsODN"
dsodn.rc.full <- "ATACCGTTATTAACATATGACAACTCAATTAAAC"
names(dsodn.rc.full) <- "dsODN.RC"
# Specify the gene-specific primers here
gsp3 <- "GTTTAATTGAGTTGTCATATGTTAATAACGGT"
names(gsp3) <- "gsp3"
gsp5 <- "ATACCGTTATTAACATATGACAACTCAATTAA"
names(gsp5) <- "gsp5"


# Detecting the dsODN -----------------------------------------------------
coi <- str_sub(reads$seq, 1, 40)










reads$dsodn <- NA
reads$dsodn.start <- NA
reads$dsodn.end <- NA
# First round - detecting the primer by a simple text search.
for (i in seq(length(dsodn.full))) {
  primer.search <- unlist(str_locate(coi, dsodn.full[i]))
  whch <- !is.na(primer.search[,"start"])
  reads$dsodn[whch] <- names(dsodn.full[i])
  reads$dsodn.start[whch] <- primer.search[whch,"start"]
  reads$dsodn.end[whch] <- primer.search[whch,"end"]
  print(dsodn.full[i])
}
rm(primer.search, i)
# V primer was found in
sum(!is.na(reads$dsodn)) # During the first round.


reads$dsodn.rc <- NA
reads$dsodn.rc.start <- NA
reads$dsodn.rc.end <- NA
for (i in seq(length(dsodn.rc.full))) {
  primer.search <- unlist(str_locate(coi, dsodn.rc.full[i]))
  whch <- !is.na(primer.search[,"start"])
  reads$dsodn.rc[whch] <- names(dsodn.rc.full[i])
  reads$dsodn.rc.start[whch] <- primer.search[whch,"start"]
  reads$dsodn.rc.end[whch] <- primer.search[whch,"end"]
  print(dsodn.rc.full[i])
}
rm(primer.search, i)
sum(!is.na(reads$dsodn.rc)) 

# Second round - looking for a primer allowing for one mutation.
reads.no.dsodn <- reads %>% filter(is.na(dsodn))
coi <- str_sub(reads.no.dsodn$seq, 1, 40)
for (i in seq(length(dsodn.full))) {
  dsodn.i <- create.mutants.lib(dsodn.full[i])
  primer.search <- unlist(str_locate(coi, dsodn.i))
  whch <- !is.na(primer.search[,"start"])
  reads.no.dsodn$dsodn[whch] <- names(dsodn.full[i])
  reads.no.dsodn$dsodn.start[whch] <- primer.search[whch,"start"]
  reads.no.dsodn$dsodn.end[whch] <- primer.search[whch,"end"]
  print(dsodn.full[i])
}
# dsODN was found in
sum(!is.na(reads.no.dsodn$dsodn)) # During the second round.
# Gluing back the dataframe.
reads <- reads %>% filter(!is.na(dsodn)) %>%
  bind_rows(reads.no.dsodn)
sum(!is.na(reads$dsodn))
rm(reads.no.dsodn)

reads.no.dsodn.rc <- reads %>% filter(is.na(dsodn.rc))
coi <- str_sub(reads.no.dsodn.rc$seq, 1, 40)
for (i in seq(length(dsodn.rc.full))) {
  dsodn.rc.i <- create.mutants.lib(dsodn.rc.full[i])
  primer.search <- unlist(str_locate(coi, dsodn.rc.i))
  whch <- !is.na(primer.search[,"start"])
  reads.no.dsodn.rc$dsodn.rc[whch] <- names(dsodn.rc.full[i])
  reads.no.dsodn.rc$dsodn.rc.start[whch] <- primer.search[whch,"start"]
  reads.no.dsodn.rc$dsodn.rc.end[whch] <- primer.search[whch,"end"]
  print(dsodn.rc.full[i])
}
sum(!is.na(reads.no.dsodn.rc$dsodn.rc)) # During the second round.
# Gluing back the dataframe.
reads <- reads %>% filter(!is.na(dsodn.rc)) %>%
  bind_rows(reads.no.dsodn.rc)
sum(!is.na(reads$dsodn.rc))



rm(reads.no.dsodn.rc, primer.search, coi, dsodn.full, dsodn.i, dsodn.rc.full,
   dsodn.rc.i, i, read, sample, whch)



write.table(reads, "_____.txt",
            quote = FALSE, sep = "\t", dec = ".", row.names = FALSE)



# Filtering the reads with the ODN found ----------------------------------
# Specify the samples amplified with GSP5 or GSP3 here
reads.gsp5 <- reads %>% filter(sample %in% paste0("S", seq(1,8)))
reads.gsp3 <- reads %>% filter(sample %in% paste0("S", seq(9,16)))

rm(reads)

reads.gsp5 <- reads.gsp5 %>% filter(!is.na(dsodn))
reads.gsp3 <- reads.gsp3 %>% filter(!is.na(dsodn.rc))


# Extracting the sequences after the dsODN --------------------------------

reads.gsp5 <- reads.gsp5 %>%
  mutate(seq.to.blast = str_sub(seq, start = dsodn.end + 1))
reads.gsp3 <- reads.gsp3 %>%
  mutate(seq.to.blast = str_sub(seq, start = dsodn.rc.end + 1))

reads.gsp5 <- reads.gsp5 %>%
  mutate(seq.blast.length = str_length(seq.to.blast))
reads.gsp3 <- reads.gsp3 %>%
  mutate(seq.blast.length = str_length(seq.to.blast))

library(ggridges)

reads <- reads.gsp5 %>% bind_rows(reads.gsp3)

ggplot(reads)+
  geom_density_ridges(aes(x = seq.blast.length, y = sample, fill = sample),
                      alpha = 0.7)+
  scale_y_discrete(limits = rev(paste0("S", seq(19,36))))+
  theme_classic()+
  theme(legend.position = "none")

ggsave("____.svg", device = "svg", width = 4, height = 4)

# Collapsing the reads into clones ----------------------------------------
# The reads are collapsed into unique sequences to make it easier to process
# them further. The information about the reads number assigned to each clone is
# retained and used downstream for the quantification.
# Cropping the reads from dsODN end.
clones <- reads %>% group_by(seq.to.blast, sample) %>%
  summarise(read.count = n(),
            dsodn = nth(dsodn, 1),
            dsodn.rc = nth(dsodn.rc, 1))

clones$clone.id <- seq(nrow(clones))

clones <- clones %>% mutate(seq.length = str_length(seq.to.blast))

ggplot(clones)+
  geom_density_ridges(aes(x = seq.length, y = sample, fill = sample),
                      alpha = 0.7)+
  scale_y_discrete(limits = rev(paste0("S", seq(1,16))))+
  theme_classic()+
  theme(legend.position = "none")

ggsave("__LengthDistributionClones.svg", device = "svg", width = 4, height = 4)

# Filtering out short sequences
clones <- clones %>%
  filter(seq.length > 20)

rm(reads)

# Exporting the sequences as FASTA file
write.fasta(as.list(clones$seq.to.blast),
            as.list(clones$clone.id),
            "_____clonestoblast.fa")

# Run minimap2 with this FASTA file and translate the resulting SAM file to
# 6-column BED file


write.table(clones, "_______clones.txt", sep = "\t", dec = ",",
            row.names = FALSE, quote = FALSE)

write.table(reads, "______reads.txt", sep = "\t", dec = ",",
            row.names = FALSE, quote = FALSE)
rm(reads)





