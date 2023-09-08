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


# Importing the clones alignment ------------------------------------------
# Specify the name of the BED file here
clones.aln <- read.table("______bed.txt", header = FALSE,
                         sep = " ", stringsAsFactors = FALSE)

clones <- read.table("_______clones.txt", header = TRUE,
                     sep = "\t", stringsAsFactors = FALSE)


colnames(clones.aln) <- c("chrom", "start", "end", "clone.id", "score",
                          "strand")

unique(clones.aln$chrom)

clones.aln <- clones.aln %>%
  filter(chrom %in% c("chrX", "chrY", "chrM", paste0("chr", seq(1,22))))

# Filtering out clones with ambiguous alignment
clones.aln.agg <- clones.aln.agg %>%
  left_join(clones %>% select(clone.id, read.count), by = "clone.id")

clones.aln <- clones.aln %>%
  filter(!clone.id %in% clones.aln.agg$clone.id[clones.aln.agg$aln.count > 1])

clones <- clones %>% left_join(clones.aln, by = "clone.id")

rm(clones.aln.agg, clones.aln)

clones <- clones %>%
  filter(!is.na(chrom))


write.table(clones, "_____clones_aligned.txt", sep = "\t", dec = ",",
            row.names = FALSE, quote = FALSE)



# Specify the sgRNA here
sgrna1 <- "CTTGCCCCACAGGGCAGTAA*GG"
sgrna2 <- "TCCACTCCTGATGCTGTTAT*GG"
target.size <- 23

# Specify the sample assignment here
ctrl.samples <- c()
treatment.samples <- c()
inhibitor.samples <- c()

# Specify the path to GRCh38 genes list from NCBI
hg38.genes <-
  read.table("hg38_genes_NCBI.txt", header = TRUE,
             sep = "\t", stringsAsFactors = FALSE)
hg38.genes <- hg38.genes %>%
  filter(!str_detect(chrom, "_") & str_detect(name, "NM")) %>%
  group_by(name2) %>%
  summarise(chrom = nth(chrom, 1), start = min(txStart),
            end = max(txEnd)) %>%
  mutate(gene.name = name2) %>% select(-name2)

# Symmetrical reads search ------------------------------------------------
clones <- clones %>%
  mutate(start.round = round(start/10)*10,
         end.round = round(end/10)*10)

clones.plus <- clones %>%
  filter(strand == "+") %>%
  group_by(chrom,strand,start.round) %>%
  summarise(samples.count.plus = n_distinct(sample),
            samples.plus = paste0(unique(sample), collapse = "|"),
            clone.count.plus = n(),
            read.count.plus = sum(read.count),
            seq.length.plus = mean(seq.length),
            starts.plus = paste0(unique(start), collapse = "|"),
            ends.plus = paste0(unique(end), collapse = "|")) %>%
  mutate(focal.point = start.round) %>%
  select(-start.round)

clones.minus <- clones %>%
  filter(strand == "-") %>%
  group_by(chrom,strand,end.round) %>%
  summarise(samples.count.minus = n_distinct(sample),
            samples.minus = paste0(unique(sample), collapse = "|"),
            clone.count.minus = n(),
            read.count.minus = sum(read.count),
            seq.length.minus = mean(seq.length),
            starts.minus = paste0(unique(start), collapse = "|"),
            ends.minus = paste0(unique(end), collapse = "|")) %>%
  mutate(focal.point = end.round) %>%
  select(-end.round)

clones.symmetry <- clones.plus %>%
  inner_join(clones.minus, by = c("chrom", "focal.point")) %>%
  mutate(clones.n.total = clone.count.plus + clone.count.minus,
         reads.n.total = read.count.plus + read.count.minus) %>%
  select(-strand.x,-strand.y) %>%
  arrange(desc(reads.n.total)) %>%
  ungroup()

clones.symmetry$common.pos <- NA
clones.symmetry$common.pos.n <- NA
clones.symmetry$samples.total <- NA
clones.symmetry$samples.n.total <- NA

for(i in seq(nrow(clones.symmetry))){
  # Overlapping the start/end coordinates
  starts <- clones.symmetry$starts.plus[i] %>% str_split("\\|") %>% unlist()
  ends <- clones.symmetry$ends.minus[i] %>% str_split("\\|") %>% unlist()
  common.pos.intersect <- intersect(starts,ends)
  clones.symmetry$common.pos[i] <- paste0(common.pos.intersect, collapse = "|")
  clones.symmetry$common.pos.n[i] <- length(common.pos.intersect)
  
  # Calculating in how many samples and in which ones these symmetrical reads are found
  samples.plus <- clones.symmetry$samples.plus[i] %>% str_split("\\|") %>% unlist()
  samples.minus <- clones.symmetry$samples.minus[i] %>% str_split("\\|") %>% unlist()
  clones.symmetry$samples.total[i] <- paste0(unique(c(samples.plus, samples.minus)), collapse = "|")
  clones.symmetry$samples.n.total[i] <- length(unique(c(samples.plus, samples.minus)))
}
rm(i, starts, ends, common.pos.intersect, samples.plus, samples.minus)

# Some info
# Common focal points found in...
print(paste0("Common focal points found for ",
             sum(clones.symmetry$common.pos != ""),
             " regions.",
             collapse = ""))

print(paste0("A single focal point found for ",
             sum(clones.symmetry$common.pos.n == 1),
             " regions.",
             collapse = ""))


clones.symmetry$target.seq <- NA
clones.symmetry$sgrna <- NA
clones.symmetry$min.ed <- NA
clones.symmetry$cut.pos <- NA

# Sequence extraction and ED measurement
for(i in seq(nrow(clones.symmetry))){
  target.chrom <- clones.symmetry$chrom[i]
  sgrna.alignments.df <- data.frame()
  if(clones.symmetry$common.pos.n[i] >= 1){
    common.positions <- clones.symmetry$common.pos[i] %>%
      str_split("\\|") %>% unlist() %>% as.integer()
    eds.tmp <- list()
    for(k in common.positions){
      target.start <- k
      
      # Checking the plus strand targeting
      target.seq <- as.character(Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg38,
                                                    target.chrom,
                                                    target.start-16,
                                                    target.start-16+target.size-1))
      sgrna.alignments.df <- sgrna.alignments.df %>%
        bind_rows(data.frame(target.seq = target.seq,
                             sgrna = "sgrna1",
                             strand = "plus",
                             cut.pos = target.start,
                             ed = StrDist(target.seq, sgrna1, method = "hamming")[1]))
      
      sgrna.alignments.df <- sgrna.alignments.df %>%
        bind_rows(data.frame(target.seq = target.seq,
                             sgrna = "sgrna2",
                             strand = "plus",
                             cut.pos = target.start,
                             ed = StrDist(target.seq, sgrna2, method = "hamming")[1]))
      

      
      # Checking the minus strand targeting
      target.seq <- as.character(Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg38,
                                                    target.chrom,
                                                    target.start-5,
                                                    target.start-5+target.size-1))
      target.seq <- target.seq %>% rc()
      
      sgrna.alignments.df <- sgrna.alignments.df %>%
        bind_rows(data.frame(target.seq = target.seq,
                             sgrna = "sgrna1",
                             strand = "minus",
                             cut.pos = target.start,
                             ed = StrDist(target.seq, sgrna1, method = "hamming")[1]))
      
      sgrna.alignments.df <- sgrna.alignments.df %>%
        bind_rows(data.frame(target.seq = target.seq,
                             sgrna = "sgrna2",
                             strand = "minus",
                             cut.pos = target.start,
                             ed = StrDist(target.seq, sgrna2, method = "hamming")[1]))

      
    }
    
    
  } else {
    # If no common start/end positions are found, use focal point and tiling around it
    target.start <- clones.symmetry$focal.point[i]-25
    window.size <- 50
    target.seq <- as.character(Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg38,
                                                  target.chrom,
                                                  target.start,
                                                  target.start+window.size-1))

    for(j in seq(str_length(target.seq)-(target.size-1))){
      target.seq.tmp <- str_sub(target.seq, start = j, end = j + target.size - 1)

      sgrna.alignments.df <- sgrna.alignments.df %>%
        bind_rows(data.frame(target.seq = target.seq.tmp,
                             sgrna = "sgrna1",
                             strand = "plus",
                             cut.pos = target.start+j+15,
                             ed = StrDist(target.seq.tmp, sgrna1, method = "hamming")[1]))
      sgrna.alignments.df <- sgrna.alignments.df %>%
        bind_rows(data.frame(target.seq = target.seq.tmp,
                             sgrna = "sgrna2",
                             strand = "plus",
                             cut.pos = target.start+j+15,
                             ed = StrDist(target.seq.tmp, sgrna2, method = "hamming")[1]))
      sgrna.alignments.df <- sgrna.alignments.df %>%
        bind_rows(data.frame(target.seq = target.seq.tmp %>% rc(),
                             sgrna = "sgrna1",
                             strand = "minus",
                             cut.pos = target.start+j+4,
                             ed = StrDist(target.seq.tmp %>% rc(), sgrna1, method = "hamming")[1]))
      
      sgrna.alignments.df <- sgrna.alignments.df %>%
        bind_rows(data.frame(target.seq = target.seq.tmp %>% rc(),
                             sgrna = "sgrna2",
                             strand = "minus",
                             cut.pos = target.start+j+4,
                             ed = StrDist(target.seq.tmp %>% rc(), sgrna2, method = "hamming")[1]))

    }
    
  }

  
  sgrna.alignments.df <- sgrna.alignments.df %>%
    filter(ed == min(ed))
  clones.symmetry$target.seq[i] <- paste0(sgrna.alignments.df$target.seq, collapse = "|")
  clones.symmetry$sgrna[i] <- paste0(sgrna.alignments.df$sgrna, collapse = "|")
  clones.symmetry$cut.pos[i] <- paste0(sgrna.alignments.df$cut.pos, collapse = "|")
  clones.symmetry$min.ed[i] <- sgrna.alignments.df$ed[1]

  rm(target.seq.tmp, sgrna.tmp, j, k, min.ed, which.sgrna, eds.tmp,
     target.seq, ed.plus.sgrna1, ed.plus.sgrna2, ed.minus.sgrna1,
     ed.minus.sgrna2, target.start, common.positions)
  
  
  
  if(i %in% seq(1, nrow(clones.symmetry), 100)){
    print(paste0(round(i*100/nrow(clones.symmetry),1), "% done"))
  }
  
}
rm(i, target.chrom, target.size, window.size, sgrna.alignments.df)
beep()
# Some info
# Common focal points found in...
print(paste0("Less than 11 nt ED found for ",
             sum(clones.symmetry$min.ed < 11),
             " regions.",
             collapse = ""))

print(paste0("Number of common-focal-point regions with less than 11 nt ED: ",
             sum(clones.symmetry$min.ed < 11 & clones.symmetry$common.pos.n >= 1),
             collapse = ""))

print(paste0("Number of single-focal-point regions with less than 11 nt ED: ",
             sum(clones.symmetry$min.ed < 11 & clones.symmetry$common.pos.n == 1),
             collapse = ""))


write.table(clones.symmetry, "___SymmtericalReads.txt",
            quote = FALSE, sep = "\t", dec = ".",
            row.names = FALSE)

clones.agg.samples.plus <- clones %>%
  filter(strand == "+") %>%
  group_by(chrom,strand,start.round,sample) %>%
  summarise(clone.count.plus = n(),
            read.count.plus = sum(read.count)) %>%
  mutate(focal.point = start.round) %>%
  ungroup() %>%
  select(-start.round)

clones.agg.samples.minus <- clones %>%
  filter(strand == "-") %>%
  group_by(chrom,strand,end.round,sample) %>%
  summarise(clone.count.minus = n(),
            read.count.minus = sum(read.count)) %>%
  mutate(focal.point = end.round) %>%
  ungroup() %>%
  select(-end.round)

clones.agg.samples <- clones.agg.samples.plus %>%
  full_join(clones.agg.samples.minus, by = c("sample", "chrom", "focal.point")) %>%
  mutate(clone.count.plus = replace_na(clone.count.plus, 0),
         clone.count.minus = replace_na(clone.count.minus, 0),
         read.count.plus = replace_na(read.count.plus, 0),
         read.count.minus = replace_na(read.count.minus, 0),
         clones.n = clone.count.plus + clone.count.minus,
         reads.n = read.count.plus + read.count.minus) %>%
  select(chrom,focal.point,sample,reads.n)

clones.agg.samples.ctrl <- clones.agg.samples %>%
  mutate(sample.type = ifelse(sample %in% inhibitor.samples, "trt.w.inh",
                              ifelse(sample %in% ctrl.samples, "ctrl", "trt.wo.inh"))) %>%
  group_by(chrom, focal.point, sample.type) %>%
  summarise(reads.n = sum(reads.n, na.rm = TRUE)) %>%
  pivot_wider(names_from = "sample.type", values_from = "reads.n")



clones.symmetry.simple <- clones.symmetry %>%
  select(chrom,focal.point,cut.pos,min.ed,target.seq,sgrna,
         samples.n.total,samples.total,clones.n.total,reads.n.total) %>%
  left_join(clones.agg.samples.ctrl, by = c("chrom","focal.point")) %>%
  mutate_at(c("trt.w.inh", "trt.wo.inh", "ctrl"), ~replace_na(.,0))


clones.symmetry.simple <- clones.symmetry.simple %>%
  arrange(min.ed)


clones.symmetry.simple$gene <- NA

for(i in seq(nrow(clones.symmetry.simple))){
  hg38.genes.this.chrom <- hg38.genes %>% filter(chrom == clones.symmetry.simple$chrom[i])
  cut <- clones.symmetry.simple$cut.pos[i] %>% str_split("\\|") %>% unlist() %>% as.integer()
  genes <- c()
  for(k in cut){
    genes <- c(genes, hg38.genes.this.chrom$gene.name[hg38.genes.this.chrom$start <= k &
                                              hg38.genes.this.chrom$end >= k])

  }
  
  if(length(genes) > 0){
    clones.symmetry.simple$gene[i] <- paste0(genes, collapse = "|")
  }
  rm(genes, hg38.genes.this.chrom)
}


write.table(clones.symmetry.simple, "__SymmetryReadsSimplified.txt",
            quote = FALSE, sep = "\t", dec = ".",
            row.names = FALSE)

