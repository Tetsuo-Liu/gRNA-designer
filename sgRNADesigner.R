#This program was created to design gRNAs that introduce mutations specifically at the 5' splice donor site.
#This program is intended for use with SpCas9(PAM:NGG)
#Enter your target gene name on line 11 and it will be executed.
#If your gene of interest has several mRNA variants, you can also select a specific variant in line 17.

library(biomaRt)
library(dplyr)
library(stringr)
library(Biostrings)

target.gene <- "FOXP3" #Substitute the name of your target gene with mgi_symbol.

#Convert target.gene from mgi_symbol to refseq_mrna
mus <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
target.gene.refseq <- getBM(attributes = "refseq_mrna", filters = "mgi_symbol", values = target.gene, mart = mus)
target.gene.refseq
target.gene.refseq <- target.gene.refseq[1,] #If target.gene has several mRNA variant, target.gene.refseq will contain more than one refseq_mrna. Select one desirable refseq_mrna.   

#Determine whether the target gene is a plus-stranded or minus-stranded 
strand <-  getBM(attributes = "strand", filters = "refseq_mrna", values = target.gene.refseq, mart = mus)
strand <- strand[,1]

#Obtain target.gene.sequence
target.gene.sequence <- getSequence(id = target.gene.refseq, type = "refseq_mrna", seqType = "gene_exon_intron", mart = mus) 
target.gene.sequence <- target.gene.sequence[,1]

#Obtain the sequence of each exon
exon.sequence <- getBM(attributes = c("exon_chrom_end", "gene_exon"), filters = "refseq_mrna", values = target.gene.refseq, mart = mus)
exon <- c()
for (i in 1:length(exon.sequence[,1])) {
  exon <- c(exon, i)
}
exon <- as.data.frame(exon)
if (strand == 1) {
  exon.sequence <- exon.sequence %>%
    arrange(exon_chrom_end) %>%
    select(-exon_chrom_end) %>%
    bind_cols(exon) %>%
    select(exon, gene_exon)
} else if (strand == -1) {
  exon.sequence <- exon.sequence %>%
    arrange(desc(exon_chrom_end)) %>%
    select(-exon_chrom_end) %>%
    bind_cols(exon) %>%
    select(exon, gene_exon)
}

#Search appropriate PAM sequences
#Sense strand
position.sense <- data.frame(matrix(rep(NA, 2), nrow=1))[numeric(0), ]
colnames(position.sense) <- c("start", "end") 
for (i in 1:(length(exon.sequence[,1])-1)) { 　
  a <- as.data.frame(str_locate(target.gene.sequence, exon.sequence[i,2]))
  position.sense <- full_join(position.sense, a)
}
position.sense <- position.sense %>%
  mutate(PAM_start = end + 2) %>%
  mutate(PAM_end = end +11) %>%
  mutate(PAM_range = str_sub(target.gene.sequence, start = PAM_start, end = PAM_end)) %>%
  mutate("PAM?" = str_detect(PAM_range, "GG")) 
#Antisense strand
position.antisense <- position.sense %>%
  mutate(PAM_start = end -7) %>%
  mutate(PAM_end = end +2) %>%
  mutate(PAM_range = str_sub(target.gene.sequence, start = PAM_start, end = PAM_end)) %>%
  mutate("PAM?" = str_detect(PAM_range, "CC"))

#Calculate appropriate sgRNA sequence
#sense strand
PAM.sense <- data.frame(matrix(rep(NA, 3), nrow=1))[numeric(0), ] 
colnames(PAM.sense) <- c("exon", "start", "end") 
for (i in 1:length(position.sense[,6])) {
  if (position.sense[i,6] == TRUE) {
    b <- as.data.frame(str_locate_all(position.sense[i,5], "GG")) %>%
      mutate(exon = i)
    PAM.sense <- full_join(PAM.sense, b)
  }
}
PAM.sense <- PAM.sense %>%
  mutate(PAM_position = position.sense[exon,3] + start -2) %>%
  mutate(PAM = str_sub(target.gene.sequence, start = PAM_position, end = PAM_position +2)) %>%
  mutate(sgRNA = str_sub(target.gene.sequence, start = PAM_position -20, end = PAM_position -1)) %>%
  mutate(length = 20) %>%
  mutate("%GC" = str_count(sgRNA, "G|C")/length *100) %>%
  mutate(Orientation = "Sense") %>%
  select(c(-start, -end, -PAM_position))
#Antisense strand
PAM.antisense <- data.frame(matrix(rep(NA, 3), nrow=1))[numeric(0), ] 
colnames(PAM.antisense) <- c("exon", "start", "end")
for (i in 1:length(position.antisense[,6])) {
  if (position.antisense[i,6] == TRUE) {
    c <- as.data.frame(str_locate_all(position.antisense[i,5], "CC")) %>%
      mutate(exon = i)
    PAM.antisense <- full_join(PAM.antisense, c)
  }
}
PAM.antisense <- PAM.antisense %>%
  mutate(PAM_position = position.antisense[exon,3] + start -1) %>%
  mutate(PAM = str_sub(target.gene.sequence, start = PAM_position, end = PAM_position +2)) %>%
  mutate(sgRNA = str_sub(target.gene.sequence, start = PAM_position +3, end = PAM_position +22)) %>%
  mutate(length = 20) %>%
  mutate("%GC" = str_count(sgRNA, "G|C")/length *100) %>%
  mutate(Orientation = "Antisense") %>%
  select(c(-start, -end, -PAM_position))
if (any(position.antisense[,6]) != FALSE) {
  rev_PAM <- c()
  rev_sgRNA <- c()
  for (i in 1:length(PAM.antisense[,1])) {
    d <- as.character(reverseComplement(DNAString(PAM.antisense[i,2])))
    e <- as.character(reverseComplement(DNAString(PAM.antisense[i,3])))
    rev_PAM <- c(rev_PAM, d)
    rev_sgRNA <- c(rev_sgRNA, e)
  }
  rev <- tibble(rev_PAM, rev_sgRNA)
  colnames(rev) <- c("PAM", "sgRNA")
  PAM.antisense <- PAM.antisense %>%
    select(-PAM,-sgRNA) %>% 
    bind_cols(rev) %>%
    select(exon, PAM, sgRNA, everything())
}

#Combine PAM.sense and PAM.antisense, generate Result
Result <- full_join(PAM.sense, PAM.antisense)
Result
