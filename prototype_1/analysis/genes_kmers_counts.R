### combine kmers, genes and counts ###
sam <- read.table("/Users/hoehnemi/Downloads/curhack/bowtie_output_nh.sam", fill = T, sep = "\t")
kmer_counts <- read.table("/Users/hoehnemi/Downloads/curhack/output.txt")
colnames(kmer_counts) <- c("kmer","counts")

kmer_counts$name <- 1:nrow(kmer_counts)  
kmer_counts$name <- gsub("^", "kmer",kmer_counts$name)

# extract algined kmers
sam_algined <- sam[sam$V3 !="*" ,c(1,3,10)]
sam_algined <- sam_algined[seq(1, nrow(sam_algined), 2), ]
#t <- as.data.frame(table(sam_algined$V6))
# rename to chr regions as one name
sam_algined$V3 <- gsub(":", "_",sam_algined$V3)
sam_algined$V3<- gsub("-", "_",sam_algined$V3)
colnames(sam_algined) <- c("name","chr", "kmer")


# load gtf bed file
genes_bed <- read.table("/Users/hoehnemi/Downloads/curhack/genes_hg38_109_both.bed")
colnames(genes_bed) <- c("chr_name","start","end","gene_symbol", "strand")
genes_bed$chr <- paste(genes_bed[,1],genes_bed[,2],genes_bed[,3],sep='_')

# merge kmer regions and counts 
kmer_combi <- merge(kmer_counts, sam_algined, by = "name")

# kmer combi and genes
genes_small <- genes_bed[, c(4,6)]
kmer_genes<- merge(kmer_combi, genes_small, by = "chr")

library(dplyr)
count.matrix <- kmer_genes %>%
  group_by(gene_symbol) %>%
  summarise(count=sum(counts))


kallisto_si1 <- read.csv("/Users/hoehnemi/Downloads/curhack/gene_expr_ctr_siRNA_1_109.csv", sep = ";")
gtf_bed <- read.table("/Users/hoehnemi/Downloads/curhack/genes_hg38_109_both.bed", sep = "\t")

# only chr 1
gtf_chr1 <- gtf_bed[gtf_bed$V1 =="chr1",]

kallisto_si1_chr1 <- kallisto_si1[kallisto_si1$gene_symbol %in% gtf_chr1$V4,]
# combine 
combi <- merge(count.matrix, kallisto_si1_chr1, by = "gene_symbol",all=T)
combi[is.na(combi)] <- 0
colnames(combi) <- c("gene_symbol","counts_kmers","gene_id","avg_len","counts_kallisto","rpkm","tpm")
plot(combi$counts_kmers, combi$counts_kallisto, ylim = c(0,5000), xlim = c(0,5000))


library("ggpubr")
ggscatter(combi, x = "counts_kmers", y = "counts_kallisto", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman")+
  ylim(0,5000)+
  xlim(0,5000)
ggsave("/Users/hoehnemi/Downloads/curhack/correlation_kmer_kallisto.pdf")

