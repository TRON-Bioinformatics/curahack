#!/bin/bash

# This pipeline will use a kmer based approach to store RNA seq read data in a very compact format
# which can then be applied to different annotations


# `Homo_sapiens.GRCh38.dna.chromosome.1.fa` and `Homo_sapiens.GRCh38.109.gtf` can be obtained from ensembl FTP:
# wget https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.1.fa.gz
# wget https://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz

# Generate genomic kmers and save them into a FASTA file -> ids include genomic positions
python src/generate_genomic_kmers.py -i ref_data/Homo_sapiens.GRCh38.dna.chromosome.1.fa -s 15 -k 31 -o genomic_kmers_chr1.fasta

# Generate kmers from FASTQ file and map them against genomic kmers and create a count table
# This will be your intermediate file to store long-term
python src/generate_kmer_count_table_from_fastq.py -i ref_data/Control_siRNA_1.fastq -g results/genomic_kmers_chr1.fasta -k31 -o results/kmer_counts_chr1_k31_Control_siRNA_1.csv

# Generate exonic regions table from GTF file to overlap with the genomic coordinates
python src/gtf_to_coord.py -i ref_data/Homo_sapiens.GRCh38.109.gtf -o ref_data/exons_chr1_v109.csv

# Overlap kmer counts with exonic regions to create gene counts
python src/find_overlap.py -i results/kmer_counts_chr1_k31_Control_siRNA_1.csv -c ref_data/exons_chr1_v109.csv -o results/results.csv