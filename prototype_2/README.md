# Prototype 2

Our hypothesis is that the change in consecutive gene annotation releases is rather low.
To measure how much the gene annotation changes between two consecutive releases, we propose to compare genes on the
hash of the corresponding exon coordinates across all of its transcripts.

Given a gene, gene_1, with two transcripts and three exons each with start and end coordinates such as: 
- Transcript 1: 10-20, 30-40, 50-60
- Transcript 2: 10-20, 35-40, 50-60

```
hash(gene_1) = hash("10-20,10-20,30-40,35-40,50-60,50-60")
```

By basing the comparison on the hash of the exon coordinates, we can compare gene annotations across different releases
independently of gene identifiers and potentially between different sources of gene annotations, ie: Ensembl and NCBI.

By this means we obtain two products:
1. The list of genes that have not changed between two releases. Their respective gene identifiers in each release
2. A measure of how much the gene annotation has changed between two consecutive releases based on the Jaccard index.


## Usage

```
diff_gtf.py [-h] [--output OUTPUT] gtf_file gtf_file_2
```


## Limitations

- This approach only works on gene annotations based on the same reference genome.
- Current solution has only been tested with GTF from Ensembl and may require further adaptation for other sources of gene annotations.
- This approach could be extended to more than two releases by comparing the hash of the exon coordinates across all releases.
- Performance could be improved, current implementation takes around 30 minutes to process the differences between two human GTF files. Hash calculation could use a fasta implementation, internal data structures could be optimized, parallelization could be added.=
