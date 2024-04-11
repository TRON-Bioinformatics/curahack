# Prototype 2

Our hypothesis is that the change in consecutive gene annotation releases is rather low.
To measure how much the gene annotation changes between two releases, we propose to compare genes on the
hash of the corresponding exon coordinates across all of its transcripts.

## Baseline approach

Map genes based on the gene identifier has several limitations:
- Gene identifiers change
- Identical gene models may be referred by more than one gene identifier
- Different annotation providers use different identifiers

## Hypothesis

Comparing two or more gene annotations versions we will be able to identify those genes that have identical gene and transcripts models and hence:
1. quantify the difference between two or more gene annotations
2. identify a subset of genes that be jointly analysed using different gene annotations without the need for upstream reanalysis


## Proposed approach

Two genes are considered identical when the start and end coordinates of the gene and every exon of every transcript are identical.
This approach is identifier agnostic.

- Gene coordinates are matched by an exact match
- Exon coordinates are matched using an aggregated hash after sorting all of them


Given a gene, gene_1, with two transcripts and three exons each with start and end coordinates such as: 
- Transcript 1: 10-20, 30-40, 50-60
- Transcript 2: 10-20, 35-40, 50-60

```
hash(gene_1) = hash("10-20,10-20,30-40,35-40,50-60,50-60")
```

By basing the comparison on the hash of the exon coordinates, we can compare gene annotations across different releases
independently of gene identifiers and potentially between different sources of gene annotations, ie: Ensembl and NCBI.


## Usage

```
diff_gtf.py [-h] [--output OUTPUT] gtf_file gtf_file_2
```


## Limitations

- This approach only works on gene annotations based on the same reference genome.
- Current solution has only been tested with GTF from Ensembl and may require further adaptation for other sources of gene annotations.
- This approach could be extended to more than two releases by comparing the hash of the exon coordinates across all releases.
- Performance could be improved, current implementation takes around 30 minutes to process the differences between two human GTF files. Hash calculation could use a fasta implementation, internal data structures could be optimized, parallelization could be added.=
