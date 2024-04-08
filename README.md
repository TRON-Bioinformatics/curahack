# curahack - anotation free expression analysis from RNAseq challenge
This repository contains test data for the curatime curahack hackathon. for further information on the curatime cluster and related events, please see here:
https://curatime.org/


# Description of the challenge
Gene expression analysis from bulk RNAseq data is a very well established routine process: Starting from raw sequencing data (fastq), reads are either aligned to the human genome (e.g. using star) and based on their location subsequently assigned to a given gene model (e.g. ensemble) or directly assigned using an alignment free process (e.g. kallisto). All these processes have in common that a specific gene model is required and in case a different gene model was used, the data is not directly comparable for downstream analysis. Thus, going back to the original raw sequencing data is often required to ensure comparability.


# Pain point
While the human genome reference stays the same over long periods of time, gene annotations are frequently updated and different databases exist (e.g. ucsc, ensemble, refseq, …). Ensemble, which is frequently used by many researchers, exists already in version 112 (March 2024), with several updates troughout the year; thus, frequent reanalysis of the original data is required. This is the case for in-house generated data but also for public available data, which often requires downloading of large raw sequencing files and storage of those files for long periods to enable updates on the analysis and consistency. 

For example the GDC data portal contains 44TB of raw bulk RNAseq data, that is also available preprocessed as gene and transcript level expression data available. However, when comparing to other datasets in order to assure same annotation, reprocessing of such large ammounts of data can often not be avoided, despite the presence of preprocessed expression values.

# Desired state
An annotation free preprocessing step that maintains quantitative information without the need to use a specific annotation as an intermediate step would be highly desirable. The intermediate data should allow downstream gene level expression analysis and produce count and TPM values that are comparable to standard analysig. At the same time file size should be drastically reduced compared to (gzipped) bam or fastq files.

Imaginable would be a compact file that stores information on read mappings after alignment to a genome, e.g. wig file formats have been proposed in that regard before. Also, an alignment free, quantitative kmer based approach would be imaginable and potentially other approaches exist that could be explored during the hackathon.

# Benchmark 
We provide gene level expression tables for different version of ensemble annotation for a dataset of 9 samples for 3 conditions (each condition with 3 biological replica). These serve as gold standard gene level expression quantification for each annotation. In order to allow easier handling of files reads were filtered ony for reads mapping to chromosome 1, hence only genes from chromosome 1 can be quantified using the provided test dataset. The dataset can also be used to identify differential expressed genes for the 3 provided conditions. Any preprocessing solution that you develope should allow to assess  gene level expression and maintain differences between the 3 biological conditions. At the same time the pre-processed data file should be as small as possible.

Results will be judged according to how well gene level quantification results after preprocessing with your proposed tool (set of tools) correlate with the provided gene level expression data and whether differences between the 3 conditions are maintained. Another crucial benchmarking parameter is the filesize of the preprocessed file that should be as small as possible.

# Existing approaches
Please find here a possible approach that you may find usefull as a starting point for the proposed challenge. Feel free to utilize some or all tools presented here. You could try to optimize settings for some of the tools outlined here, replace parts or come up with a completely different approach.

Previously the use of kmers was proposed for annotation free quantifcation on gene and transcript level, as outlined in the folowing publication:
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8221386/
Please be aware that we are not affilated with the authors, it should just serve as suggestion. Within the publication the authors propose to extract transcript/gene specific kmer signatures (Kmerator) which can later be quantified from RNAseq fastq files (counttags). The generation of such kmer signatures can be done using the provided tool Kmerator (https://github.com/Transipedia/kmerator). The output of kmerator is a set of fasta files that contain relvant kmer sequences. It would be consivable to use the tool and optimize input paramezers (e.g. kmer-length) for reduced file size and gene expression level accuracy. You could also try to reduce size, by limiting the number of kmers per gene. You could also write a completly novel allgorythm for the definition of most informative kmers per gene that are outputed in a fasta file.

After generation of fasta files that contain information on the kmer sequences, kmers are subsequently quantified from fastq files using the tool counttags (https://github.com/Transipedia/countTags), which provides kmer abundance as counttable. The count table would represent preprocessed intermeditate files that should be investigate wehther gene level expression can accurately determined.

# Related links

- star repository: https://github.com/alexdobin/STAR
- kallisto repository: https://github.com/pachterlab/kallisto
- ensembl: http://www.ensembl.org/info/genome/genebuild/index.html
- kmerator: https://github.com/Transipedia/kmerator
- counttags: https://github.com/Transipedia/countTags
- original raw data used for generation of test data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE217893
- gdc data protal: https://portal.gdc.cancer.gov
