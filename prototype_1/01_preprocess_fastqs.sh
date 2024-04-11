# runs QC on input FASTQs
fastqc downsampled_fastq/*.fastq.gz

# runs Trimmomatic on input FASTQs cropping to 72 bp, min length 72 bp and sliding window of 2 bp below 32 Phred score
trimmomatic SE \
	downsampled_fastq/Control_siRNA_1.fastq.gz \
	trimmed_fastq/Control_siRNA_1.trimmomatic.fastq.gz \
	CROP:72 MINLEN:72 SLIDINGWINDOW:32:2

trimmomatic SE \
	downsampled_fastq/Control_siRNA_2.fastq.gz \
	trimmed_fastq/Control_siRNA_2.trimmomatic.fastq.gz \
	CROP:72 MINLEN:72 SLIDINGWINDOW:32:2

trimmomatic SE \
	downsampled_fastq/Control_siRNA_3.fastq.gz \
	trimmed_fastq/Control_siRNA_3.trimmomatic.fastq.gz \
	CROP:72 MINLEN:72 SLIDINGWINDOW:32:2

trimmomatic SE \
	downsampled_fastq/STAT5A_siRNA_1.fastq.gz \
	trimmed_fastq/STAT5A_siRNA_1.trimmomatic.fastq.gz \
	CROP:72 MINLEN:72 SLIDINGWINDOW:32:2

trimmomatic SE \
	downsampled_fastq/STAT5A_siRNA_2.fastq.gz \
	trimmed_fastq/STAT5A_siRNA_2.trimmomatic.fastq.gz \
	CROP:72 MINLEN:72 SLIDINGWINDOW:32:2

trimmomatic SE \
	downsampled_fastq/STAT5A_siRNA_3.fastq.gz \
	trimmed_fastq/STAT5A_siRNA_3.trimmomatic.fastq.gz \
	CROP:72 MINLEN:72 SLIDINGWINDOW:32:2

trimmomatic SE \
	downsampled_fastq/STAT5B_siRNA_1.fastq.gz \
	trimmed_fastq/STAT5B_siRNA_1.trimmomatic.fastq.gz \
	CROP:72 MINLEN:72 SLIDINGWINDOW:32:2

trimmomatic SE \
	downsampled_fastq/STAT5B_siRNA_2.fastq.gz \
	trimmed_fastq/STAT5B_siRNA_2.trimmomatic.fastq.gz \
	CROP:72 MINLEN:72 SLIDINGWINDOW:32:2

trimmomatic SE \
	downsampled_fastq/STAT5B_siRNA_3.fastq.gz \
	trimmed_fastq/STAT5B_siRNA_3.trimmomatic.fastq.gz \
	CROP:72 MINLEN:72 SLIDINGWINDOW:32:2

# runs QC again after trimming and filters
fastqc trimmed_fastq/*.fastq.gz
