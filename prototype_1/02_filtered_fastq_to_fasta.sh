\##################################

\####### Filtered FASTQ to FASTA ########

\##################################

dsk -file Control_siRNA_1.trimmomatic.fastq -kmer-size 72

dsk2ascii -file Control_siRNA_1.trimmomatic.h5 -out output.txt

python convert.py output.txt -> converted_output.fasta