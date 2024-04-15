bowtie2-build -x reference.fa reference_index
bowtie2 -x reference_index.fa -f converted_dsk_output.fa -S bowtie_output.sam
