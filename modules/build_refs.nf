process GENOMIC_KMERS {
    cpus 1
    memory "40g"
    tag "GENOME_FASTA"
    publishDir "${params.outdir}", mode: 'copy'

    conda ("${baseDir}/environments/biopython.yml")

    input:
      path(fastq)

    output:
      path("genomic_kmers.fa"), emit: genomic_kmers

    script:
    """
    generate_genomic_kmers.py \
        --input_fasta ${params.ref_genome} \
        --step_size ${params.step_size} \
        --kmer_size ${params.kmer_size} \
        --output_fasta genomic_kmers.fa
    """
}

process ANNOTATE_EXONS {
    cpus 1
    memory "2g"
    tag "GTF"
    publishDir "${params.outdir}", mode: 'copy'

    input:
      path(fastq)

    output:
      path("exon_coordinates.csv"), emit: exon_coord

    script:
    """
    gtf_to_coord.py \
        --input_gtf ${params.ref_gtf} \
        --output_csv exon_coordinates.csv
    """
}