process COUNT_KMERS {
    cpus 1
    memory "8g"
    tag "${name}"
    publishDir "${params.outdir}", mode: 'copy'

    conda ("${baseDir}/environments/biopython.yml")

    input:
      tuple val(name), path(fastq)

    output:
      tuple val("${name}"), path("kmer_counts.csv"), emit: kmer_counts

    script:
    """
    count_kmers.py \
        --input_fastq ${fastq} \
        --genomic_kmers ${params.genomic_kmers} \
        --output_csv kmer_counts.csv
    """
}

process ANNOTATE_KMERS {
    cpus 1
    memory "8g"
    tag "${name}"
    publishDir "${params.outdir}", mode: 'copy'

    input:
      tuple val(name), path(kmer_counts)

    output:
      tuple val("${name}"), path("results.csv"), emit: annotated_kmers

    script:
    """
    find_overlap.py \
        --input_counts ${kmer_counts} \
        --exon_coords ${params.exon_coords} \
        --output_csv kmer_counts.csv
    """
}