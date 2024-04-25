#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { COUNT_KMERS ; ANNOTATE_KMERS } from './modules/quantify_kmers.nf'

def helpMessage() {
    log.info params.help_message
}

if (params.help) {
    helpMessage()
    exit 0
}

if (!params.outdir) {
    log.error "--outdir is required"
    exit 1
}

if (!params.genomic_kmers) {
    log.error "--genomic_kmers is required"
    exit 1
}

if (!params.exon_coords) {
    log.error "--exon_coords is required"
    exit 1
}

// checks required inputs
if (params.input_table) {
    Channel
        .fromPath(params.input_table)
        .splitCsv(header: ['name', 'fastq'], sep: ";")
        .map{ row-> tuple(row.name, row.fastq) }
        .set { input_files }
} else {
    exit 1, "Input file not specified!"
}

workflow QUANTIFY_KMERS {

    take:
    input_files

    main:
    COUNT_KMERS(
        input_files
    )
    ANNOTATE_KMERS(
        COUNT_KMERS.out.kmer_counts
    )

    emit:
    ANNOTATE_KMERS.out.annotated_kmers

}

workflow {
    QUANTIFY_KMERS(input_files)
}

// create nextflow pipeline with fastq as input parameter
// nextflow run main.nf -profile conda --input_table ref_data/input_table.csv --outdir <output_dir>