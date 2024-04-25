#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GENOMIC_KMERS ; ANNOTATE_EXONS } from './modules/build_refs.nf'
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

// checks required inputs
if (params.mode == "build_references") {
    if (params.ref_genome && params.ref_gtf) {
        log.info "Building references"
    } else {
        exit 1, "Reference genome and GTF file not specified!"
    }
} else if (params.mode == "quantify_kmers") {
    if (params.input_table) {
        Channel
            .fromPath(params.input_table)
            .splitCsv(header: ['name', 'fastq'], sep: ";")
            .map{ row-> tuple(row.name, row.fastq) }
            .set { input_files }
    } else {
        exit 1, "Input file not specified!"
    }
}

workflow BUILD_REFERENCES {

    take:
    genome_fasta
    gtf_file

    main:
    GENOMIC_KMERS(genome_fasta)
    ANNOTATE_EXONS(gtf_file)

    emit:
    genomic_kmers = GENOMIC_KMERS.out.genomic_kmers
    exon_coord = ANNOTATE_EXONS.out.exon_coord

}


workflow QUANTIFY_KMERS {

    take:
    input_files

    main:
    COUNT_KMERS(input_files)
    ANNOTATE_KMERS(COUNT_KMERS.out.kmer_counts)

    emit:
    ANNOTATE_KMERS.out.annotated_kmers

}

workflow {
    if (params.mode == "build_references"){
        BUILD_REFERENCES(
            Channel.fromPath(params.ref_genome).combine(
            Channel.fromPath(params.ref_gtf)
        ))
    } else if (params.mode == "quantify_kmers"){
        QUANTIFY_KMERS(input_files)
    } else {
        log.error "Invalid mode"
        exit 1
    }
}

// create nextflow pipeline with fastq as input parameter
// nextflow run main.nf -profile conda --input_table ref_data/input_table.csv --outdir <output_dir>