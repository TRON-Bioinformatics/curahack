#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GENOMIC_KMERS ; ANNOTATE_EXONS } from './modules/build_refs.nf'

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
if (!params.ref_genome) {
    log.error "--ref_genome is required"
    exit 1
}

if (!params.ref_gtf) {
    log.error "--ref_gtf is required"
    exit 1
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

workflow {
    BUILD_REFERENCES(
        file(params.ref_genome),
        file(params.ref_gtf)
    )
}

// create nextflow pipeline with fastq as input parameter
// nextflow run main.nf -profile conda --input_table ref_data/input_table.csv --outdir <output_dir>