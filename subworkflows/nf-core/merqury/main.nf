//
// Kmer-based assembly evaluation
//

include { MERYL_COUNT     } from '../../../modules/nf-core/meryl/count/main'
include { MERYL_UNIONSUM  } from '../../../modules/nf-core/meryl/unionsum/main'
include { MERQURY_MERQURY } from '../../../modules/nf-core/merqury/merqury/main'

workflow MERQURY {

    take:
    fastq  // channel: [ val(meta), [ fastq ] ]
    kvalue // integer: >0

    main:

    ch_versions = Channel.empty()

    // Make sure FASTQ items are not empty, and named appropriately.
    // It's possible that the genomeqc workflow could pass tuples some of which have empty fastqs.
    fastq
        | map{meta, fq -> fq ? [meta, file(fq)] : [meta, fq]}
        | filter { meta, fq -> fq && fq.name =~ /(\.fastq|\.fq|\.fastq\.gz|\.fq\.gz)$/ }
        | set {fastq}

    // MODULE: MERYL_COUNT
    MERYL_COUNT(
        fastq,
        kvalue
    )
    ch_meryl_db = MERYL_COUNT.out.meryl_db
    ch_versions = ch_versions.mix(MERYL_COUNT.out.versions.first())

    // MODULE: MERYL_UNIONSUM
    MERYL_UNIONSUM(
        ch_meryl_db,
        kvalue
    )
    ch_meryl_union = MERYL_UNIONSUM.out.meryl_db
    ch_versions = ch_versions.mix(MERYL_UNIONSUM.out.versions.first())

    // MODULE: MERQURY_MERQURY
    ch_meryl_union
        | join(fastq)
        | set {ch_merqury_inputs}
    MERQURY_MERQURY ( ch_merqury_inputs )
    ch_merqury_qv                           = MERQURY_MERQURY.out.assembly_qv
    ch_merqury_stats                        = MERQURY_MERQURY.out.stats
    ch_merqury_spectra_cn_fl_png            = MERQURY_MERQURY.out.spectra_cn_fl_png
    ch_merqury_spectra_asm_fl_png           = MERQURY_MERQURY.out.spectra_asm_fl_png
    ch_hapmers_blob_png                     = MERQURY_MERQURY.out.hapmers_blob_png
    ch_merqury_outputs                      = ch_merqury_qv
                                            | mix(ch_merqury_stats)
                                            | mix(ch_merqury_spectra_cn_fl_png)
                                            | mix(ch_merqury_spectra_asm_fl_png)
                                            | mix(ch_hapmers_blob_png)
                                            | flatMap { meta, data -> data }
    ch_versions                             = ch_versions.mix(MERQURY_MERQURY.out.versions.first())


    emit:
    // TODO nf-core: edit emitted channels
    merqury  = ch_merqury_outputs  // channel: TODO [ val(meta), [ bam ] ]
    versions = ch_versions         // channel: [ versions.yml ]
}

