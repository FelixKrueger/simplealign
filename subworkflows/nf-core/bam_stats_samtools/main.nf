//
// Run SAMtools stats, flagstat and idxstats
//

include { SAMTOOLS_STATS    } from '../../../modules/nf-core/samtools/stats/main'
include { SAMTOOLS_IDXSTATS } from '../../../modules/nf-core/samtools/idxstats/main'
include { SAMTOOLS_FLAGSTAT } from '../../../modules/nf-core/samtools/flagstat/main'

workflow BAM_STATS_SAMTOOLS {
    
    take:
    ch_bam_bai    // channel: [ val(meta), path(bam), path(bai) ]
    // Felix: I believe this is not the correct input channel.ch_fasta   // channel: [ val(meta), path(fasta) ]
    // This is:
    // ch_fasta   //channel: /path/to/reference.fasta

    main:
    ch_versions = Channel.empty()

    // ch_bam_bai.view()
    // println("ch_bam_bai is:")
    // ch_bam_bai.view() //works, tested it

    SAMTOOLS_STATS ( ch_bam_bai )
    // ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions.first()))

    SAMTOOLS_FLAGSTAT ( ch_bam_bai )
    ch_versions = ch_versions.mix(SAMTOOLS_FLAGSTAT.out.versions.first())

    SAMTOOLS_IDXSTATS ( ch_bam_bai )
    ch_versions = ch_versions.mix(SAMTOOLS_IDXSTATS.out.versions.first())

    emit:
    // stats    = SAMTOOLS_STATS.out.stats       // channel: [ val(meta), path(stats) ]
    flagstat = SAMTOOLS_FLAGSTAT.out.flagstat // channel: [ val(meta), path(flagstat) ]
    idxstats = SAMTOOLS_IDXSTATS.out.idxstats // channel: [ val(meta), path(idxstats) ]

    versions = ch_versions                    // channel: [ path(versions.yml) ]
}
