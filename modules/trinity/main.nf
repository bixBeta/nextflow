runmode = params.mode

process TRINITY {

    maxForks 1
    label 'process_high'
    tag "$id"

    publishDir "trinity_assembly" , overwrite: true

    input:
        tuple val(id), path(r1_reads, stageAs: '??_r1.fq.gz'), path(r2_reads, stageAs: '??_r2.fq.gz')
        val(libtype)

    output:
        path "${id}_trinityAssemblyOutput/Trinity.fasta", emit: fasta
        path "*"

    script:

    def left    = r1_reads instanceof List ? r1_reads.join(',') : r1_reads
    def right   = r2_reads instanceof List ? r2_reads.join(',') : r2_reads
    def ss_flag = libtype ? "--SS_lib_type ${libtype}" : ""

    if (runmode == "PE" || runmode == "PES" || runmode == "PEBS") {

    """
        Trinity --seqType fq \\
            --max_memory ${task.memory.toGiga()}G \\
            --CPU ${task.cpus} \\
            ${ss_flag} \\
            --min_contig_length 400 \\
            --NO_SEQTK --output ${id}_trinityAssemblyOutput \\
            --left ${left} \\
            --right ${right} >& trinity.log
    """

    } else {
        error "Runmode ${runmode} is not supported"
    }

}


// ─────────────────────────────────────────────────────────────────────────────
// TRINITY_STATS — basic assembly QC (N50, total bases, contig count)
// ─────────────────────────────────────────────────────────────────────────────
process TRINITY_STATS {
    label 'process_low'
    publishDir "trinity_assembly"   // no overwrite — folder already exists from TRINITY

    input:
        path(fasta)

    output:
        path "trinity_stats.txt", emit: stats

    script:
    """
    TrinityStats.pl ${fasta} > trinity_stats.txt
    """
}
