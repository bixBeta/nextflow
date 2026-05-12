runmode = params.mode

process TRINITY {

    maxForks 1
    label 'process_high'
    tag "$id"

    publishDir "trinity_assembly" , overwrite: true

    input:
        tuple val(id), path(r1_reads, stageAs: 'r1/*'), path(r2_reads, stageAs: 'r2/*')
        val(libtype)

    output:
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
