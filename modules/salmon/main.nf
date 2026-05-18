// ─────────────────────────────────────────────────────────────────────────────
// SALMON_INDEX — build Salmon index from Trinity.fasta (run once)
// ─────────────────────────────────────────────────────────────────────────────
process SALMON_INDEX {
    label 'process_medium'
    publishDir "trinity_assembly/salmon_index"

    input:
        path(fasta)

    output:
        path "trinity_index", emit: index

    script:
    """
    salmon index -t ${fasta} -i trinity_index -p ${task.cpus}
    """
}


// ─────────────────────────────────────────────────────────────────────────────
// SALMON_QUANT — per-sample quantification against Trinity index
// ─────────────────────────────────────────────────────────────────────────────
process SALMON_QUANT {
    label 'process_medium'
    tag "$id"
    publishDir "trinity_assembly/salmon_quant", mode: 'copy'

    input:
        tuple val(id), path(reads)
        path(index)

    output:
        path "${id}", emit: quant_dir

    script:
    def r1    = reads instanceof List ? reads[0] : reads
    def r2    = reads instanceof List && reads.size() > 1 ? reads[1] : null
    def mates = r2 ? "-1 ${r1} -2 ${r2}" : "-r ${r1}"
    """
    salmon quant -i ${index} -l A \\
        ${mates} \\
        -p ${task.cpus} \\
        --validateMappings \\
        -o ${id}
    """
}
