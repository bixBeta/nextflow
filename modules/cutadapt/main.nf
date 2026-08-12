runmode = params.mode

process CUTADAPT {
    maxForks 8
    tag "$id"
    label 'process_high'

    publishDir "trimmed_fastqs", mode: "symlink", overwrite: true
    publishDir "cutadapt_logs", mode: "symlink", overwrite: true, pattern: "*.cutadapt.json"

    input:
        tuple val(id), path(reads)

    output:
        tuple val(id), path("*gz")       , emit: trimmed_fqs
        path("*.cutadapt.json")          , emit: cutadapt_json
        path "versions.yml"              , emit: versions

    script:

    if ( runmode == "SE3PL" )

        """
        R1_raw=${reads}
        R1_trimmed=${id}_val_1.fq.gz

        cutadapt -m 50 -O 20 \
            -a "polyA=A{20}" \
            -a "QUALITY=G{20}" \
            -n 2 \
            \${R1_raw} | \
        cutadapt -m 50 -O 3 \
            --nextseq-trim=10 \
            -a "r1adapter=A{18}AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC;min_overlap=3;max_error_rate=0.100000" \
            - | \
        cutadapt -m 50 -O 20 \
            -g "r1adapter=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC;min_overlap=20" \
            --discard-trimmed \
            --json=${id}.cutadapt.json \
            -o \${R1_trimmed} \
            -

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            cutadapt: \$(cutadapt --version)
        END_VERSIONS
        """

    else {
        error "Runmode ${runmode} is not supported by CUTADAPT"
        exit 0
    }
}
