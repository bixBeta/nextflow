runmode = params.mode

process FASTPM {
    maxForks 8
    tag "$id"
    label 'process_high'
    
    publishDir "trimmed_fastqs", mode: "symlink", overwrite: true
    publishDir "fastp_logs", mode: "symlink", overwrite: true, pattern:"*.fastp.json"

    input:
        tuple val(id), path(reads)
    
    output:
        tuple val(id), path("*gz")              , emit: trimmed_fqs
        path("*.fastp.json")                    , emit: fastp_json
        path "versions.yml"                     , emit: versions

    script:

    if ( runmode == "SE" || runmode == "SES" || runmode == "SEBS" || runmode == "SEB" ){
        
        """
        fastp \
        -z 4 -w 16 \
        --length_required 50 --qualified_quality_phred 20 \
        --trim_poly_g \
        -i ${reads} \
        -o ${id}_val_1.fq.gz \
        -h ${id}.fastp.html \
        -j ${id}.fastp.json

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastp: \$(fastp --version 2>&1 | head -1 | sed 's/fastp //')
        END_VERSIONS
        """

    }

    else if ( runmode == "PE" || runmode == "PES" || runmode == "PEBS" || runmode == "PEB" ){

        """
            fastp \
            -z 4 -w 16 \
            --length_required 50 --qualified_quality_phred 20 \
            --trim_poly_g \
            -i ${reads[0]} \
            -I ${reads[1]} \
            -o ${id}_val_1.fq.gz \
            -O ${id}_val_2.fq.gz \
            -h ${id}.fastp.html \
            -j ${id}.fastp.json

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastp: \$(fastp --version 2>&1 | head -1 | sed 's/fastp //')
        END_VERSIONS
        """


    } else {

        error "Runmode ${runmode} is not supported"
        exit 0
    }



}

