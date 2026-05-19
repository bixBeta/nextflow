screen      = params.screen
runmode     = params.mode

process SCREENM {

    label 'process_screen'
    tag "$id"

    publishDir "fq_screen" , overwrite: true

    input:
        tuple val(id), path(trimmed)
        path(screen_conf)

    output:
        path "*"                , emit: screen_out
        path "versions.yml"    , emit: versions



    script:
    

    if ( runmode == "SE" || runmode == "SES" || runmode == "SEBS" ){

    """
     fastq_screen --conf ${screen_conf} ${trimmed[0]}
     mv *screen.txt ${id}_R1_screen.txt
     mv *screen.html ${id}_R1_screen.html

     cat <<-END_VERSIONS > versions.yml
     "${task.process}":
         fastq_screen: \$(fastq_screen --version 2>&1 | head -1 | sed 's/FastQ Screen v//')
     END_VERSIONS
    """

    }

    else if ( runmode == "PE" || runmode == "PES" || runmode == "PEBS" ){

    """
     fastq_screen --conf ${screen_conf} ${trimmed[0]}
     mv *screen.txt ${id}_R1_screen.txt
     mv *screen.html ${id}_R1_screen.html

     cat <<-END_VERSIONS > versions.yml
     "${task.process}":
         fastq_screen: \$(fastq_screen --version 2>&1 | head -1 | sed 's/FastQ Screen v//')
     END_VERSIONS
    """

    }  else {

        error "Runmode ${runmode} is not supported"
        exit 0 
    }  

}