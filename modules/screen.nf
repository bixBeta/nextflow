screen = params.screen

process SCREENM {

    label 'process_screen'
    tag "$id"

    publishDir "fq_screen" , overwrite: true

    input:
        tuple val(id), path(trimmed)


    output:
        path '*'

    
    script:

        if ( screen ) {

            fastq_screen --conf ${baseDir}/screen.conf ${trimmed}

        }



}