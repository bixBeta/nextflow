screen = params.screen

process SCREENM {

    label 'process_medium'
    tag "$id"

    publishDir "fq_screen" , overwrite: true

    input:
        tuple val(id), path(fq_screen)
        path screenConf


    output:
        path '*'

    
    script:

        if ( screen ) {

            fastq_screen --conf ${projectDir}/screen.conf ${fq_screen}

        }



}