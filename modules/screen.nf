screen = params.screen

process SCREENM {

    label 'process_screen'
    tag "$id"

    publishDir "fq_screen" , overwrite: true

    input:
        tuple val(id), path(trimmed)


    output:
        path '*'

    when: 
        screen == true    

    script:
    
     """
    
     fastq_screen --conf ${baseDir}/screen.conf ${trimmed}
            
     """
       



}