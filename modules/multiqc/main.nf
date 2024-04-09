process MQC {

    label 'process_mqc'
    
    publishDir "Reports", mode: "move", overwrite: false
    input:

        path "*"              

    output:
        path "*html"                    , emit: mqc_out  

    when:
        
    script:

    """
       multiqc -n ${params.id}.star.multiqc.report --config ${baseDir}/multiqc_config.yaml -m star .

    """

}


process MQC2 {

    label 'process_mqc'

    publishDir "Reports", mode: "move", overwrite: true
    
    input:

        path "*"              

    output:
        path "*html"                    , emit: mqc_out2  

    when:
        
    script:

    """
       multiqc -n ${params.id}.starSplit.multiqc.report --config multiqc_config.yaml -m star .

    """

}


process MQCSCREENM {

    label 'process_mqc'

    publishDir "Reports", mode: "move", overwrite: true
    
    input:

        path "*"              

    output:
        path "*html"                    , emit: mqc_out_screen  

    when:
        
    script:

    """
       multiqc -n ${params.id}.fq.screen.multiqc.report --config multiqc_config.yaml .

    """

}