// mqcgenome =  params.genome 

process MQC {

    label 'process_mqc'
    
    publishDir "Reports", mode: "move", overwrite: true
    input:

        path "*"
        path(conf)
        path(logo)  
        val(mqcgenome)           

    output:
        path "*html"                    , emit: mqc_out  

    when:
        
    script:

    """
       export  MQC_GENOME=${mqcgenome} 
       multiqc -n ${params.id}.star.multiqc.report --config ${conf} --cl-config "custom_logo: ${logo}" -m star .

    """

}


process MQC2 {

    label 'process_mqc'

    publishDir "Reports", mode: "move", overwrite: true
    
    input:

        path "*"              
        path(conf)
        path(logo)  
        val(mqcgenome)
        
    output:
        path "*html"                    , emit: mqc_out2  

    when:
        
    script:

    """
       export  MQC_GENOME=${mqcgenome}
       multiqc -n ${params.id}.starSplit.multiqc.report --config ${conf} --cl-config "custom_logo: ${logo}" -m star .

    """

}


process MQCSCREENM {

    label 'process_mqc'

    publishDir "Reports", mode: "move", overwrite: true
    
    input:

        path "*"     
        path(conf)         
        path(logo)  
        
    output:
        path "*html"                    , emit: mqc_out_screen  

    when:
        
    script:

    """
       multiqc -n ${params.id}.fq.screen.multiqc.report --config ${conf} --cl-config "custom_logo: ${logo}" .

    """

}