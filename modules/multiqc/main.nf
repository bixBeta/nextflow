// mqcgenome =  params.genome 

process MQC {

    label 'process_mqc'

    publishDir "Reports/${mqcgenome}", mode: "move", overwrite: true
    input:

        path "*"
        path(conf)
        path(logo)
        val(mqcgenome)
        path(mqc_versions)

    output:
        path "*html"                    , emit: mqc_out
        path "versions.yml"            , emit: versions

    when:

    script:

    """
       export  MQC_GENOME=${mqcgenome}
       multiqc -n ${params.id}.star.multiqc.report --config ${conf} --cl-config "custom_logo: ${logo}"  .

       cat <<-END_VERSIONS > versions.yml
       "${task.process}":
           multiqc: \$(multiqc --version 2>&1 | sed 's/multiqc, version //')
       END_VERSIONS
    """

}


process MQC2 {

    label 'process_mqc'

    publishDir "Reports/${mqcgenome}", mode: "move", overwrite: true

    input:

        path "*"
        path(conf)
        path(logo)
        val(mqcgenome)
        path(mqc_versions)

    output:
        path "*html"                    , emit: mqc_out2

    when:

    script:

    """
       export  MQC_GENOME=${mqcgenome}
       multiqc -n ${params.id}.starSplit.multiqc.report --config ${conf} --cl-config "custom_logo: ${logo}"  .

    """

}


process MQC3 {

    label 'process_mqc'

    publishDir "Reports/trinity", mode: "move", overwrite: true

    input:
        val(salmon_quant_path)
        path(trinity_conf)
        path(logo)
        path(mqc_versions)

    output:
        path "*html", emit: mqc_out3

    script:
    """
    multiqc -n ${params.id}.trinity.multiqc.report \\
        --config ${trinity_conf} \\
        --cl-config "custom_logo: ${logo}" \\
        ${salmon_quant_path}/*/
    """

}


process MQCSCREENM {

    label 'process_mqc'

    publishDir "Reports", mode: "move", overwrite: true

    input:

        path "*"
        path(conf)
        path(logo)
        path(mqc_versions)

    output:
        path "*html"                    , emit: mqc_out_screen

    when:

    script:

    """
       multiqc -n ${params.id}.fq.screen.multiqc.report --config ${conf} --cl-config "custom_logo: ${logo}" .

    """

}