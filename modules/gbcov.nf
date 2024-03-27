gbcovRun = params.gbcov 

process GBCOV1M {
    tag "chr subset - $chromosome"
    label "process_medium"

    // publishDir "GBCOV", mode: "copy", overwrite: true
    
    input:
        path bam
        val chromosome
    
    output:
        path "*.chr*.bam"           , emit: sub_bam
        path "*.chr*.bam.bai"       , emit: sub_bam_index
    
    script:

    if (gbcovRun)
        """
        BASE=`basename \$(echo ${bam}) .Aligned.sortedByCoord.out.bam `

        mv ${bam} \${BASE}.bam

        samtools index \${BASE}.bam

        samtools view -b \${BASE}.bam ${chromosome} > \${BASE}.chr${chromosome}.bam

        samtools index  \${BASE}.chr${chromosome}.bam
        
        """



}



process GBCOV2M {

    tag "$pin"
    label "process_medium"

    publishDir "GBCOV", mode: "copy", overwrite: true
    
    input:
        val pin 
        path "bed"
        val "gbcov"
        
    
    output:
        path "*"             , emit: gbcov_png

    
    script:
    if (gbcovRun)

        b = gbcov.join(",")
        println("samples used for genebody coverage -- ", b )
        """
        geneBody_coverage.py -r ${bed} -i ${b} -o ${pin} 
        """

}
