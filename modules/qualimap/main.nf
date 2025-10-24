
process QUALIMAP {

    maxForks 4
    tag "$bam"
    label "process_qualimap"

   
    publishDir "STATS/QUALIMAP_RES",             mode: "symlink", overwrite: true, pattern: "*bamqc"
    

    input:
         path(bam)

    output:
        path("*bamqc")                   , emit: bamqc_out
        
    script:

        
        """
        id=`basename ${bam} .Aligned.sortedByCoord.out.bam`
            qualimap bamqc -bam ${bam} \\
            -nt 16 -c \\
            --java-mem-size=30G \\
            -outdir \${id}.bamqc
        
        sed -i "s/bam file = \${id}.dupMarked.bam/bam file = \${id}.bam/g" \${id}.bamqc/genome_results.txt

        """

}
