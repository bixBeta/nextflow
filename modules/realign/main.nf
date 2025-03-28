runmode         = params.mode
splitmap        = params.genome2
strandedness    = params.strand


process STARM2 {
    maxForks 1
    tag "$id"
    label 'process_high'
    
    publishDir "STAR_BAMS2"   , mode: "symlink", overwrite: true , pattern: "*.bam"
    publishDir "STAR_COUNTS2" , mode: "symlink", overwrite: true , pattern: "*ReadsPerGene.out.tab"

    input:
        tuple val(id), path(unmapped)
        path genome2
        val splitname
        val gkey 
        
    output:
        path "*ReadsPerGene.out.tab"                                        , emit: read_per_gene_tab2 
        path "*Log.final.out"                                               , emit: log_final2
        path "*Log.out"                                                     , emit: log_out2
        path "*Log.progress.out"                                            , emit: log_progress2
        path "*SJ.out.tab"                                                  , emit: sj_out_tab2
        path "*_val_*.fq*"                      , optional:true             , emit: unmapped2
        path "*bam"                                                         , emit: bam_sorted2

    script:

    if (runmode == "SES" & splitmap != null )

        """
            STAR \
            --runThreadN ${task.cpus} \
            --genomeDir  ${genome2} \
            --readFilesIn ${unmapped} \
            --readFilesCommand gunzip -c \
            --outSAMstrandField intronMotif \
            --outFilterIntronMotifs RemoveNoncanonical \
            --outSAMtype BAM SortedByCoordinate \
            --outReadsUnmapped Fastx \
            --outFileNamePrefix ${id}-non-${gkey}-mappedTo-${splitname}. \
            --limitBAMsortRAM 61675612266 \
            --quantMode GeneCounts
        
        mv *.out.mate1 ${id}.non.${gkey}.non.${splitname}_val_1.fq
        gzip *_val_1.fq 

        """

    else if (runmode == "SEBS" & splitmap != null )

        """
            STAR \
            --runThreadN ${task.cpus} \
            --genomeDir  ${genome2} \
            --readFilesIn ${unmapped} \
            --readFilesCommand gunzip -c \
            --outSAMstrandField intronMotif \
            --outFilterIntronMotifs RemoveNoncanonical \
            --outSAMtype BAM SortedByCoordinate \
            --outReadsUnmapped Fastx \
            --outFileNamePrefix ${id}-non-${gkey}-mappedTo-${splitname}. \
            --limitBAMsortRAM 61675612266 \
            --quantMode GeneCounts \
            --alignIntronMax 1 \
            --alignMatesGapMax 45000 

        mv *.out.mate1 ${id}.non.${gkey}.non.${splitname}_val_1.fq
        gzip *_val_1.fq 

        """

    else if (params.mode == "PES" & splitmap != null )
   
        """
            STAR \
            --runThreadN ${task.cpus} \
            --genomeDir ${genome2} \
            --readFilesIn ${unmapped[0]} ${unmapped[1]} \
            --readFilesCommand gunzip -c \
            --outSAMstrandField intronMotif \
            --outFilterIntronMotifs RemoveNoncanonical \
            --outSAMtype BAM SortedByCoordinate \
            --outFileNamePrefix ${id}-non-${gkey}-mappedTo-${splitname}. \
            --limitBAMsortRAM 61675612266 \
            --quantMode GeneCounts \
            --outReadsUnmapped Fastx 

        mv *.out.mate1 ${id}.non.${gkey}.non.${splitname}_val_1.fq
        mv *.out.mate2 ${id}.non.${gkey}.non.${splitname}_val_2.fq

        gzip *_val_1.fq 
        gzip *_val_2.fq 

        """
    else if (params.mode == "PEBS" & splitmap != null  )

        """
             STAR \
            --runThreadN ${task.cpus} \
            --genomeDir ${genome2} \
            --readFilesIn ${unmapped[0]} ${unmapped[1]} \
            --readFilesCommand gunzip -c \
            --outSAMstrandField intronMotif \
            --outFilterIntronMotifs RemoveNoncanonical \
            --outSAMtype BAM SortedByCoordinate \
            --outFileNamePrefix ${id}-non-${gkey}-mappedTo-${splitname}. \
            --limitBAMsortRAM 61675612266 \
            --quantMode GeneCounts \
            --outReadsUnmapped Fastx \
            --alignIntronMax 1 \
            --alignMatesGapMax 45000                   

        mv *.out.mate1 ${id}.non.${gkey}.non.${splitname}_val_1.fq
        mv *.out.mate2 ${id}.non.${gkey}.non.${splitname}_val_2.fq

        gzip *_val_1.fq 
        gzip *_val_2.fq 

        """

    else {
        error "Invalid alignment mode: ${runmode} or star run not specified. To run star use --star in the run command "
        exit 0
    } 
}





process COUNTSM2 {

    publishDir "STAR_COUNTS2/rawCounts2" , mode: "symlink", overwrite: true , pattern: "*.rawCounts"

    input:
        path(counts)

    output:
        path "*.rawCounts"            ,         emit: raw_counts

    script:

        if ( strandedness == 2 )

        """
        BASE=`basename ${counts} .ReadsPerGene.out.tab`
        awk 'NR > 4 {print \$1 "\t" \$4}' ${counts} > \$BASE.rawCounts

        """

        else if ( strandedness == 1 )

        """
        BASE=`basename ${counts} .ReadsPerGene.out.tab`
        awk 'NR > 4 {print \$1 "\t" \$3}' ${counts} > \$BASE.rawCounts
        
        """

        else 

        """
        BASE=`basename ${counts} .ReadsPerGene.out.tab`
        awk 'NR > 4 {print \$1 "\t" \$2}' ${counts} > \$BASE.rawCounts
        
        """


}