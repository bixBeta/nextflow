runmode     = params.mode
splitmap    = params.genome2



process STARM2 {
    maxForks 1
    tag "$id"
    label 'process_high'
    
    // publishDir "$baseDir/STAR_OUT", mode: "copy", overwrite: false
    
    input:
        tuple val(id), path(unmapped)
        path genome2

    output:
        path "*ReadsPerGene.out.tab"                                        , emit: read_per_gene_tab2 
        path "*Log.final.out"                                               , emit: log_final2
        path "*Log.out"                                                     , emit: log_out2
        path "*Log.progress.out"                                            , emit: log_progress2
        path "*SJ.out.tab"                                                  , emit: sj_out_tab2
        path "*.out.mate*"                      , optional:true             , emit: unmapped2
        path "*bam"                                                         , emit: bam_sorted2

    script:

    if (runmode == "SES" & splitmap != null )

        """
            STAR \
            --runThreadN ${task.cpus} \
            --genomeDir  ${genome2} \
            --readFilesIn ${trimmed} \
            --readFilesCommand gunzip -c \
            --outSAMstrandField intronMotif \
            --outFilterIntronMotifs RemoveNoncanonical \
            --outSAMtype BAM SortedByCoordinate \
            --outReadsUnmapped Fastx \
            --outFileNamePrefix ${id}.splitMappedTo.${splitmap} \
            --limitBAMsortRAM 61675612266 \
            --quantMode GeneCounts
        
        gzip *.out.mate1 

        """

    else if (runmode == "SEBS" & splitmap != null )

        """
            STAR \
            --runThreadN ${task.cpus} \
            --genomeDir  ${genome2} \
            --readFilesIn ${trimmed} \
            --readFilesCommand gunzip -c \
            --outSAMstrandField intronMotif \
            --outFilterIntronMotifs RemoveNoncanonical \
            --outSAMtype BAM SortedByCoordinate \
            --outReadsUnmapped Fastx \
            --outFileNamePrefix ${id}.splitMappedTo.${splitmap} \
            --limitBAMsortRAM 61675612266 \
            --quantMode GeneCounts \
            --alignIntronMax 1 \
            --alignMatesGapMax 45000 

        gzip *.out.mate1

        """

    else if (params.mode == "PES" & splitmap != null )
   
        """
            STAR \
            --runThreadN ${task.cpus} \
            --genomeDir ${genome2} \
            --readFilesIn ${trimmed[0]} ${trimmed[1]} \
            --readFilesCommand gunzip -c \
            --outSAMstrandField intronMotif \
            --outFilterIntronMotifs RemoveNoncanonical \
            --outSAMtype BAM SortedByCoordinate \
            --outFileNamePrefix ${id}.splitMappedTo.${splitmap} \
            --limitBAMsortRAM 61675612266 \
            --quantMode GeneCounts \
            --outReadsUnmapped Fastx 

        gzip *.out.mate1
        gzip *.out.mate2

        """
    else if (params.mode == "PEBS" & splitmap != null  )

        """
             STAR \
            --runThreadN ${task.cpus} \
            --genomeDir ${genome2} \
            --readFilesIn ${trimmed[0]} ${trimmed[1]} \
            --readFilesCommand gunzip -c \
            --outSAMstrandField intronMotif \
            --outFilterIntronMotifs RemoveNoncanonical \
            --outSAMtype BAM SortedByCoordinate \
            --outFileNamePrefix ${id}.splitMappedTo.${splitmap} \
            --limitBAMsortRAM 61675612266 \
            --quantMode GeneCounts \
            --outReadsUnmapped Fastx \
            --alignIntronMax 1 \
            --alignMatesGapMax 45000                   

        gzip *.out.mate1
        gzip *.out.mate2

        """

    else {
        error "Invalid alignment mode: ${runmode} or star run not specified. To run star use --star in the run command "
        exit 0
    } 
}


