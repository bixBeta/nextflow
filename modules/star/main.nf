runmode         = params.mode
gkey            = params.genome
strandedness    = params.strand


process STARM {
    maxForks 1
    tag "$id"
    label 'process_high'
    
    publishDir "STAR_BAMS"   , mode: "symlink", overwrite: true , pattern: "*.bam"
    publishDir "STAR_COUNTS" , mode: "symlink", overwrite: true , pattern: "*ReadsPerGene.out.tab"
    
    input:
        tuple val(id), path(trimmed)
        path genome
        

    output:
        path "*ReadsPerGene.out.tab"                                        , emit: read_per_gene_tab 
        path "*Log.final.out"                                               , emit: log_final
        path "*Log.out"                                                     , emit: log_out
        path "*Log.progress.out"                                            , emit: log_progress
        path "*SJ.out.tab"                                                  , emit: sj_out_tab
        path "*bam"                                                         , emit: bam_sorted
        tuple val(id), path("*_val_*.fq*")  ,     optional:true             , emit: unmapped
    
    script:

    if (runmode == "SE" )
        """
            STAR \
            --runThreadN ${task.cpus} \
            --genomeDir  ${genome} \
            --readFilesIn ${trimmed} \
            --readFilesCommand gunzip -c \
            --outSAMstrandField intronMotif \
            --outFilterIntronMotifs RemoveNoncanonical \
            --outSAMtype BAM SortedByCoordinate \
            --outFileNamePrefix ${id}. \
            --limitBAMsortRAM 61675612266 \
            --quantMode GeneCounts

        """
 
 
    else if (runmode == "SES"  )

        """
            STAR \
            --runThreadN ${task.cpus} \
            --genomeDir  ${genome} \
            --readFilesIn ${trimmed} \
            --readFilesCommand gunzip -c \
            --outSAMstrandField intronMotif \
            --outFilterIntronMotifs RemoveNoncanonical \
            --outSAMtype BAM SortedByCoordinate \
            --outReadsUnmapped Fastx \
            --outFileNamePrefix ${id}. \
            --limitBAMsortRAM 61675612266 \
            --quantMode GeneCounts
        
        BASE=`basename ${gkey}`
        mv *.out.mate1 ${id}.non.\${BASE}_val_1.fq
        gzip *_val_1.fq
        
        """

    else if (runmode == "SEBS"  )

        """
            STAR \
            --runThreadN ${task.cpus} \
            --genomeDir  ${genome} \
            --readFilesIn ${trimmed} \
            --readFilesCommand gunzip -c \
            --outSAMstrandField intronMotif \
            --outFilterIntronMotifs RemoveNoncanonical \
            --outSAMtype BAM SortedByCoordinate \
            --outReadsUnmapped Fastx \
            --outFileNamePrefix ${id}. \
            --limitBAMsortRAM 61675612266 \
            --quantMode GeneCounts \
            --alignIntronMax 1 \
            --alignMatesGapMax 45000 

        BASE=`basename ${gkey}`
        mv *.out.mate1 ${id}.non.\${BASE}_val_1.fq
        gzip *_val_1.fq

        """

    else if (params.mode == "PE"  )
       
        """
            STAR \
            --runThreadN ${task.cpus} \
            --genomeDir ${genome} \
            --readFilesIn ${trimmed[0]} ${trimmed[1]} \
            --readFilesCommand gunzip -c \
            --outSAMstrandField intronMotif \
            --outFilterIntronMotifs RemoveNoncanonical \
            --outSAMtype BAM SortedByCoordinate \
            --outFileNamePrefix ${id}. \
            --limitBAMsortRAM 61675612266 \
            --quantMode GeneCounts

        """

    else if (params.mode == "PES"  )
   
        """
            STAR \
            --runThreadN ${task.cpus} \
            --genomeDir ${genome} \
            --readFilesIn ${trimmed[0]} ${trimmed[1]} \
            --readFilesCommand gunzip -c \
            --outSAMstrandField intronMotif \
            --outFilterIntronMotifs RemoveNoncanonical \
            --outSAMtype BAM SortedByCoordinate \
            --outFileNamePrefix ${id}. \
            --limitBAMsortRAM 61675612266 \
            --quantMode GeneCounts \
            --outReadsUnmapped Fastx 


        BASE=`basename ${gkey}`
        mv *.out.mate1 ${id}.non.\${BASE}_val_1.fq
        mv *.out.mate2 ${id}.non.\${BASE}_val_2.fq

        gzip *_val_1.fq
        gzip *_val_2.fq
       
        """
    
    
    else if (params.mode == "PEBS"  )

        """
             STAR \
            --runThreadN ${task.cpus} \
            --genomeDir ${genome} \
            --readFilesIn ${trimmed[0]} ${trimmed[1]} \
            --readFilesCommand gunzip -c \
            --outSAMstrandField intronMotif \
            --outFilterIntronMotifs RemoveNoncanonical \
            --outSAMtype BAM SortedByCoordinate \
            --outFileNamePrefix ${id}. \
            --limitBAMsortRAM 61675612266 \
            --quantMode GeneCounts \
            --outReadsUnmapped Fastx \
            --alignIntronMax 1 \
            --alignMatesGapMax 45000     


        BASE=`basename ${gkey}`
        mv *.out.mate1 ${id}.non.\${BASE}_val_1.fq
        mv *.out.mate2 ${id}.non.\${BASE}_val_2.fq

        gzip *_val_1.fq
        gzip *_val_2.fq

        """

    else {
        error "Invalid alignment mode: ${runmode} "
        exit 0
    } 
}



process COUNTSM {

    publishDir "STAR_COUNTS/rawCounts" , mode: "symlink", overwrite: true , pattern: "*.rawCounts"

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