runmode         = params.mode
strandedness    = params.strand


process STARM {
    maxForks 3
    tag "$id"
    label 'process_high'
    
    publishDir "STAR_BAMS/${mqcgenome}"   , mode: "symlink", overwrite: true , pattern: "*.bam"
    publishDir "STAR_COUNTS/${mqcgenome}" , mode: "symlink", overwrite: true , pattern: "*ReadsPerGene.out.tab"
    
    input:
        tuple val(id), path(trimmed)
        path genome
        val(mqcgenome)
        

    output:
        path "*ReadsPerGene.out.tab"                                        , emit: read_per_gene_tab 
        path "*Log.final.out"                                               , emit: log_final
        path "*Log.out"                                                     , emit: log_out
        path "*Log.progress.out"                                            , emit: log_progress
        path "*SJ.out.tab"                                                  , emit: sj_out_tab
        path "*bam"                                                         , emit: bam_sorted
        tuple val(id), path("*_val_*.fq*")  ,     optional:true             , emit: unmapped
        path "versions.yml"                                                 , emit: versions

    script:

    if (runmode == "SE3PL" )
        """
        STAR \
            --runThreadN ${task.cpus} \
            --genomeDir ${genome} \
            --readFilesIn ${trimmed} \
            --readFilesCommand zcat \
            --outSAMtype BAM SortedByCoordinate \
            --outFilterType BySJout \
            --outFilterMultimapNmax 20 \
            --alignSJoverhangMin 8 \
            --alignSJDBoverhangMin 1 \
            --outFilterMismatchNmax 999 \
            --outFilterMismatchNoverReadLmax 0.04 \
            --alignIntronMin 20 \
            --alignIntronMax 1000000 \
            --outFileNamePrefix ${id}. \
            > ${id}.star.log 2>&1

        samtools index ${id}.Aligned.sortedByCoord.out.bam

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            STAR: \$(STAR --version)
            samtools: \$(samtools --version | head -1 | sed 's/samtools //')
        END_VERSIONS
        """

    else if (runmode == "SE" )
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

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            STAR: \$(STAR --version)
        END_VERSIONS
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

        BASE=`basename ${mqcgenome}`
        mv *.out.mate1 ${id}.non.\${BASE}_val_1.fq
        gzip *_val_1.fq

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            STAR: \$(STAR --version)
        END_VERSIONS
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

        BASE=`basename ${mqcgenome}`
        mv *.out.mate1 ${id}.non.\${BASE}_val_1.fq
        gzip *_val_1.fq

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            STAR: \$(STAR --version)
        END_VERSIONS
        """

    else if (runmode == "SEB"  )

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
            --quantMode GeneCounts \
            --alignIntronMax 1 \
            --alignMatesGapMax 45000

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            STAR: \$(STAR --version)
        END_VERSIONS
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

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            STAR: \$(STAR --version)
        END_VERSIONS
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

        BASE=`basename ${mqcgenome}`
        mv *.out.mate1 ${id}.non.\${BASE}_val_1.fq
        mv *.out.mate2 ${id}.non.\${BASE}_val_2.fq

        gzip *_val_1.fq
        gzip *_val_2.fq

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            STAR: \$(STAR --version)
        END_VERSIONS
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

        BASE=`basename ${mqcgenome}`
        mv *.out.mate1 ${id}.non.\${BASE}_val_1.fq
        mv *.out.mate2 ${id}.non.\${BASE}_val_2.fq

        gzip *_val_1.fq
        gzip *_val_2.fq

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            STAR: \$(STAR --version)
        END_VERSIONS
        """

    else if (params.mode == "PEB"  )

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
            --alignIntronMax 1 \
            --alignMatesGapMax 45000

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            STAR: \$(STAR --version)
        END_VERSIONS
        """


    else {
        error "Invalid alignment mode: ${runmode} "
        exit 0
    } 
}



process COUNTSM {

    publishDir "STAR_COUNTS/${mqcgenome}/rawCounts" , mode: "symlink", overwrite: true , pattern: "*.rawCounts"

    input:
        path(counts)
        val(mqcgenome)

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
