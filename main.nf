nextflow.enable.dsl=2

params.sheet            = "sample-sheet.csv"
params.outdir           = "$projectDir/STAR_OUT"
params.reads            = "$workDir/fastqs/*_*{1,2}.f*.gz"
params.help             = false
params.listGenomes      = false
params.genome           = null
params.mode             = "PE"
params.id               = "TREx_ID"
params.gbcov            = false
params.chromosub        = "10"
params.genome2          = null
params.splitname        = "na"
params.screen           = false
params.strand           = 2

runmode = params.mode
pin = channel.value(params.id)

if( params.help ) {

log.info """
R  N  A  -  S  E  Q      W  O  R  K  F  L  O  W  -  @bixBeta
=======================================================================================================================================================================
Usage:
    nextflow run https://github.com/bixbeta/nextflow -r main < args ... >

Args:
    * --listGenomes    : Get extended list of genomes available for this pipeline
    * --id             : TREx Project ID 
    * --sheet          : sample-sheet.csv < default: looks for a file named sample-sheet.csv in the project dir >

        -------------------------------------------
        Sample Sheet Example:    
        label   fastq1          fastq2
        SS1     SS1_R1.fastq.gz SS1_R2.fastq.gz
        SS2     SS2_R1.fastq.gz SS2_R2.fastq.gz  
        .
        .
        . etc.
        -------------------------------------------

    * --genome          : Genome index. Use --listGenomes flag to see all available genomes. Also supports a path value for starIndex dir. 
    * --genome2         : Secondary Genome index. This will align the --genome subtracted reads to --genome2 index
    * --screen          : Invokes the fastq_screen step. See the screen.conf file here <https://github.com/bixBeta/nextflow/blob/main/screen.conf> for more details 
    * --gbcov           : Runs GeneBodyCoverage Program on sub-setted bams.
    * --chromosub       : Subset bams to specified chromosome name. < defaults to chromosome 10 >
    * --splitname       : A string that will be used to denote --genome2 e.g. "GRC100011A", "Cat_custom" etc. 

"""

    exit 0
}




log.info """
R  N  A -  S  E  Q      W  O  R  K  F  L  O  W  -  @bixBeta  
=========================================================================================================================
trexID       : ${params.id}
reads        : ${params.reads}
genome       : ${params.genome}
mode         : ${params.mode}
genome2      : ${params.genome2}
screen       : ${params.screen}
gbcov        : ${params.gbcov}
chromosub    : ${params.chromosub}

"""

// STAR index MAP
genomeDir = [
hg38				:"/workdir/genomes/Homo_sapiens/hg38/UCSC/hg38.star",
mm10				:"/workdir/genomes/Mus_musculus/mm10/UCSC/mm10.star",
GRCh38				:"/workdir/genomes/Homo_sapiens/hg38/ENSEMBL/GRCh38.star",
GRCm38				:"/workdir/genomes/Mus_musculus/mm10/ENSEMBL/GRCm38.star",
cat					:"/workdir/genomes/Felis_catus/Felis_catus9.0/Ensembl/genomeDir",
chicken				:"/workdir/genomes/Gallus_gallus/Galgal5/ENSEMBL/galgal5.star",
chicken6			:"/workdir/genomes/Gallus_gallus/Galgal6/ENSEMBL/gtf.102/star.index.102",
horse				:"/workdir/genomes/Equus_caballus/EquCab3/ENSEMBL/genebuild-102/star.index.102",
horse2				:"/workdir/genomes/Equus_caballus/EquCab2/ENSEMBL/EquCab2.star.index",
ATCC_13047			:"/workdir/genomes/Enterobacter_cloacae/ATCC_13047/custom/ATCC_13047.GTF",
grape				:"/workdir/genomes/Vitis_vinifera/GCA_000003745.2/ENSEMBL/Vitis_vinifera.12X.43.bed12",
rat					:"/workdir/genomes/Rattus_norvegicus/rn6/ENSEMBL/rat.star",
ercc				:"/workdir/genomes/contaminants/ERCC_spikeIns/ercc.star",
lonchura			:"/workdir/genomes/Lonchura_striata/LonStrDom1/ENSEMBL/lonchura.star",
goose				:"/workdir/genomes/Anser_brachyrhynchus/ASM259213v1/ENSEMBL/goose.star",
ehv8				:"/workdir/genomes/FastQ_Screen_Genomes/EHV8/ehv8.star",
erdman				:"/workdir/genomes/Mycobacterium_tuberculosis/Erdman_GCA_000350205.1/ENSEMBL/genomeDir",
TB					:"/workdir/genomes/Mycobacterium_tuberculosis/CDC1551_Ensembl/cdc1551.star",
TB2					:"/workdir/genomes/Mycobacterium_tuberculosis/Ensembl_GCA_000668235/GCA_000668235.star",
maize4				:"/workdir/genomes/Zea_mays/B73_RefGen_v4/ENSEMBL/star.maize",
maize3				:"/workdir/genomes/Zea_mays/B73_RefGen_v3/NCBI/genomeDir",
finch				:"/workdir/genomes/Taeniopygia_guttata/taeGut3.2.4/ENSEMBL/UPDATED.ANNOTS/star.index.updated",
finch2				:"/workdir/genomes/Geospiza_fortis_ground_finch/GeoFor_1.0/NCBI/star.index",
green				:"/workdir/genomes/Chlorocebus_sabaeus/CHlSab1/ENSEMBL/genomeDir",
dog					:"/workdir/genomes/Canis_familiaris/canFam3/ENSEMBL/CanFam3.1.101/genomeDir",
faba				:"/workdir/genomes/Vicia_faba/VfEP_Reference-Unigene/NCBI/genomeDir-wo-gff",
aphid				:"/workdir/genomes/Acyrthosiphon_pisum/pea_aphid_22Mar2018_4r6ur/NCBI/genomeDir",
cholera				:"/workdir/genomes/Vibrio_cholerae/N16961/NCBI/genomeDir",
BG8					:"/workdir/genomes/Methylomicrobium_album_BG8/ASM21427v3/NCBI/custom3-CDS2exon/genomeDir_fixSpaces",
aedes				:"/workdir/genomes/Aedes_aegypti/AaegL5.0/NCBI/AaegL5.0.star",
aedesVB				:"/workdir/genomes/Aedes_aegypti/AaegL5.0/Vectorbase/AaegyptiLVP_AGWG_release49/genomeDir.vectorBase",
DC3000c				:"/workdir/genomes/Pseudomonas_syringae/Tomato_DC3000/Custom_Zichu/starIndex_8.gtf_simplenames",
Theileria			:"/workdir/genomes/Theileria/Tannulata_ASM322v1/NCBI/starIndex",
cow					:"/workdir/genomes/Bos_taurus/ARS-UCD1.2_GCA_002263795.2/ENSEMBL/star.index",
salmonella			:"/workdir/genomes/Salmonella_enterica/ASM21085v2/NCBI/genomeDir",
salmonella2			:"/workdir/genomes/Salmonella_enterica/ASM2216v1/NCBI/genomeDir",
EA273				:"/workdir/genomes/Erwinia_amylovora/GCF_000091565.1/ncbi/genomeDir",
ddSmed				:"/workdir/genomes/Schmidtea_mediterranea/dd_Smed_v6/NCBI/genomeDir",
yeast				:"/workdir/genomes/Saccharomyces_cerevisiae/R64-1-1_GCA_000146045.2/ENSEMBL/star.index",
TAIR10				:"/workdir/genomes/Arabidopsis_thaliana/TAIR10/ENSEMBL/genomeDir",
crow				:"/workdir/genomes/Corvus_moneduloides/bCorMon1/NCBI/genomeDir",
orbicella			:"/workdir/genomes/Orbicella_faveolata/GCA_002042975.1/ncbi/genomeDir",
bacillus			:"/workdir/genomes/Bacillus_subtilis/GCA_000009045/ENSEMBL/genomeDir",
pao1				:"/workdir/genomes/Pseudomonas_aeruginosa/PAO1/NCBI/genomeDir",
grape				:"/workdir/genomes/Vitis_vinifera/GCA_000003745.2/ENSEMBL/Vitis_vinifera.star",
cowUMD				:"/workdir/genomes/Bos_taurus/UMD_3.1.1/NCBI/genomeDir/",
guppy				:"/workdir/genomes/Poecilia_reticulata/GCF_000633615.1/NCBI/genomeDir",
dm6					:"/workdir/genomes/Drosophila_melanogaster/dm6/ENSEMBL/genomeDir",
tomato				:"/workdir/genomes/Solanum_lycopersicum/custom/genomeDir",
macaca				:"/workdir/genomes/Macaca_fascicularis/GCF_000364345.1_Macaca_fascicularis_5.0/NCBI/genomeDir"]

// BED12's for RSEQC geneBody coverage process 
bed12 = [
hg38				:"/workdir/genomes/Homo_sapiens/hg38/UCSC/genes.bed12",
mm10				:"/workdir/genomes/Mus_musculus/mm10/UCSC/BED12/mm10.ucsc.bed12",
GRCh38				:"/workdir/genomes/Homo_sapiens/hg38/ENSEMBL/Homo_sapiens.GRCh38.bed12",
GRCm38				:"/workdir/genomes/Mus_musculus/mm10/ENSEMBL/Mus_musculus.GRCm38.bed12",
cat					:"/workdir/genomes/Felis_catus/Felis_catus9.0/Ensembl/Felis_catus.Felis_catus_9.0.95.bed12",
chicken				:"/workdir/genomes/Gallus_gallus/Galgal5/ENSEMBL/Gallus_gallus.Gallus_gallus-5.0.bed12",
chicken6			:"/workdir/genomes/Gallus_gallus/Galgal6/ENSEMBL/gtf.102/Gallus_gallus.GRCg6a.102.bed12",
horse				:"/workdir/genomes/Equus_caballus/EquCab3/ENSEMBL/genebuild-102/Equus_caballus.EquCab3.0.102.bed12",
horse2				:"/workdir/genomes/Equus_caballus/EquCab2/ENSEMBL/Equus_caballus.EquCab2.94.BED12",
ATCC_13047			:"/workdir/genomes/Enterobacter_cloacae/ATCC_13047/GCF_000025565.1_ASM2556v1_genomic.bed12",
grape				:"/workdir/genomes/Vitis_vinifera/GCA_000003745.2/ENSEMBL/Vitis_vinifera.star",
rat					:"/workdir/genomes/Rattus_norvegicus/rn6/ENSEMBL/Rattus_norvegicus.Rnor_6.0.bed12",
lonchura			:"/workdir/genomes/Lonchura_striata/LonStrDom1/ENSEMBL/Lonchura_striata_domestica.LonStrDom1.bed12",
goose				:"/workdir/genomes/Anser_brachyrhynchus/ASM259213v1/ENSEMBL/Anser_brachyrhynchus.ASM259213v1.bed12",
maize4				:"/workdir/genomes/Zea_mays/B73_RefGen_v4/ENSEMBL/Zea_mays.B73_RefGen_v4.bed12",
finch				:"/workdir/genomes/Taeniopygia_guttata/taeGut3.2.4/ENSEMBL/UPDATED.ANNOTS/Taeniopygia_guttata.bTaeGut1_v1.p.bed12",
finch2				:"/workdir/genomes/Geospiza_fortis_ground_finch/GeoFor_1.0/NCBI/GCF_000277835.1_GeoFor_1.0_genomic.bed12",
dog					:"/workdir/genomes/Canis_familiaris/canFam3/ENSEMBL/CanFam3.1.101/Canis_lupus_familiaris.CanFam3.1.101.bed12",
aedes				:"/workdir/genomes/Aedes_aegypti/AaegL5.0/NCBI/GCF_002204515.2_AaegL5.0_genomic.bed12",
aedesVB				:"/workdir/genomes/Aedes_aegypti/AaegL5.0/Vectorbase/AaegyptiLVP_AGWG_release49/VectorBase-49_AaegyptiLVP_AGWG.BED12",
cow					:"/workdir/genomes/Bos_taurus/ARS-UCD1.2_GCA_002263795.2/ENSEMBL/Bos_taurus.ARS-UCD1.2.BED12",
BG8					:"/workdir/genomes/Methylomicrobium_album_BG8/ASM21427v3/NCBI/test/ncbi-genomes-2020-08-27/GCF_000214275.2_ASM21427v3_genomic.bed12",
ddSmed				:"/workdir/genomes/Schmidtea_mediterranea/dd_Smed_v6/NCBI/dd_Smed_v6.pcf.bed12",
cowUMD				:"/workdir/genomes/Bos_taurus/UMD_3.1.1/NCBI/GCF_000003055.6_Bos_taurus_UMD_3.1.1_genomic.bed12",
dm6					:"/workdir/genomes/Drosophila_melanogaster/dm6/ENSEMBL/Drosophila_melanogaster.BDGP6.32.bed12",
tomato				:"/workdir/genomes/Solanum_lycopersicum/custom/ITAG4.0_gene_models.bed12",
yeast				:"/workdir/genomes/Saccharomyces_cerevisiae/R64-1-1_GCA_000146045.2/ENSEMBL/Saccharomyces_cerevisiae.R64-1-1.bed12"]


// genomeKey.list = (genomeDir.collectMany{ k,v -> (v == params.genome) ? [k] : []} as String[])
// genomeKey = genomeKey.list[0]



if( params.listGenomes) {
    
    println("")
    log.info """
    Available STAR Indices
    =========================================================================================================================
    """
    .stripIndent()

    printMap = { a, b -> println "$a ----------- $b" }
    genomeDir.each(printMap)

    log.info """
    BED12's for GeneBodyCov Process 
    =========================================================================================================================
    """
    .stripIndent()

    printMap = { a, b -> println "$a ----------- $b" }
    bed12.each(printMap)

    exit 0
}

include {   FASTPM                   } from './modules/fastp'
include {   STARM ; COUNTSM          } from './modules/star'
include {   GBCOV1M ; GBCOV2M        } from './modules/gbcov'
include {   STARM2 ; COUNTSM2         } from './modules/realign'
include {   MQC ; MQC2 ; MQCSCREENM  } from './modules/multiqc'
include {   SCREENM                  } from './modules/screen'

ch_sheet = channel.fromPath(params.sheet)

ch_mqc_conf = channel.fromPath("${projectDir}/multiqc_config.yaml")
ch_mqc_logo = channel.fromPath("${projectDir}/img/trex-extended-logo.png")

ch_screen_conf = channel.fromPath("${projectDir}/screen.conf")
splitName = params.splitname

if (genomeDir.containsKey(params.genome)){  // allows a user to pass a STAR index path via --genome parameter

    genome = genomeDir[params.genome]

} else {

    genome = params.genome
}

genome_ch = channel.value(genome)
// genome_key_ch = channel.value(params.genome)

if (bed12.containsKey(params.genome)){  // allows a user to pass a STAR index path via --genome parameter

    bed = bed12[params.genome]

} else {

    bed = null
    
}

bed_ch = channel.value(bed)


if (genomeDir.containsKey(params.genome2)){  // allows a user to pass a STAR index path via --genome parameter

    genome2 = genomeDir[params.genome2]

    splitName_ch = channel.value(params.genome2)

} else {

    genome2 = params.genome2
    splitName_ch = channel.value(params.splitname)
}

genome_ch2 = channel.value(genome2)



/* ---------------------------------------------------------------------------------------------------------
SINGLE END NOVA/NEXT-seq Workflow 
------------------------------------------------------------------------------------------------------------ */

workflow SINGLE {

    meta_ch = ch_sheet
                |  splitCsv( header:true )
                |  map { row -> [row.label, [file("fastqs/" + row.fastq1)]] }
                |  view
    

    FASTPM(meta_ch)
        .set { fastp_out }



    SCREENM(fastp_out, ch_screen_conf)

    screen_out_ch = SCREENM.out 
                        | collect
    
    MQCSCREENM(screen_out_ch, ch_mqc_conf, ch_mqc_logo)

    if( params.genome != null ){
    STARM(fastp_out, genome_ch)

    bam_ch = STARM.out.bam_sorted 
                | collect
                | flatten

    unmapped_ch = STARM.out.unmapped

    // chromo_sub = channel.value(params.gbcov)
     chromo_sub = channel.value(params.chromosub)
    }

    if ( params.gbcov & bed != null ) {
        GBCOV1M(bam_ch, chromo_sub)
            .set { gbcov1 }

        gbcov1_ch = GBCOV1M.out.sub_bam
                        //.concat(gbcov1.sub_bam_index)
                        .collect(flat : false)                    
                        // .map { it -> [it + it] }
                        
                        //.view()                    
        

        GBCOV2M(pin, bed_ch, gbcov1_ch)
    }

    if( params.genome != null ){
        mqc_ch1 = STARM.out.read_per_gene_tab
                .concat(STARM.out.log_final)
                .collect()
                //.view()
    

        MQC(mqc_ch1, ch_mqc_conf, ch_mqc_logo)

        COUNTSM(STARM.out.read_per_gene_tab.collect().flatten())
    }


    if( params.genome2 != null ){

        STARM2(unmapped_ch, genome_ch2, splitName_ch )

        mqc_ch2 = STARM2.out.read_per_gene_tab2
                        .concat(STARM2.out.log_final2)
                        .concat(mqc_ch1)
                        .collect()
                        //.view()

        // mqc_ch3 = mqc_ch1.concat(mqc_ch2).view()
    

        MQC2(mqc_ch2, ch_mqc_conf, ch_mqc_logo)
        COUNTSM2(STARM2.out.read_per_gene_tab2.collect().flatten())
    }






}


/* ---------------------------------------------------------------------------------------------------------
PAIRED END NOVA/NEXT-seq Workflow
------------------------------------------------------------------------------------------------------------ */


workflow PAIRED {
    meta_ch = ch_sheet
                |  splitCsv( header:true )
                |  map { row -> [row.label, [file("fastqs/" + row.fastq1), file("fastqs/" + row.fastq2)]] }
                |  view


    FASTPM(meta_ch)
        .set { fastp_out }
        
    // fastp_out.view()

    screen_input_ch = FASTPM.out.trimmed_reads 
                        | view
                                   

   SCREENM(screen_input_ch, ch_screen_conf)

   screen_out_ch = SCREENM.out 
                       | collect

   MQCSCREENM(screen_out_ch, ch_mqc_conf, ch_mqc_logo)
    
    if( params.genome != null ){
        STARM(fastp_out, genome_ch)
    
    bam_ch = STARM.out.bam_sorted 
                | collect
                | flatten
    unmapped_ch = STARM.out.unmapped
    
    // chromo_sub = channel.value(params.gbcov)
       chromo_sub = channel.value(params.chromosub)
    
    }

    if ( params.gbcov  & bed != null ) {
        GBCOV1M(bam_ch, chromo_sub)
            .set { gbcov1 }

        gbcov1_ch = GBCOV1M.out.sub_bam
                        .collect(flat : false)                    
                        //.view()                    
        
        GBCOV2M(pin, bed_ch, gbcov1_ch)
    }

    if( params.genome != null ){
        mqc_ch1 = STARM.out.read_per_gene_tab
                .concat(STARM.out.log_final)
                .collect()
                //.view()


        MQC(mqc_ch1, ch_mqc_conf, ch_mqc_logo)
        COUNTSM(STARM.out.read_per_gene_tab.collect().flatten())

    }

    if( params.genome2 != null ){

        STARM2(unmapped_ch, genome_ch2, splitName_ch )

        mqc_ch2 = STARM2.out.read_per_gene_tab2
                        .concat(STARM2.out.log_final2)
                        .concat(mqc_ch1)
                        .collect()
                        //.view()

        // mqc_ch3 = mqc_ch1.concat(mqc_ch2).view()
    

        MQC2(mqc_ch2, ch_mqc_conf, ch_mqc_logo)
        COUNTSM2(STARM2.out.read_per_gene_tab2.collect().flatten())
    }

}


workflow {

    if ( params.mode == "SE" || params.mode == "SES" || params.mode == "SEBS" ){

        SINGLE()
    } 

    // else if ( params.mode == "UNMS" ){

    //     // SINGLE_SPLIT()
    // }

    else if ( params.mode == "PE" || params.mode == "PES" || params.mode == "PEBS" ){

        PAIRED()
    }
    
    // else if ( params.mode == "UNMP" ){

    //     // PAIRED_SPLIT()
    // }

    else {

        error "Invalid alignment mode provided"
        exit 0

    }

}



/* ---------------------------------------------------------------------------------------------------------
GLOBAL PROCESSES
------------------------------------------------------------------------------------------------------------ */


