screen      = params.screen
runmode     = params.mode

process SCREENM {

    label 'process_screen'
    tag "$id"

    publishDir "fq_screen" , overwrite: true

    input:
        path(r1)
        path(screen_conf)

    output:
        path "*"
 


    script:
    

    if ( runmode == "SE" || runmode == "SES" || runmode == "SEBS" ){
     
    """
     fastq_screen --conf ${screen_conf} ${r1}       
    
    """
       
    }

    else if ( runmode == "PE" || runmode == "PES" || runmode == "PEBS" ){

    """
     fastq_screen --conf ${screen_conf} ${r1}       
    
    """

    }  else {

        error "Runmode ${runmode} is not supported"
        exit 0 
    }  

}