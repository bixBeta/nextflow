include { CZID_META   } from '../../modules/czid/main.nf'
include { CZID_UPLOAD } from '../../modules/czid/main.nf'

workflow CZID {

    take:
        sheet            // path  : sample-sheet.csv
        reads_ch         // channel: tuple val(id), path(reads) — trimmed OR unmapped
        czid_project     // val
        czid_host        // val
        czid_sample_type // val

    main:
        CZID_META(sheet, czid_host, czid_sample_type)

        reads_ch
            .view { id, reads ->
                def r1 = reads instanceof List ? reads[0] : reads
                def r2 = reads instanceof List && reads.size() > 1 ? reads[1] : null
                def safe_id = id.replaceAll('_', '-')
                "CZID_INPUT >> sample: ${id}\n  safe_id: ${safe_id}\n  R1: ${r1}\n  R2: ${r2 ?: 'N/A'}\n  -> ${safe_id}_R1.fastq.gz / ${r2 ? safe_id + '_R2.fastq.gz' : 'N/A'}"
            }

        CZID_UPLOAD(reads_ch, czid_project)

}
