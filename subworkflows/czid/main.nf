include { CZID_META    } from '../../modules/czid/main.nf'
include { CZID_PROJECT } from '../../modules/czid/main.nf'
include { CZID_UPLOAD  } from '../../modules/czid/main.nf'

workflow CZID {

    take:
        sheet            // path  : sample-sheet.csv
        reads_ch         // channel: tuple val(id), path(reads) — trimmed OR unmapped
        czid_project     // val
        czid_host        // val
        czid_sample_type // val

    main:
        CZID_META(sheet, czid_host, czid_sample_type)
        CZID_PROJECT(czid_project)
        CZID_UPLOAD(reads_ch, CZID_PROJECT.out.project_ready)

}
