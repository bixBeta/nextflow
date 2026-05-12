// ─────────────────────────────────────────────────────────────────────────────
// CZID_META  — auto-generate a minimal CZ ID metadata CSV from sample-sheet
// ─────────────────────────────────────────────────────────────────────────────
process CZID_META {
    label 'process_low'

    input:
        path(sheet)
        val(czid_host)
        val(czid_sample_type)

    output:
        path "czid_metadata.csv", emit: meta_csv

    script:
    """
    python3 - <<'EOF'
import csv, datetime

with open("${sheet}") as f:
    samples = [row["label"] for row in csv.DictReader(f)]

with open("czid_metadata.csv", "w", newline="") as f:
    fields = ["Sample Name", "Host Organism", "Sample Type",
              "Collection Date", "Collection Location"]
    w = csv.DictWriter(f, fieldnames=fields)
    w.writeheader()
    for s in samples:
        w.writerow({
            "Sample Name"        : s,
            "Host Organism"      : "${czid_host}",
            "Sample Type"        : "${czid_sample_type}",
            "Collection Date"    : str(datetime.date.today()),
            "Collection Location": "not collected",
        })
EOF
    """
}


// ─────────────────────────────────────────────────────────────────────────────
// CZID_UPLOAD  — upload one sample's reads to CZ ID metagenomics
// Auth: CZID_CLI_SECRET env var (Tower secret) — czid-cli reads it automatically
// ─────────────────────────────────────────────────────────────────────────────
process CZID_UPLOAD {
    tag "$id"
    label 'process_low'
    errorStrategy 'ignore'

    input:
        tuple val(id), path(reads)
        val(czid_project)

    script:
    def r1      = reads instanceof List ? reads[0] : reads
    def r2_flag = reads instanceof List && reads.size() > 1 ? "--input-file-2 ${reads[1]}" : ""
    """
    czid metagenomics upload-sample \
        --project-name "${czid_project}" \
        --sample-name  "${id}"          \
        --input-file   ${r1}            \
        ${r2_flag}
    """
}
