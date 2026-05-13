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
// Renames files to <safe_id>_R1.fastq.gz format before upload:
//   - underscores in sample label replaced with hyphens (CZ ID fails on
//     multiple underscores in the filename)
//   - suffix normalised to .fastq.gz (_R1/_R2 required by CZ ID)
// ─────────────────────────────────────────────────────────────────────────────
process CZID_UPLOAD {
    tag "$id"
    label 'process_low'
    errorStrategy 'ignore'
    secret 'CZID_CLI_SECRET'

    input:
        tuple val(id), path(reads)
        val(czid_project)

    script:
    def r1 = reads instanceof List ? reads[0] : reads
    def r2 = reads instanceof List && reads.size() > 1 ? reads[1] : null
    """
    export CZID_CLI_ACCEPTED_USER_AGREEMENT=Y

    # Redirect HOME to work dir so czid-cli doesn't try to write
    # to ~/.config/czid-cli/ which is read-only in Singularity
    export HOME=\$(pwd)
    mkdir -p \${HOME}/.config/czid-cli
    printf 'secret: %s\\naccepted_user_agreement: Y\\n' "\${CZID_CLI_SECRET}" > \${HOME}/.config/czid-cli/config.yaml

    # Replace underscores in sample label with hyphens so the final
    # filename has exactly one underscore: <safe_id>_R1.fastq.gz
    SAFE_ID=\$(echo "${id}" | tr '_' '-')

    ln -sf ${r1} \${SAFE_ID}_R1.fastq.gz
    """ + (r2 ? "ln -sf ${r2} \${SAFE_ID}_R2.fastq.gz\n" : "") + """
    czid metagenomics upload-sample \
        \${SAFE_ID}_R1.fastq.gz \
        """ + (r2 ? "\${SAFE_ID}_R2.fastq.gz \\\n        " : "") + """--project  "${czid_project}" \
        --sample-name "\${SAFE_ID}" \
        --sequencing-platform Illumina
    """
}
