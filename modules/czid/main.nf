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
            "Sample Name"       : s,
            "Host Organism"     : "${czid_host}",
            "Sample Type"       : "${czid_sample_type}",
            "Collection Date"   : str(datetime.date.today()),
            "Collection Location": "not collected",
        })
EOF
    """
}


// ─────────────────────────────────────────────────────────────────────────────
// CZID_PROJECT  — create CZ ID project if it does not already exist
// ─────────────────────────────────────────────────────────────────────────────
process CZID_PROJECT {
    label 'process_low'
    errorStrategy 'ignore'

    input:
        val(czid_project)

    output:
        val(czid_project), emit: project_ready

    script:
    """
    python3 - <<'EOF'
import os, sys, requests

token   = os.environ.get("CZID_TOKEN", "")
project = "${czid_project}"
base    = "https://czid.org"
headers = {
    "Authorization" : f"Bearer {token}",
    "Content-Type"  : "application/json",
    "Accept"        : "application/json",
}

resp = requests.get(f"{base}/projects.json", headers=headers, timeout=30)
if resp.status_code != 200:
    print(f"Warning: could not list projects (status {resp.status_code}). Will attempt creation.")
    projects = []
else:
    projects = resp.json()

exists = any(p.get("name") == project for p in projects)

if exists:
    print(f"CZ ID project already exists: {project}")
else:
    r = requests.post(
        f"{base}/projects.json",
        headers=headers,
        json={"project": {"name": project, "public_access": 0}},
        timeout=30,
    )
    if r.status_code in (200, 201):
        print(f"Created CZ ID project: {project}")
    else:
        print(f"Warning: project creation returned status {r.status_code}: {r.text}")
EOF
    """
}


// ─────────────────────────────────────────────────────────────────────────────
// CZID_UPLOAD  — upload one sample's reads to CZ ID metagenomics
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
    czid login --token \$CZID_TOKEN
    czid metagenomics upload-sample \
        --project-name "${czid_project}" \
        --sample-name  "${id}"          \
        --input-file   ${r1}            \
        ${r2_flag}
    """
}
