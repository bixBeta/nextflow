process DUMP_VERSIONS {
    label 'process_low'
    publishDir "pipeline_info", mode: 'copy'

    input:
        path versions

    output:
        path "software_versions_mqc.yml", emit: mqc_yml
        path "software_versions.yml",     emit: yml

    script:
    """
    #!/usr/bin/env python3
    import yaml, glob

    versions = {}
    for f in glob.glob("*.yml"):
        if f.startswith("software_versions"):
            continue
        with open(f) as fh:
            data = yaml.safe_load(fh)
            if isinstance(data, dict):
                versions.update(data)

    # Write plain versions YAML
    with open("software_versions.yml", "w") as fh:
        yaml.dump(versions, fh, default_flow_style=False)

    # Write MultiQC-compatible custom data YAML
    rows = "".join(
        f'        <dt>{tool}</dt><dd><samp>{ver}</samp></dd>\\n'
        for process_versions in versions.values()
        for tool, ver in (process_versions.items() if isinstance(process_versions, dict) else {}.items())
    )
    mqc = {
        "id": "software_versions",
        "section_name": "Software Versions",
        "section_href": "https://github.com/bixBeta/nextflow",
        "plot_type": "html",
        "description": "Versions collected at run time.",
        "data": f"<dl class=\\"dl-horizontal\\">\\n{rows}    </dl>"
    }
    with open("software_versions_mqc.yml", "w") as fh:
        yaml.dump(mqc, fh, default_flow_style=False)
    """
}
