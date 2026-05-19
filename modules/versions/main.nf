process DUMP_VERSIONS {
    label 'process_low'
    publishDir "pipeline_info", mode: 'copy'

    input:
        path(versions, stageAs: "versions_?.yml")

    output:
        path "software_versions_mqc.yml", emit: mqc_yml
        path "software_versions.yml",     emit: yml

    script:
    """
    #!/usr/bin/env python3
    import yaml, glob

    # Collect all per-process versions
    all_versions = {}
    for f in glob.glob("versions_*.yml"):
        with open(f) as fh:
            data = yaml.safe_load(fh)
            if isinstance(data, dict):
                all_versions.update(data)

    # Write full process-level versions YAML for reference
    with open("software_versions.yml", "w") as fh:
        yaml.dump(all_versions, fh, default_flow_style=False)

    # Flatten and deduplicate by tool name (last-seen wins)
    tools = {}
    for process_versions in all_versions.values():
        if isinstance(process_versions, dict):
            for tool, ver in process_versions.items():
                tools[tool] = str(ver).strip()

    # Write MultiQC custom-content table (one row per tool)
    mqc = {
        "id": "pipeline_software_versions",
        "section_name": "Software Versions",
        "description": "Versions of software tools used in this pipeline run.",
        "section_href": "https://github.com/bixBeta/nextflow",
        "plot_type": "table",
        "pconfig": {
            "id": "pipeline_software_versions_table",
            "title": "Software Versions",
            "no_violin": True,
        },
        "headers": {
            "Version": {"description": "Software version used in this run", "scale": False, "format": "{}"}
        },
        "data": {tool: {"Version": ver} for tool, ver in sorted(tools.items())}
    }
    with open("software_versions_mqc.yml", "w") as fh:
        yaml.dump(mqc, fh, default_flow_style=False, allow_unicode=True)
    """
}
