/*
 * Group-scoped isolate comparative analyses.
 */

process PREP_COMPARISON_GROUPS {
    label 'process_low'

    input:
    val assembly_records
    val gff_records
    val faa_records
    val comparison_manifest
    val samples_include

    output:
    path 'comparison_groups/*', emit: groups
    path 'versions.yml',        emit: versions

    script:
    def assemblies = assembly_records.collect { meta, fasta -> "${meta.id}\t${meta.group}\t${fasta}" }.join('\n')
    def gffs       = gff_records.collect      { meta, gff   -> "${meta.id}\t${gff}" }.join('\n')
    def faas       = faa_records.collect      { meta, faa   -> "${meta.id}\t${faa}" }.join('\n')
    def assemblies_json = groovy.json.JsonOutput.toJson(assemblies)
    def gffs_json       = groovy.json.JsonOutput.toJson(gffs)
    def faas_json       = groovy.json.JsonOutput.toJson(faas)
    """
    mkdir -p comparison_groups

    python - <<'PY'
    import csv
    import os
    import re
    import shutil
    import zlib
    from pathlib import Path

    Path("sample_assemblies.tsv").write_text(${assemblies_json})
    Path("sample_gffs.tsv").write_text(${gffs_json})
    Path("sample_faas.tsv").write_text(${faas_json})

    project_dir = Path("${projectDir}").resolve()
    comparison_manifest = "${comparison_manifest}"
    samples_include = "${samples_include}"

    def truthy(value):
        return str(value or "").strip().lower() in {"1", "true", "t", "yes", "y"}

    def read_tsv(path, ncols):
        rows = []
        p = Path(path)
        if not p.exists():
            return rows
        for line in p.read_text().splitlines():
            if not line.strip():
                continue
            parts = line.rstrip("\\n").split("\\t")
            while len(parts) < ncols:
                parts.append("")
            rows.append(parts[:ncols])
        return rows

    def link_or_copy(src, dest):
        src = Path(src).resolve()
        dest = Path(dest)
        dest.parent.mkdir(parents=True, exist_ok=True)
        if dest.exists() or dest.is_symlink():
            dest.unlink()
        try:
            os.symlink(src, dest)
        except OSError:
            shutil.copy2(src, dest)

    def clean_id(value):
        cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", value)
        return cleaned or "entry"

    def resolve_input(value):
        if not value:
            return ""
        p = Path(value)
        if not p.is_absolute():
            p = project_dir / p
        return str(p)

    sample_rows = read_tsv("sample_assemblies.tsv", 3)
    gffs = {sid: path for sid, path in read_tsv("sample_gffs.tsv", 2)}
    faas = {sid: path for sid, path in read_tsv("sample_faas.tsv", 2)}
    samples = {
        sid: {
            "comparison_group": group or "all",
            "entry_type": "sample",
            "id": sid,
            "fasta": fasta,
            "gff": gffs.get(sid, ""),
            "faa": faas.get(sid, ""),
            "parsnp_reference": False,
        }
        for sid, group, fasta in sample_rows
    }

    include_ids = None
    if samples_include and samples_include != "null":
        include_path = Path(resolve_input(samples_include))
        if include_path.exists():
            include_ids = {
                line.strip().split()[0]
                for line in include_path.read_text().splitlines()
                if line.strip() and not line.startswith("#")
            }

    groups = {}
    manifest_path = Path(resolve_input(comparison_manifest)) if comparison_manifest and comparison_manifest != "null" else None
    if manifest_path and manifest_path.exists():
        with manifest_path.open(newline="") as handle:
            reader = csv.DictReader(handle)
            required = {"comparison_group", "entry_type", "id"}
            missing = required - set(reader.fieldnames or [])
            if missing:
                raise SystemExit(f"comparison_manifest missing columns: {', '.join(sorted(missing))}")
            for row in reader:
                group = row["comparison_group"].strip()
                entry_type = row["entry_type"].strip().lower()
                entry_id = row["id"].strip()
                if not group or not entry_id:
                    raise SystemExit("comparison_manifest rows require comparison_group and id")
                if entry_type == "sample":
                    if entry_id not in samples:
                        raise SystemExit(f"comparison_manifest sample id not found in samplesheet: {entry_id}")
                    entry = dict(samples[entry_id])
                    entry["comparison_group"] = group
                elif entry_type == "reference":
                    fasta = resolve_input((row.get("fasta") or "").strip())
                    if not fasta:
                        raise SystemExit(f"reference row {entry_id} requires fasta")
                    entry = {
                        "comparison_group": group,
                        "entry_type": "reference",
                        "id": entry_id,
                        "fasta": fasta,
                        "gff": resolve_input((row.get("gff") or "").strip()),
                        "faa": resolve_input((row.get("faa") or "").strip()),
                        "parsnp_reference": truthy(row.get("parsnp_reference")),
                    }
                else:
                    raise SystemExit(f"entry_type must be sample or reference: {entry_type}")
                groups.setdefault(group, []).append(entry)
    else:
        for sid, entry in samples.items():
            if include_ids is not None and sid not in include_ids:
                continue
            groups.setdefault(entry["comparison_group"], []).append(entry)

    if not groups:
        raise SystemExit("No isolate comparison groups were created")

    for group, entries in groups.items():
        n_refs = sum(1 for e in entries if e.get("parsnp_reference"))
        if n_refs > 1:
            raise SystemExit(f"comparison group {group} has more than one parsnp_reference=true row")

        group_id = clean_id(group)
        out = Path("comparison_groups") / group_id
        fastas = out / "fastas"
        gffs_dir = out / "gffs"
        faas_dir = out / "faas"
        chewie = out / "chewbbaca_fastas"
        for d in (fastas, gffs_dir, faas_dir, chewie):
            d.mkdir(parents=True, exist_ok=True)

        selected_reference = next((e for e in entries if e.get("parsnp_reference")), None)
        if selected_reference is None:
            selected_reference = next((e for e in entries if e["entry_type"] == "reference"), entries[0])

        with (out / "entries.tsv").open("w") as handle:
            handle.write("comparison_group\\tentry_type\\tid\\tfasta\\tgff\\tfaa\\tparsnp_reference\\n")
            for e in entries:
                eid = clean_id(e["id"])
                fasta_dest = fastas / f"{eid}.fasta"
                link_or_copy(e["fasta"], fasta_dest)
                gff_dest = ""
                faa_dest = ""
                if e.get("gff"):
                    gff_dest = gffs_dir / f"{eid}.gff3"
                    link_or_copy(e["gff"], gff_dest)
                if e.get("faa"):
                    faa_dest = faas_dir / f"{eid}.faa"
                    link_or_copy(e["faa"], faa_dest)

                hashed = eid
                if len(hashed) > 37:
                    hashed = f"g{zlib.crc32(e['id'].encode()):010d}"
                link_or_copy(fasta_dest, chewie / f"{hashed}.fasta")

                handle.write("\\t".join([
                    group_id,
                    e["entry_type"],
                    e["id"],
                    str(fasta_dest),
                    str(gff_dest),
                    str(faa_dest),
                    "true" if e is selected_reference else "false",
                ]) + "\\n")

        with (out / "parsnp_reference.txt").open("w") as handle:
            handle.write(str((fastas / f"{clean_id(selected_reference['id'])}.fasta").resolve()) + "\\n")

        with (out / "genome_hash_map.tsv").open("w") as handle:
            handle.write("original_id\\tchewbbaca_id\\n")
            for e in entries:
                eid = clean_id(e["id"])
                hashed = eid if len(eid) <= 37 else f"g{zlib.crc32(e['id'].encode()):010d}"
                handle.write(f"{e['id']}\\t{hashed}\\n")
    PY

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """.stripIndent()

}

process FASTANI {
    tag   { meta.id }
    label 'process_high'

    input:
    tuple val(meta), path(group_dir)

    output:
    tuple val(meta), path("${meta.id}.fastani.tsv"), emit: results
    tuple val(meta), path("${meta.id}.fastani.matrix"), emit: matrix, optional: true
    path 'versions.yml', emit: versions

    script:
    def args = task.ext.args ?: '--matrix'
    """
    find ${group_dir}/fastas -type f -name '*.fasta' | sort > reference_list.txt
    fastANI ${args} \\
        --queryList reference_list.txt \\
        --refList reference_list.txt \\
        -o ${meta.id}.fastani.tsv \\
        -t ${task.cpus}
    if [[ -f ${meta.id}.fastani.tsv.matrix ]]; then
        mv ${meta.id}.fastani.tsv.matrix ${meta.id}.fastani.matrix
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastani: \$(fastANI --version 2>&1 | head -1)
    END_VERSIONS
    """

    stub:
    """
    echo -e "query\\tref\\tani\\tfragments_mapped\\tfragments_total" > ${meta.id}.fastani.tsv
    touch ${meta.id}.fastani.matrix
    echo '"${task.process}": {fastani: stub}' > versions.yml
    """
}

process PARSNP {
    tag   { meta.id }
    label 'process_high'
    label 'process_long'

    input:
    tuple val(meta), path(group_dir)

    output:
    tuple val(meta), path("${meta.id}_parsnp/parsnp.aln"), emit: alignment
    tuple val(meta), path("${meta.id}_parsnp"),            emit: outdir
    path 'versions.yml',                                   emit: versions

    script:
    def args = task.ext.args ?: '-c'
    """
    outdir=${meta.id}_parsnp
    ref=\$(cat ${group_dir}/parsnp_reference.txt)
    parsnp \\
        -r "\$ref" \\
        -d ${group_dir}/fastas/*.fasta \\
        ${args} \\
        -p ${task.cpus} \\
        -o \$outdir
    harvesttools -i \$outdir/parsnp.ggr -M \$outdir/parsnp.aln || touch \$outdir/parsnp.aln

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        parsnp: \$(parsnp --version 2>&1 | head -1)
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${meta.id}_parsnp
    echo ">stub" > ${meta.id}_parsnp/parsnp.aln
    echo "ACGT" >> ${meta.id}_parsnp/parsnp.aln
    echo '"${task.process}": {parsnp: stub}' > versions.yml
    """
}

process GUBBINS {
    tag   { meta.id }
    label 'process_high'
    label 'process_long'

    input:
    tuple val(meta), path(alignment)

    output:
    tuple val(meta), path("${meta.id}.filtered_polymorphic_sites.fasta"), emit: polymorphic_sites, optional: true
    tuple val(meta), path("${meta.id}.*"), emit: outputs
    path 'versions.yml', emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    run_gubbins.py ${args} --threads ${task.cpus} -p ${meta.id} ${alignment}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gubbins: \$(run_gubbins.py --version 2>&1 | head -1)
    END_VERSIONS
    """

    stub:
    """
    echo ">stub" > ${meta.id}.filtered_polymorphic_sites.fasta
    echo "ACGT" >> ${meta.id}.filtered_polymorphic_sites.fasta
    echo '"${task.process}": {gubbins: stub}' > versions.yml
    """
}

process PANAROO_RUN {
    tag   { meta.id }
    label 'process_high'
    label 'process_long'

    input:
    tuple val(meta), path(group_dir)

    output:
    tuple val(meta), path("${meta.id}_panaroo/core_gene_alignment.aln"), emit: core_alignment, optional: true
    tuple val(meta), path("${meta.id}_panaroo/gene_presence_absence.csv"), emit: gene_presence_absence, optional: true
    tuple val(meta), path("${meta.id}_panaroo"), emit: outdir
    path 'versions.yml', emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    panaroo \\
        -i ${group_dir}/gffs/*.gff3 \\
        --out_dir ${meta.id}_panaroo \\
        -t ${task.cpus} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        panaroo: \$(panaroo --version 2>&1 | head -1)
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${meta.id}_panaroo
    echo ">stub" > ${meta.id}_panaroo/core_gene_alignment.aln
    echo "ACGT" >> ${meta.id}_panaroo/core_gene_alignment.aln
    echo "Gene,stub" > ${meta.id}_panaroo/gene_presence_absence.csv
    echo '"${task.process}": {panaroo: stub}' > versions.yml
    """
}

process CHEWBACCA_RUN {
    tag   { meta.id }
    label 'process_high'
    label 'process_long'

    input:
    tuple val(meta), path(group_dir)
    val training_file
    val thresholds

    output:
    tuple val(meta), path("${meta.id}_chewbbaca"), emit: outdir
    tuple val(meta), path("${meta.id}_chewbbaca/genome_hash_map.tsv"), emit: hash_map
    path 'versions.yml', emit: versions

    script:
    def training = training_file ? "--training-file ${training_file}" : ''
    def threshold_args = thresholds ?: '0.90 0.95 0.99 1'
    """
    mkdir -p ${meta.id}_chewbbaca
    cp ${group_dir}/genome_hash_map.tsv ${meta.id}_chewbbaca/genome_hash_map.tsv
    chewie CreateSchema \\
        --input-files ${group_dir}/chewbbaca_fastas \\
        --output-directory ${meta.id}_chewbbaca/schema \\
        --cpu-cores ${task.cpus} \\
        ${training}
    chewie AlleleCall \\
        --input-files ${group_dir}/chewbbaca_fastas \\
        --schema-directory ${meta.id}_chewbbaca/schema/schema_seed \\
        --output-directory ${meta.id}_chewbbaca/alleles \\
        --cpu-cores ${task.cpus} \\
        ${training}
    chewBBACA.py RemoveGenes \\
        -i ${meta.id}_chewbbaca/alleles/results_alleles.tsv \\
        -g ${meta.id}_chewbbaca/alleles/paralogous_counts.tsv \\
        -o ${meta.id}_chewbbaca/alleles/results_alleles_NoParalogs.tsv
    chewie ExtractCgMLST \\
        --input-file ${meta.id}_chewbbaca/alleles/results_alleles_NoParalogs.tsv \\
        --output-directory ${meta.id}_chewbbaca/cgMLST \\
        --threshold ${threshold_args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        chewbbaca: \$(chewie --version 2>&1 | head -1)
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${meta.id}_chewbbaca/cgMLST
    cp ${group_dir}/genome_hash_map.tsv ${meta.id}_chewbbaca/genome_hash_map.tsv
    echo -e "FILE\\tAllele1" > ${meta.id}_chewbbaca/cgMLST/cgMLST.tsv
    echo '"${task.process}": {chewbbaca: stub}' > versions.yml
    """
}

process IQTREE {
    tag   { meta.id }
    label 'process_high'

    input:
    tuple val(meta), path(alignment)

    output:
    tuple val(meta), path("${meta.id}.treefile"), emit: tree, optional: true
    tuple val(meta), path("${meta.id}.*"),        emit: outputs
    path 'versions.yml',                          emit: versions

    script:
    def args = task.ext.args ?: '-fast -m GTR'
    """
    iqtree -s ${alignment} -pre ${meta.id} -nt ${task.cpus} ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iqtree: \$(iqtree --version 2>&1 | head -1)
    END_VERSIONS
    """

    stub:
    """
    echo "(stub);" > ${meta.id}.treefile
    echo '"${task.process}": {iqtree: stub}' > versions.yml
    """
}
