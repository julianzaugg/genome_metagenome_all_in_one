/*
 * Collect all versions.yml files emitted during the run and merge them into
 * a single software_versions.tsv published to pipeline_info/.
 */

process DUMP_SOFTWARE_VERSIONS {
    label 'process_single'

    input:
    path(versions, stageAs: '?_*')

    output:
    path 'software_versions.tsv', emit: tsv
    path 'versions.yml',          emit: versions

    script:
    """
    merge_versions.py $versions > software_versions.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    echo -e "tool\tversion\tprocess" > software_versions.tsv
    echo -e "example_tool\t1.0.0\tSTUB:DUMMY" >> software_versions.tsv
    echo '"${task.process}": {python: stub}' > versions.yml
    """
}
