process SALMON_MERGE {
    label 'process_single'

    conda "conda-forge::python=3.9"
    container "biocontainers/python:3.9--1"

    input:
    path quant_dirs // all per-sample salmon result folders

    output:
    path "gene_tpm.tsv"   , emit: tpm
    path "gene_counts.tsv", emit: counts
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    salmon_merge.py ${quant_dirs}
    """

    stub:
    """
    touch gene_tpm.tsv gene_counts.tsv
    """
}