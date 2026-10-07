process FEATURECOUNTS_TPM {
    label 'process_single'

    conda "conda-forge::python=3.12"
    container "docker.io/library/python:3.12-slim"

    input:
    path counts // all per-sample *.featureCounts.tsv files

    output:
    path "gene_tpm.tsv"   , emit: tpm
    path "gene_counts.tsv", emit: counts
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    featurecounts_to_tpm.py ${counts}
    """

    stub:
    """
    touch gene_tpm.tsv gene_counts.tsv
    """
}