#!/usr/bin/env python3
"""
Merge per-sample Salmon gene-level results (quant.genes.sf) into gene x sample tables:
  gene_counts.tsv  estimated reads/fragments per gene (Salmon NumReads, not whole numbers)
  gene_tpm.tsv     TPM as computed by Salmon
Usage: salmon_merge.py <sample1_dir> [<sample2_dir> ...]
"""
import os
import sys


def read_quant(path):
    """Return {gene: (tpm, numreads)} from one quant.genes.sf file."""
    values = {}
    with open(path) as fh:
        next(fh)  # header: Name  Length  EffectiveLength  TPM  NumReads
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            values[fields[0]] = (float(fields[3]), float(fields[4]))
    return values


def write_table(path, genes, samples, values, fmt):
    with open(path, "w") as out:
        out.write("gene_id\t" + "\t".join(samples) + "\n")
        for g in genes:
            out.write(g + "\t" + "\t".join(fmt.format(values[s][g]) for s in samples) + "\n")


def main(dirs):
    samples, data = [], {}
    for d in sorted(dirs):
        sample = os.path.basename(os.path.normpath(d))  # salmon result folder is named after the sample
        samples.append(sample)
        data[sample] = read_quant(os.path.join(d, "quant.genes.sf"))

    genes = sorted(set().union(*(data[s] for s in samples)))
    tpm = {s: {g: data[s].get(g, (0.0, 0.0))[0] for g in genes} for s in samples}
    counts = {s: {g: data[s].get(g, (0.0, 0.0))[1] for g in genes} for s in samples}

    write_table("gene_counts.tsv", genes, samples, counts, "{:.3f}")
    write_table("gene_tpm.tsv", genes, samples, tpm, "{:.4f}")


if __name__ == "__main__":
    main(sys.argv[1:])
