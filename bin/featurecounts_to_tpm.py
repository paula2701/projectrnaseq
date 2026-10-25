#!/usr/bin/env python3
"""
Merge per-sample featureCounts tables into gene x sample tables:
  gene_counts.tsv  raw read/fragment counts
  gene_tpm.tsv     TPM = (count / length_kb) / sum(count / length_kb) * 1e6
Usage: featurecounts_to_tpm.py <sample1.featureCounts.tsv> [<sample2...>]
"""
# 
import os
import sys


def read_featurecounts(path):
    """Return {gene: length} and {gene: count} from one featureCounts output file."""
    lengths, counts = {}, {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("Geneid"):
                continue  # skip program comment line and header
            # columns: Geneid  Chr  Start  End  Strand  Length  <count>
            fields = line.rstrip("\n").split("\t")
            gene, length, count = fields[0], int(fields[5]), int(fields[-1])
            lengths[gene] = length
            counts[gene] = count
    return lengths, counts


def write_table(path, genes, samples, values, fmt):
    with open(path, "w") as out:
        out.write("gene_id\t" + "\t".join(samples) + "\n")
        for g in genes:
            out.write(g + "\t" + "\t".join(fmt.format(values[s][g]) for s in samples) + "\n")


def main(files):
    samples, lengths, counts = [], {}, {}
    for f in sorted(files):
        sample = os.path.basename(f).replace(".featureCounts.tsv", "")
        sample_lengths, sample_counts = read_featurecounts(f)
        samples.append(sample)
        lengths.update(sample_lengths)  # same GTF -> same lengths for every sample
        counts[sample] = sample_counts

    genes = sorted(lengths)
    tpm = {}
    for s in samples:
        rate = {g: counts[s][g] / (lengths[g] / 1000) for g in genes}  # reads per kilobase
        total = sum(rate.values())
        tpm[s] = {g: (rate[g] / total * 1e6 if total else 0.0) for g in genes}

    write_table("gene_counts.tsv", genes, samples, counts, "{}")
    write_table("gene_tpm.tsv", genes, samples, tpm, "{:.4f}")


if __name__ == "__main__":
    main(sys.argv[1:])