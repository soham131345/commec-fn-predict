"""
Copyright (c) 2026 Soham Bhole. All rights reserved.

Demo generator for the Commec plugin.

Takes a toxin/virulence protein (from the test set, or any FASTA you pass),
reverse-translates it to DNA with randomized synonymous codons, optionally
applies amino-acid drift so it is distant from anything in a database, and
writes the FASTAs you need for the demo:

  demo/original.fasta    - natural-identity protein, codon-randomized DNA
  demo/recoded.fasta     - same protein, synonymously recoded (identical
                           protein, zero DNA-level similarity tricks needed)
  demo/drifted.fasta     - drifted protein (default 40% substitutions),
                           beyond reliable homology detection

Demo flow (needs commec installed, Linux/macOS or WSL):
  commec screen demo/drifted.fasta -d <commec_db_dir> -o demo/drifted
  python scripts/commec_fn_plugin.py --input-fasta demo/drifted.fasta \
      --commec-json demo/drifted.output.json --model-dir models/kmer \
      -o demo/drifted.augmented.json

Expected story: commec homology steps find nothing (PASS), the function-
prediction block in the augmented JSON says FLAG. Screenshot both.
"""

import argparse
import csv
import os
import random

CODONS = {
    "A": ["GCT", "GCC", "GCA", "GCG"], "R": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
    "N": ["AAT", "AAC"], "D": ["GAT", "GAC"], "C": ["TGT", "TGC"],
    "Q": ["CAA", "CAG"], "E": ["GAA", "GAG"], "G": ["GGT", "GGC", "GGA", "GGG"],
    "H": ["CAT", "CAC"], "I": ["ATT", "ATC", "ATA"], "L": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
    "K": ["AAA", "AAG"], "M": ["ATG"], "F": ["TTT", "TTC"],
    "P": ["CCT", "CCC", "CCA", "CCG"], "S": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
    "T": ["ACT", "ACC", "ACA", "ACG"], "W": ["TGG"], "Y": ["TAT", "TAC"],
    "V": ["GTT", "GTC", "GTA", "GTG"], "*": ["TAA"],
}
AA = "ACDEFGHIKLMNPQRSTVWY"


def read_fasta(path):
    seqs, header, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if header:
                    seqs[header] = "".join(buf)
                header, buf = line[1:].split()[0], []
            else:
                buf.append(line)
    if header:
        seqs[header] = "".join(buf)
    return seqs


def to_dna(protein, rng):
    return "".join(rng.choice(CODONS[a]) for a in protein if a in CODONS) + "TAA"


def drift(protein, rate, rng):
    out = []
    for c in protein:
        if c in AA and rng.random() < rate:
            out.append(rng.choice([a for a in AA if a != c]))
        else:
            out.append(c)
    return "".join(out)


def wrap(s, w=60):
    return "\n".join(s[i:i + w] for i in range(0, len(s), w))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--protein-fasta", help="FASTA containing the protein to use")
    src.add_argument("--from-test", nargs=2,
                     metavar=("TEST_FASTA", "TEST_LABELS"),
                     help="pick the first positive from the test split")
    p.add_argument("-o", "--outdir", required=True)
    p.add_argument("--drift-rate", type=float, default=0.4)
    p.add_argument("--seed", type=int, default=42)
    args = p.parse_args()

    rng = random.Random(args.seed)
    if args.protein_fasta:
        name, prot = next(iter(read_fasta(args.protein_fasta).items()))
    else:
        seqs = read_fasta(args.from_test[0])
        labels = {}
        with open(args.from_test[1]) as fh:
            for row in csv.DictReader(fh):
                labels[row["accession"]] = int(row["label"])
        name = next(a for a in seqs if labels.get(a) == 1)
        prot = seqs[name]

    os.makedirs(args.outdir, exist_ok=True)
    drifted = drift(prot, args.drift_rate, rng)
    identity = sum(a == b for a, b in zip(prot, drifted)) / len(prot) * 100

    files = {
        "original.fasta": ("original_" + name, to_dna(prot, rng)),
        "recoded.fasta": ("recoded_" + name, to_dna(prot, rng)),
        "drifted.fasta": (f"drifted{args.drift_rate}_" + name, to_dna(drifted, rng)),
    }
    for fname, (hdr, dna) in files.items():
        with open(os.path.join(args.outdir, fname), "w") as fh:
            fh.write(f">{hdr}\n{wrap(dna)}\n")

    with open(os.path.join(args.outdir, "proteins.fasta"), "w") as fh:
        fh.write(f">original_{name}\n{wrap(prot)}\n")
        fh.write(f">drifted_{name}\n{wrap(drifted)}\n")

    print(f"Protein: {name} ({len(prot)} aa)")
    print(f"Drifted protein identity to original: {identity:.1f}%")
    print(f"Wrote 3 DNA FASTAs + proteins.fasta to {args.outdir}/")
    print("Next: screen demo/drifted.fasta with commec, then run the plugin "
          "on the same file. Compare the verdicts.")


if __name__ == "__main__":
    main()
