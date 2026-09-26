"""
Copyright (c) 2026 Soham Bhole. All rights reserved.

Subsample a FASTA to N random records.

Usage: python subsample_fasta.py input.fasta output.fasta N [--seed 13]
"""

import argparse
import random


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("input")
    p.add_argument("output")
    p.add_argument("n", type=int)
    p.add_argument("--seed", type=int, default=13)
    args = p.parse_args()

    records, header, buf = [], None, []
    with open(args.input) as fh:
        for line in fh:
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(buf)))
                header, buf = line, []
            else:
                buf.append(line)
    if header is not None:
        records.append((header, "".join(buf)))

    n = min(args.n, len(records))
    chosen = random.Random(args.seed).sample(records, n)
    with open(args.output, "w") as fh:
        for h, s in chosen:
            fh.write(h)
            for i in range(0, len(s), 60):
                fh.write(s[i:i + 60] + "\n")
    print(f"Wrote {n} of {len(records)} records to {args.output}")


if __name__ == "__main__":
    main()
