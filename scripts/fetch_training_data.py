"""
Copyright (c) 2026 Soham Bhole. All rights reserved.
"""

"""
Fetch labelled protein training data for homology-independent toxin/virulence
function prediction. Uses only PUBLIC, already-annotated sequences (UniProt
Tox-Prot / keyword annotations, optionally a user-supplied VFDB FASTA).

Outputs:
  <outdir>/positives.fasta   - toxin / virulence-factor proteins
  <outdir>/negatives.fasta   - reviewed proteins excluding those keywords
  <outdir>/labels.csv        - accession,label(1|0),source

Requires: requests (pip install requests). Internet access to rest.uniprot.org.
"""

import argparse
import csv
import sys
import time

import requests

UNIPROT_STREAM = "https://rest.uniprot.org/uniprotkb/stream"

POSITIVE_QUERY = '((keyword:"Toxin") OR (keyword:"Virulence")) AND (reviewed:true)'

NEGATIVE_QUERY = (
    '(reviewed:true) AND (length:[50 TO 1500])'
    ' NOT (keyword:"Toxin") NOT (keyword:"Virulence")'
    ' NOT (keyword:"Enterotoxin") NOT (keyword:"Ion channel impairing toxin")'
)

def fetch_fasta(query: str, max_records: int, batch: int = 500) -> dict:
    """Stream FASTA from UniProt REST API; returns {accession: sequence}."""
    seqs = {}
    cursor = None
    while len(seqs) < max_records:
        size = min(batch, max_records - len(seqs))
        params = {"query": query, "format": "fasta", "size": size}
        if cursor:
            params["cursor"] = cursor
        for attempt in range(4):
            try:
                r = requests.get(UNIPROT_STREAM, params=params, timeout=120)
                if r.status_code == 200:
                    break
            except requests.RequestException:
                pass
            time.sleep(2 ** attempt)
        else:
            sys.exit("UniProt request failed repeatedly; check connectivity.")
        header, seq = None, []
        for line in r.text.splitlines():
            if line.startswith(">"):
                if header:
                    seqs.setdefault(header, "".join(seq))
                header = line[1:].split("|")[1] if "|" in line else line[1:].split()[0]
                seq = []
            else:
                seq.append(line.strip())
        if header:
            seqs.setdefault(header, "".join(seq))
        if len(seqs) >= max_records:
            break
        link = r.headers.get("Link", "")
        if "rel=\"next\"" in link:
            cursor = link.split("cursor=")[1].split(">")[0].rstrip('";')
        else:
            break
        print(f"  fetched {len(seqs)} so far...", file=sys.stderr)
    if len(seqs) > max_records:
        keep = set(list(seqs)[:max_records])
        seqs = {a: s for a, s in seqs.items() if a in keep}
    return seqs

def read_fasta(path: str) -> dict:
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

def write_fasta(seqs: dict, path: str):
    with open(path, "w") as fh:
        for acc, s in seqs.items():
            fh.write(f">{acc}\n")
            for i in range(0, len(s), 60):
                fh.write(s[i:i + 60] + "\n")

def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-o", "--outdir", required=True)
    p.add_argument("--max-positives", type=int, default=10000)
    p.add_argument("--max-negatives", type=int, default=20000)
    p.add_argument("--vfdb-fasta", default=None,
                   help="Optional local VFDB protein FASTA to merge into positives "
                        "(download VFDB_setA_pro.fas / VFDB_setB_pro.fas from mgc.ac.cn/VFs).")
    args = p.parse_args()

    import os
    os.makedirs(args.outdir, exist_ok=True)

    print("Fetching positives (Tox-Prot + virulence keywords)...", file=sys.stderr)
    pos = fetch_fasta(POSITIVE_QUERY, args.max_positives)
    sources = {a: "uniprot_toxprot_kw" for a in pos}
    if args.vfdb_fasta:
        vfdb = read_fasta(args.vfdb_fasta)
        for acc, s in vfdb.items():
            if acc not in pos:
                pos[acc] = s
                sources[acc] = "vfdb"
    print(f"  {len(pos)} positives.", file=sys.stderr)

    print("Fetching negatives (reviewed, toxin/virulence-excluded)...", file=sys.stderr)
    neg = fetch_fasta(NEGATIVE_QUERY, args.max_negatives)

    pos_seqs = set(pos.values())
    neg = {a: s for a, s in neg.items() if s not in pos_seqs}
    print(f"  {len(neg)} negatives.", file=sys.stderr)

    write_fasta(pos, f"{args.outdir}/positives.fasta")
    write_fasta(neg, f"{args.outdir}/negatives.fasta")
    with open(f"{args.outdir}/labels.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["accession", "label", "source"])
        for a in pos:
            w.writerow([a, 1, sources.get(a, "uniprot")])
        for a in neg:
            w.writerow([a, 0, "uniprot_swissprot"])
    print(f"Done. Wrote {len(pos)} positives, {len(neg)} negatives to {args.outdir}")

if __name__ == "__main__":
    main()
