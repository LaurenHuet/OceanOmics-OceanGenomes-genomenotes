#!/usr/bin/env python3

import sys

INPUT  = sys.argv[1]
OUTPUT = sys.argv[2]

EXPECTED_COLS = [
    "qseqid", "staxids", "bitscore", "qseqid2", "sseqid", "pident", "length",
    "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue"
]

def parse_seqkit_sliding(qseqid):
    """
    Parse seqkit sliding ID:
    SUPER_1_sliding:6000001-6100000
    -> (SUPER_1, 6000001)
    """
    if "_sliding:" not in qseqid:
        return qseqid, None

    base, coords = qseqid.split("_sliding:")
    start, _ = coords.split("-")
    return base, int(start)

with open(INPUT) as fin, open(OUTPUT, "w") as fout:
    # write BlobTools-required header
    fout.write("\t".join(EXPECTED_COLS) + "\n")

    for line in fin:
        line = line.rstrip("\n")

        # skip header if present
        if line.startswith("qseqid\t"):
            continue

        fields = line.split("\t")
        if len(fields) != 14:
            raise ValueError(f"Expected 14 columns, got {len(fields)}:\n{line}")

        qseqid = fields[0]
        qstart = int(fields[9])
        qend   = int(fields[10])

        scaffold, window_start = parse_seqkit_sliding(qseqid)

        if window_start is not None:
            # lift coordinates from window-local (1-based) to genome
            fields[9]  = str(window_start + qstart - 1)
            fields[10] = str(window_start + qend   - 1)

        # rewrite seq IDs to scaffold name
        fields[0] = scaffold
        fields[3] = scaffold

        fout.write("\t".join(fields) + "\n")
