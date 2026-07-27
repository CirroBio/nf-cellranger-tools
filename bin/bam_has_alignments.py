#!/usr/bin/env python3
# Exit 0 if the BAM file named in $1 holds at least one alignment record,
# 1 if it holds only a header, and 2 if it cannot be parsed.
# Samples for which cellranger called no cells produce a header-only BAM,
# which bamtofastq cannot convert.

import gzip
import struct
import sys

BAM_MAGIC = b"BAM\1"
MALFORMED = 2

fp = sys.argv[1]


def take(handle, n):
    buf = handle.read(n)
    if len(buf) != n:
        print(f"Truncated BAM: {fp}", file=sys.stderr)
        sys.exit(MALFORMED)
    return buf


def take_int32(handle):
    return struct.unpack("<i", take(handle, 4))[0]


try:
    # BAM is BGZF, which the gzip module reads as concatenated members
    with gzip.open(fp, "rb") as handle:
        if take(handle, 4) != BAM_MAGIC:
            print(f"Not a BAM file: {fp}", file=sys.stderr)
            sys.exit(MALFORMED)

        # Skip the SAM text header
        take(handle, take_int32(handle))

        # Skip the reference sequence dictionary
        for _ in range(take_int32(handle)):
            take(handle, take_int32(handle))  # name
            take(handle, 4)  # l_ref

        # Whatever remains is the start of the first alignment record
        has_alignments = bool(handle.read(1))

# Never report a BAM we could not read as one that simply holds no reads
except (OSError, EOFError, struct.error) as e:
    print(f"Could not read {fp}: {e}", file=sys.stderr)
    sys.exit(MALFORMED)

sys.exit(0 if has_alignments else 1)
