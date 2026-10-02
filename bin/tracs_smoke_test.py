#!/usr/bin/env python
"""Call each TRACS C++ kernel on tiny inputs.

TRACS compiles its extension with -O3 -march=<target>, so a build for a newer CPU
dies with SIGILL (exit 132) on import or inside these kernels. Used by
TRACS_PREFLIGHT.
"""
import os
import tempfile

import numpy as np
from TRACS import calculate_posteriors, pairsnp, trans_dist

with tempfile.TemporaryDirectory() as tmp:
    msa = os.path.join(tmp, "msa.fasta")
    with open(msa, "w") as fh:
        fh.write(">a\nACGTACGTACGTACGTACGT\n>b\nACGTACGTACGAACGTACGT\n>c\nACGTTCGTACGTACGTACGA\n")
    pairsnp(fasta=[msa], n_threads=1, dist=1000000, filter=False)

calculate_posteriors(np.random.default_rng(0).random((1000, 4)) * 10, [1.0, 0.5, 0.1, 0.05], False, 0.01)
trans_dist([0, 1, 5], [0.0, 0.0, 0.0], 10.0, 0.1, 1e-6)
print("TRACS C++ extension OK on this CPU")
