#!/usr/bin/env python3
"""Review v1.2 decoded full-portrait frames against calculated radiance."""
from pathlib import Path
import sys

import review_kn_video as review
from render_kn_accretion import OUTPUT_DIR, source_hash


if __name__=='__main__':
    # The inherited compositor is identical; only the physical source/cache
    # and provenance function differ from the v1.1 review entry point.
    review.OUTPUT_DIR=OUTPUT_DIR
    review.source_hash=source_hash
    review.ENTRYPOINT=Path(__file__).resolve()
    review.main()
