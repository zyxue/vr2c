
CIGAR_TABLE = """# http://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment.cigartuples
BAM_CMATCH = 0                  # M
BAM_CINS = 1                    # I
BAM_CDEL = 2                    # D
BAM_CREF_SKIP = 3               # N
BAM_CSOFT_CLIP = 4              # S
BAM_CHARD_CLIP = 5              # H
BAM_CPAD = 6                    # P
BAM_CEQUAL = 7                  # =
BAM_CDIFF = 8                   # X
BAM_CBACK = 9                   # B
"""

exec(CIGAR_TABLE)
