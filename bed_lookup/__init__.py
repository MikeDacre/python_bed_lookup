"""
#============================================================================#
#                                                                            #
#          FILE: bed_lookup (python 3)                                       #
#        AUTHOR: Michael D Dacre, mike.dacre@gmail.com                       #
#  ORGANIZATION: Stanford University                                         #
#       LICENSE: MIT License, property of Stanford, use as you wish          #
#       CREATED: 2015-12-13 17:38                                            #
# Last modified: 2015-12-17 12:09                                            #
#                                                                            #
#============================================================================#
"""

# Define max length to use dictionary
_max_len = 180000000  # In bytes, this is approximately 5 million lines
# Line length dictionary:
#   1 million: ~ 35500000
#   5 million: ~ 180000000
#  10 million: ~ 366000000
# 100 million: ~ 3660000000
#
# NOTE: All calculations are based on snp files, bed files with gene models
# containing multiple exons will have different byte sizes. However, such files
# will tend to be much smaller than the above cutoff anyway, so it should not
# make much difference.

__all__ = ["_bed_lookup"]
from ._bed_lookup import BedFile
