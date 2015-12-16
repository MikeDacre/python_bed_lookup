"""
#============================================================================#
#                                                                            #
#          FILE: bed_lookup (python 3)                                       #
#        AUTHOR: Michael D Dacre, mike.dacre@gmail.com                       #
#  ORGANIZATION: Stanford University                                         #
#       LICENSE: MIT License, property of Stanford, use as you wish          #
#       CREATED: 2015-12-13 17:38                                            #
# Last modified: 2015-12-16 14:51                                            #
#                                                                            #
#============================================================================#
"""

# Define max length to use dictionary
_max_len = 1500000  # 1,000,000 is a reasonable default

__all__ = ["_bed_lookup"]
from ._bed_lookup import BedFile
