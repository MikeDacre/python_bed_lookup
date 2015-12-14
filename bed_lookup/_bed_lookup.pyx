"""
Lookup gene from a bed file
"""


class BedFile():
    """ A Lookup Object """
    def __init__(self, bedfile):
        self.dict = {}
        cdef int start, end
        with open(bedfile) as infile:
            for line in infile:
                f = line.rstrip().split('\t')
                if len(f) < 4:
                    continue
                chr   = f[0]
                start = int(f[1])
                end   = int(f[2])
                gene  = f[3]
                if chr not in self.dict:
                    self.dict[chr] = {}
                self.dict[chr][gene] = (start, end)

    def lookup(self, chromosome, location):
        """ Lookup your gene. Returns the gene name """
        cdef int l
        l = int(location)
        answer = ''
        for k, v in self.dict[chromosome].items():
            if v[0] < l < v[1]:
                answer = k
                break
        if answer:
            return answer
        else:
            return None
