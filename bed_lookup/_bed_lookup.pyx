"""
Lookup gene from a bed file with either
a cythonized dictionary or sqlite depending
on the size of the bedfile
Lookup method is transparent to the user
"""
import sys
import os
import sqlite3
import gzip
import bz2
from subprocess import check_output as sub
from collections import defaultdict
from os.path import getsize

# C++ Library Import
from libcpp.string cimport string

# Get max len from __init__.py
from . import _max_len

from . import logme
logme.MIN_LEVEL = 'info'


def gopen(infile, mode='r'):
    """ Return file handle of file regardless of zipped or not
        Text mode enforced for compatibility with python2 """
    mode   = mode[0] + 't'
    p2mode = mode
    if hasattr(infile, 'write'):
        return infile
    if isinstance(infile, str):
        if infile.endswith('.gz'):
            return gzip.open(infile, mode)
        if infile.endswith('.bz2'):
            if hasattr(bz2, 'open'):
                return bz2.open(infile, mode)
            else:
                return bz2.BZ2File(infile, p2mode)
        return open(infile, p2mode)


cdef large_file(infile):
    """ Check if file is greater than _max_len """
    size = getsize(infile)*3.7 if infile.endswith('.gz') else getsize(infile)
    return size > _max_len


cdef class Gene(object):
    cdef long start, end
    cdef public string name

    def __init__(self, values):
        self.start = int(values[1])
        self.end = int(values[2])
        self.name = values[3].encode()

    def find(self, i):
        return self.start <= i < self.end

    def __repr__(self):
        return "{}({}:{})".format(self.name, self.start, self.end)


cdef class Chrom(list):
    def add(self, aGene):
        self.append(aGene)

    def find(self, loc):
        cdef long i = int(loc)
        for gene in self:
            if gene.find(i):
                # Handle python 2/3 strings
                return gene.name.decode()
        return None

    def __repr__(self):
        astr = []
        for gene in self:
            astr += [repr(gene)]
        return '; '.join(astr)


class BedFile():
    """ An object to hold data from a bed file and allow lookup by coordinate.
        To use, create a bedfile object with your bedfile:
            b = BedFile(bedfile)
        Then search by coordinate:
            result = b.lookup('chr1', 10003021) """

    def lookup(self, chromosome, location):
        """ Lookup your gene. Returns the gene name """
        if self._type == 'sq':
            return self._lookup_sqlite(chromosome, str(location))
        elif self._type == 'dt':
            return self._lookup_dict(chromosome, location)

    def lookup_df(self, df, chrom_col, pos_col):
        """Use a pandas dataframe and return a series with the same index.

        Args:
            df (DataFrame):  A pandas dataframe
            chrom_col (str): The name of the column with the chromosome name
            pos_col (str):   The name of the column with the position

        Returns:
            Series: A pandas series with the same index as the original df.
        """
        return df.apply(self.lookup_series, args=(chrom_col, pos_col), axis=1)

    def lookup_series(self, series, chrom_col, pos_col):
        """To use with pandas DataFrame.apply().

        Note: you need to specify `axis=1` when using with df.apply().

        Args:
            series (Series): A pandas series from df.apply
            chrom_col (str): The name of the column with the chromosome name
            pos_col (str):   The name of the column with the position

        Returns:
            Series: A pandas series with the same index as the original df.
        """
        return self.lookup(series[chrom_col], series[pos_col])

    # Private functions
    def _lookup_sqlite(self, chromosome, location):
        """ Simple sqlite query """
        expr = ("SELECT name FROM '{0}' INDEXED BY '{0}_start_end' " +
                "WHERE {1} BETWEEN start AND end").format(chromosome, location)

        try:
            self._c.execute(expr)
        except sqlite3.OperationalError as e:
            if str(e).startswith('no such table'):
                logme.log(("Chromosome '{}' is not in " +
                           "the lookup table, lookup failed." +
                           "\n").format(chromosome), level='error')
                return None
            else:
                raise(e)

        answer = self._c.fetchone()
        if answer:
            return answer[0]
        else:
            logme.log(("Location '{}' on Chromosome '{}' " +
                       "is not in the lookup table, lookup failed." +
                       "\n").format(location, chromosome), level='debug')
            return None

    def _lookup_dict(self, chromosome, location):
        """ Simple dictionary query with cython for math """
        cdef int loc
        loc = int(location)
        if chromosome in self._data:
            ans = self._data[chromosome].find(location)
            if ans:
                return ans
            else:
                logme.log(("Location '{}' on Chromosome '{}' " +
                           "is not in the lookup table, lookup failed." +
                           "\n").format(location, chromosome), level='debug')
                return None
        else:
            logme.log(("Chromosome '{}' is not in " +
                       "the lookup table, lookup failed." +
                       "\n").format(chromosome), level='error')
            return None

    def _init_sqlite(self, bedfile):
        """ Initialize sqlite3 object """
        logme.log('Bedfile is large, using sqlite\n', level='info')
        db_name = bedfile if bedfile.endswith('.db') else bedfile + '.db'
        # Check if the alternate db exists if db doesn't exist
        if not os.path.exists(db_name):
            if bedfile.endswith('.gz'):
                alt_path = '.'.join(bedfile.split('.')[:-1]) + '.db'
            else:
                alt_path = bedfile + '.gz' + '.db'
            if os.path.exists(alt_path):
                db_name = alt_path
                exists = True
            else:
                exits = False
        else:
            exists = True

        # If the database already exists, use it
        if exists:
            logme.log('Using existing db, if this is ' +
                      'not what you want, delete ' + db_name + '\n',
                      level='info')
            self._conn = sqlite3.connect(db_name)
            self._c = self._conn.cursor()
            return

        # Create an sqlite database from bed file
        logme.log('Creating sqlite database, this ' +
                  'may take a long time.\n', level='info')
        self._conn = sqlite3.connect(db_name)
        self._c = self._conn.cursor()

        with gopen(bedfile) as infile:
            for line in infile:
                f = line.rstrip().split('\t')
                if len(f) < 4:
                    continue
                # Check if db exists and create if it does
                expr = ("SELECT * FROM sqlite_master WHERE name = '{}' " +
                        "and type='table';").format(f[0])
                self._c.execute(expr)
                if not self._c.fetchall():
                    exp = ("CREATE TABLE '{}' (name text, start int, " +
                           "end int);").format(f[0])
                    self._c.execute(exp)
                    self._conn.commit()
                expr = ("INSERT INTO '{}' VALUES " +
                        "('{}','{}','{}')").format(f[0], f[3], f[1], f[2])
                self._c.execute(expr)
            self._conn.commit()
            # Create indicies
            self._c.execute('''SELECT name FROM sqlite_master WHERE type='table';''')
            for i in self._c.fetchall():
                exp = ("CREATE INDEX '{0}_start_end' ON '{0}' " +
                       "(start, end)").format(i[0])
                self._c.execute(exp)
                self._conn.commit()

    def _init_dict(self, bedfile):
        self._data = defaultdict(Chrom)
        with gopen(bedfile) as infile:
            for line in infile:
                f = line.rstrip().split('\t')
                if len(f) < 4:
                    continue
                chr   = f[0]
                gene  = Gene(f)
                self._data[chr].add(gene)

    def __init__(self, bedfile):
        if bedfile.endswith('.db') or large_file(bedfile):
            # Use sqlite for files greater than max len
            self._type = 'sq'
            self._init_sqlite(bedfile)
        else:
            self._type = 'dt'
            self._init_dict(bedfile)
