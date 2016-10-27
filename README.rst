#################
Python Bed Lookup
#################

Allows very fast searching of a bed file of any size by gene/snp location.

For example:

.. code:: python

    from bed_lookup import BedFile
    b = BedFile('my_bed.bed')
    gene = b.lookup('chr3', 1000104)

This module requires cython, and should work with recent versions of
python2 and python3.

It can also be used with a pandas dataframe directly:

.. code:: python

   df['new_col'] = b.lookup_df(df, 'chrom', 'pos')

Note: with large dataframes, this function can be very slow, but there is
a nice trick to speed it up:

.. code:: python

   import numpy as np
   import pandas as pd
   from multiprocessing import Pool, cpu_count
   pool = Pool()
   b    = BedFile('my_bed.bed')
   df   = pd.read_csv('big_table.txt.gz', sep='\t', compression='gzip')
   dfs  = np.array_split(df, cpu_count())
   run  = []
   out  = []
   # Our chromsome column is 'chrom' and position column is 'pos'.
   for d in dfs:
       run.append(pool.apply_async(b.lookup_df, (d, 'chrom', 'pos')))
   for r in run:
       out.append(r.get())
   df['new_col'] = pd.concat(out)


************
Installation
************

Installation follows the standard python syntax:

.. code:: shell

    git clone https://github.com/MikeDacre/python_bed_lookup
    cd python_bed_lookup
    python setup.py build
    sudo python setup.py install

If you do not have root permission on you device, replace the last line with::

   python setup.py install --user

*****************************
Running from the command line
*****************************

There is a command line script called ``bed_location_lookup`` that will be installed
in ``/usr/bin`` if you install globally or in ``~/.local/usr/bin`` if you install for
your user only. The sytax for that script is::

    bed_location_lookup <bed_file> chr1_1000134 chr2_1859323 ....

It will work for any number of gene coordinate arguments. Be aware, that there is a
file opening delay when the script is run (for small bed files this will be very
small, but for large files it can be a few seconds). It is therefore much more
efficient to call a single instance of ``bed_location_lookup`` with a long list of
coordinates than it is to call it once per coordinate. For a large number of
coordinates this difference can be substantial.

``bed_location_lookup`` has a few other options also, to get those run::

    bed_location_lookup -h

Note: if you know the bed file is large and a database already exists, you can
get considerable speed up by passing the database file instead of the raw bed
file. e.g. pass ``bedfile.bed.db`` instead of ``bedfile.bed``. This bypasses the
file length check.

*************************************
Backend information and customization
*************************************

It makes use of a cython optimized dictionary lookup for small bed files
and sqlite for larger bed files. Which backend is being used is transparent
to the user, simply use the ``lookup()`` function as demonstrated in the
example above. The default file size cutoff is ~5 million lines in the bed
file, which results in a memory use of 1.2GB for a 5 million line long file.
The memory use scales linearly, so setting the limit at 1 million lines will
result in about 240MB of memory use. To change the file size cutoff edit the
``_max_len`` variable in ``bed_lookup/__init__.py``. Be aware that the file
size limit is actually measured in bytes, for speed purposes. A dictionary of
size to file length maps is provided in the ``__init__.py`` file, the default
should work fine on most systems.

Note that the sqlite backed is very slightly slower for lookups, however the
sqlite backend requires that a database exists already. If one does not exist
(the expected name is the bed file name followed by a ``.db``) already then one
is created, and this step can be very slow. Hypothetically this should only be
done once.

As noted above, when creating a BedFile object, a file length lookup is performed.
This lookup can be costly, particularly for gzipped files. To skip this step,
simply pass the database file to BedFile(), instead of the bedfile itself.

Note: this code will work with either plain text or gzipped files, gzipped files
will be slightly slower at load due to the overhead of decompression. For large
files where an sqlite database already exists, there will be only a very slight
delay relative to the uncompressed bed file (due to file length counting).

As the BedFile object is only generated once, any lookups after the creation of
this object will be very fast (less than a second) for *any* length of bed file.
Smaller files will obvious result in even quicker lookups.
