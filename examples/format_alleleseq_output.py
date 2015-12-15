#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# vim:fenc=utf-8 tabstop=4 expandtab shiftwidth=4 softtabstop=4
"""
#============================================================================#
#                                                                            #
#          FILE: make_sheet (python 3)                                       #
#        AUTHOR: Michael D Dacre, mike.dacre@gmail.com                       #
#  ORGANIZATION: Stanford University                                         #
#       LICENSE: MIT License, property of Stanford, use as you wish          #
#       CREATED: 2015-12-13 10:07                                            #
# Last modified: 2015-12-14 16:59                                            #
#                                                                            #
#   DESCRIPTION: Format the output of the alleleseq pipeline into gene and   #
#                snp level tab-delimited outputs sorted by tissue            #
#                                                                            #
#         USAGE: Provide a list of count files, in the same directory, there #
#                should exist the files <count_file>.FDR.txt and             #
#                <count_file>.interestingHets.txt                            #
#                                                                            #
#============================================================================#
"""
import os
import sys
import operator
from re import sub
from bed_lookup import BedFile  # github.com:MikeDacre/python_bed_lookup.git

# Defaults should be altered for your usage
tissue_lookup_file = 'tissue_lookup.txt'  # Format: dir\tcross\ttissue
master_lookup_file = 'master.lookup.v2'   # Format: gene_name\tmodel_name
bed_snp_file       = 'refseq.ucsc.ensembl.mRNA.mm9.nr.bed'


def create_tissue_lookup(tissue_file=tissue_lookup_file):
    """ Return a dictionary of sample data from the tissue
        lookup file. This should be altered for the sample
        data you want to extract. The index here is the name
        of the alleleseq output file up to the first period.
        In my case the name 123.cnt.counts.txt turns into 123 """
    tissue_lookup = {}
    with open(tissue_file) as infile:
        for line in infile:
            fields = line.rstrip().split('\t')
            assert(len(fields) == 4)
            tissue_lookup[fields[0]] = {'cross':   fields[1],
                                        'tissue':  fields[2],
                                        'failed':  fields[3]}
    return tissue_lookup


def create_master_lookup(master_file=master_lookup_file):
    """ Return a dictionary from a two column file. Used for
        replacing gene names extracted from the bed with other
        names. Used later for the 'model' parameter, which is
        what the gene-level data is indexed by """
    master_lookup = {}
    with open(master_file) as infile:
        for line in infile:
            fields = line.rstrip().split('\t')
            assert(len(fields) == 2)
            master_lookup[fields[0]] = fields[1]
    return master_lookup


def get_snp_data(count_files, tissue_file=tissue_lookup_file,
                 bed_file=bed_snp_file, master_file=master_lookup_file):
    """ Extract SNP-level data from the count files and add
        sample level information from the tissue_file (via
        create_tissue_lookup and gene names from bed_lookup
        and master_lookup (using a bed file and master gene
        name dictionary) """
    tissues = {}  # Dictionary to hold all data
    tissue_lookup = create_tissue_lookup(tissue_file)
    master_lookup = create_master_lookup(master_file)
    bed_lookup = BedFile(bed_file)
    for i in count_files:
        # Create entry for this count file
        t = i.split('.')[0]
        tissues[t] = {}
        # Lookup tissue-level data
        if t in tissue_lookup:
            tissues[t]['tissue'] = tissue_lookup[t]['tissue']
            tissues[t]['cross']  = tissue_lookup[t]['cross']
            tissues[t]['failed'] = tissue_lookup[t]['failed']
        else:
            tissues[t]['tissue'] = 'unknown'
            tissues[t]['cross']  = 'unknown'
            tissues[t]['failed'] = 'unknown'

        # Get list of hets that beat an FDR of 0.1
        hets = []
        with open(i + '.interestingHets.txt') as infile:
            for line in infile:
                if line.startswith('chrm'):
                    continue
                f = line.rstrip().split('\t')
                chr = f[0] if f[0].startswith('c') else 'chr' + f[0]
                snp = f[1]
                hets.append(chr + '_' + snp)

        # Extract SNPs from count file
        tissues[t]['snps'] = {}
        with open(i) as infile:
            for line in infile:
                if line.startswith('chrm'):
                    continue
                f = line.rstrip().split('\t')
                chr   = f[0] if f[0].startswith('c') else 'chr' + str(f[0])
                snp   = f[1]
                gene  = bed_lookup.lookup(chr, int(snp))
                gene  = gene if gene else ''
                model = master_lookup[gene] if gene else ''
                sig   = 'Y' if chr + '_' + snp in hets else 'N'
                id    = chr + '_' + snp
                tissues[t]['snps'][id] = {
                    'chr':       chr,
                    'snp':       snp,
                    'gene':      gene,
                    'model':     model,
                    'mat_gtyp':  f[7],
                    'pat_gtyp':  f[8],
                    'counts':    {'A': int(f[9]),
                                  'C': int(f[10]),
                                  'G': int(f[11]),
                                  'T': int(f[12]),
                                  'unknown': 'unknown'},
                    'win':       f[13],
                    'p':         f[15],
                    'beats_FDR': sig}

    return tissues


def snps_to_genes(snp_dict):
    """ Take a dictionary of snps by tissue from get_snp_data()
        and return a dictionary inventoried by tissue->gene
        instead of tissue->snp.
        'mat' and 'pat' data are converted to data named by
        parent and counts are indexed by parent rather than
        by base.
        In addition, a 'snp_count' entry is added to track the
        number of snps contributing to a count. """
    gene_dict = {}
    for t in snp_dict.keys():

        assert(t not in gene_dict)
        gene_dict[t] = {}
        gene_dict[t]['genes']  = {}
        gene_dict[t]['tissue'] = snp_dict[t]['tissue']
        gene_dict[t]['cross']  = snp_dict[t]['cross']
        gene_dict[t]['failed'] = snp_dict[t]['failed']

        for location, data in snp_dict[t]['snps'].items():
            # Assign genotype to single parent
            if snp_dict[t]['cross'] == 'CxB' or snp_dict[t]['cross'] == 'cxb':
                ca = data['pat_gtyp']
                b6 = data['mat_gtyp']
                ca_counts = data['counts'][ca]
                b6_counts = data['counts'][b6]
            elif snp_dict[t]['cross'] == 'BxC' or snp_dict[t]['cross'] == 'bxc':
                ca = data['mat_gtyp']
                b6 = data['pat_gtyp']
                ca_counts = data['counts'][ca]
                b6_counts = data['counts'][b6]
            elif snp_dict[t]['cross'] == 'unknown':
                ca = 'unknown'
                b6 = 'unknown'
                ca_counts = 'unknown'
                b6_counts = 'unknown'
            else:
                raise Exception

            m = data['model']  # Index by gene model name
            if m in gene_dict[t]['genes']:
                # This is a new SNP in an existing transcript
                gene_dict[t]['genes'][m]['snp_count'] += 1
                gene_dict[t]['genes'][m]['ca_counts'] = \
                    _combine_counts(gene_dict[t]['genes'][m]['ca_counts'], ca_counts)
                gene_dict[t]['genes'][m]['b6_counts'] = \
                    _combine_counts(gene_dict[t]['genes'][m]['b6_counts'], b6_counts)
                gene_dict[t]['genes'][m]['p_vals'].append(data['p'])
                gene_dict[t]['genes'][m]['beats_FDR'].append(data['beats_FDR'])
            else:
                gene_dict[t]['genes'][m] = {}
                gene_dict[t]['genes'][m]['snp_count'] = 1
                gene_dict[t]['genes'][m]['ca_counts'] = ca_counts
                gene_dict[t]['genes'][m]['b6_counts'] = b6_counts
                gene_dict[t]['genes'][m]['p_vals']    = [data['p']]
                gene_dict[t]['genes'][m]['beats_FDR'] = [data['beats_FDR']]
                gene_dict[t]['genes'][m]['gene']      = data['gene']

    return gene_dict


def print_snp_data(snp_dict, outfile=''):
    """ Print tab delimited data from output of make_data() """
    o = open(outfile, 'w') if outfile else sys.stdout
    o.write('DIR ID\tChr\tlocation\tGene\tTX\ttissue\tcross\tB6 Counts\tCAST Counts\t' +
            'p-value\tbeats FDR 0.1\tB6 gtyp\tCAST gtyp\tMat gTyp\tPat gTyp\t' +
            'A counts\tC counts\tG Counts\tT Counts\tWining Parent\tFailed\n')
    for t in snp_dict.keys():
        temp_dict = {}  # A dictionary to allow sorting for printing
        for k, d in snp_dict[t]['snps'].items():
            if snp_dict[t]['cross'] == 'CxB' or snp_dict[t]['cross'] == 'cxb':
                ca = d['pat_gtyp']
                b6 = d['mat_gtyp']
            elif snp_dict[t]['cross'] == 'BxC' or snp_dict[t]['cross'] == 'bxc':
                ca = d['mat_gtyp']
                b6 = d['pat_gtyp']
            elif snp_dict[t]['cross'] == 'unknown':
                ca = 'unknown'
                b6 = 'unknown'
            else:
                raise Exception

            print_string = [t, d['chr'], d['snp'], d['model'], d['gene'],
                            snp_dict[t]['tissue'], snp_dict[t]['cross'],
                            d['counts'][b6], d['counts'][ca], d['p'],
                            d['beats_FDR'], b6, ca, d['mat_gtyp'], d['pat_gtyp'],
                            d['counts']['A'], d['counts']['C'], d['counts']['G'],
                            d['counts']['T'], d['win'], snp_dict[t]['failed']]
            print_string = tuple([str(i) for i in print_string])
            if d['chr'] not in temp_dict:
                temp_dict[d['chr']] = {}
            temp_dict[d['chr']][int(d['snp'])] = print_string
        temp_dict = _make_sorted(temp_dict)
        for k, v in temp_dict:
            for i, print_string in v:
                o.write('\t'.join(print_string) + '\n')
    o.close()


def print_gene_data(gene_dict, outfile=''):
    """ Print tab delimited data from output of make_data() """
    o = open(outfile, 'w') if outfile else sys.stdout
    o.write('Gene\tTX\tTissue\tDIR ID\tCross\tB6 Counts\tCAST Counts\t' +
            'p-value\tbeats FDR 0.1\tSNP Count\tFailed\n')
    for t in gene_dict.keys():
        for k, d in gene_dict[t]['genes'].items():
            p_vals = ','.join(d['p_vals'])
            fdrs   = ','.join(d['beats_FDR'])
            print_string = [k, d['gene'], gene_dict[t]['tissue'], t,
                            gene_dict[t]['cross'], d['b6_counts'], d['ca_counts'],
                            p_vals, fdrs, d['snp_count'], gene_dict[t]['failed']]
            print_string = tuple([str(i) for i in print_string])
            o.write('\t'.join(print_string) + '\n')
    o.close()


def main(files, tissue_file=tissue_lookup_file, bed_file=bed_snp_file,
         master_file=master_lookup_file, outfile='', snp_outfile=''):
    """ Run everything """
    snp_dict  = get_snp_data(files, tissue_file, bed_file, master_file)
    gene_dict = snps_to_genes(snp_dict)
    if snp_outfile:
        print_snp_data(snp_dict, snp_outfile)
    print_gene_data(gene_dict, outfile)


#####################
# Don't alter these #
#####################
def _combine_counts(count1, count2):
    """ Sum two counts, but check that if one is 'unknown',
        both are 'unknown'. In those cases, return a single
        value of 'unknown'. """
    if count1 == 'unknown' or count2 == 'unknown':
        assert(count1 == 'unknown' and count2 == 'unknown')
        return 'unknown'
    assert(type(count1) == int and type(count2) == int)
    return count1 + count2


def _make_sorted(dict):
    """ Sort a dictionary for printing """
    print_dict = {}
    for k, v in dict.items():
        index = sub(r'chr', '', k)
        if index == 'X':
            index = '99'
        elif index == 'Y':
            index = '100'
        elif index == 'M' or index == 'MT':
            index = '101'
        print_dict[int(index)] = sorted(v.items(), key=operator.itemgetter(1))
    return sorted(print_dict.items(), key=operator.itemgetter(0))


###########################
# For running as a script #
###########################
if __name__ == '__main__' and '__file__' in globals():
    """ Command Line Argument Parsing """
    import argparse

    f_class = argparse.RawDescriptionHelpFormatter
    parser  = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=f_class)

    # Input files
    parser.add_argument('count_files', nargs='+',
                        help="Files created by AlleleSeq_pipeline_v1.2a/CombineSnpCounts.py")

    # Optional Files
    parser.add_argument('-t', '--tissue_file', default=tissue_lookup_file,
                        help="Tissue lookup file")
    parser.add_argument('-b', '--bed_snp_file', default=bed_snp_file,
                        help="Bed format gene lookup file")
    parser.add_argument('-m', '--master_lookup', default=master_lookup_file,
                        help="Model from gene lookup file")
    parser.add_argument('-o', '--gene_outfile', default='',
                        help="Output file, Default STDOUT")
    parser.add_argument('-s', '--snp_outfile', default='',
                        help="Also print SNP level data to this file. " +
                             "By default, SNP level data are not output.")

    args = parser.parse_args()

    # Run the script
    main(args.count_files, args.tissue_file, args.bed_snp_file, args.master_lookup,
         args.gene_outfile, args.snp_outfile)

##
# The End #
