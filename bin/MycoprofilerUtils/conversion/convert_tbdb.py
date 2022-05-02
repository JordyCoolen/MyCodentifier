#!/usr/bin/env python

import argparse
import traceback
import simplejson as json
import sys
import re
import os

INFO = "module for converting TBDB (mycobacterium tuberculosis resistance genes DB " \
       "@ https://github.com/jodyphelan/tbdb/blob/master/tbdb.csv) and adding " \
       "sequencing depth of these regions from given BAM file, ouput is JSON"

# SAMPLE JSON
# { "target_cov":
#   "Rv1630": {
#     "gene_name": "rpsA",
#      mutations: [
#       "hgvs_c": "c.1349c>t".
#       "cov": 20-20 ]
#     "overall_coverage": [
#       "2<": 3,
#       "2-10": 5,
#       "10-20": 10,
#       "20-30": 56,
#       "30-60": 190,
#       "60>": 298 ]
#       }
# }


# BED input: (chr, start, end, depth)
# NC_000962.3	0	1	48
# NC_000962.3	1	3	49
# NC_000962.3	3	5	50
# NC_000962.3	5	15	51
# NC_000962.3	15	18	52


def parse_args():
    parser = argparse.ArgumentParser(description=INFO,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--gff', '-g', required=True,
                    help='GFF file to check VCF and resistance DB')
    parser.add_argument('--resDB', '-d', required=True,
                    help='resistance db in format: gene_or_locus\\thgvs_notation\\tdrug\n')
    parser.add_argument('--acc', '-a', nargs='?',
                    help='Accession code to select corresponding resistance DB for organism')
    parser.add_argument('--makeBED', '-b', required=False, action='store_true',
                        help='creates a BED file from the complete regions from TBDB genes')
    parser.add_argument('--bedcov', '-c', required=False, type=argparse.FileType('r'),
                        default=sys.stdin, help="output of bedtools genomecov -bg, default: stdin as pipe")
    parser.add_argument('--window', '-w', required=False, type=int,
                        help='Window (up/downstream) to use for creating a BED file')
    # parser.add_argument('--coverage')

    # parse all arguments
    args = parser.parse_args()
    return args


def GFF_check(id, return_chr=False):
    """
    check gene names and locus tags with entries in the reference GFF
    :param id: gene name or locus tag
    :param return_chr: return chr, start, end for BED file
    :return: locus_tag, gene_name, success (boolean)
    """
    # read the lines and skip headers
    for line in gff_lines:
        if not line.startswith("#"):
            if line.split('\t')[2].strip() == 'gene':
                chr = line.split('\t')[0]
                start, end = line.split('\t')[3], line.split('\t')[4]
                strand = line.split('\t')[6].strip()
                info = line.split('\t')[8]
                gene_id = re.search(r'ID=([0-9a-zA-Z_]*)(;)?', info).group(1)
                locus = re.search(r'locus_tag=([0-9a-zA-Z_]*)(;)?', info).group(1)
                if 'gene=' in info:
                    gene_name = re.search(r'gene=([0-9a-zA-Z_]*)(;)?', info).group(1)
                else:
                    gene_name = ''
                if id == locus or id == gene_name or id == gene_id:
                    if return_chr:
                        return chr, start, end, strand, True
                    return locus, gene_name, start, end, strand, True
    return '', '', '', '', '', False


def compare_genome_cov(gen_start, gen_end, mut=None, strand=None, window=150):
    # if no BED (genome cov file) was given, return empty strings
    if not genome_cov:
        return '', ''
    mut_end = None
    mutation = True if mut and strand else False
    if mutation:
        if mut.group('prot'):
            prot_pos = gen_start + (int(mut.group('prot')) * 3) + 1
            mut_start, mut_end = prot_pos, prot_pos + 3
        elif mut.group('indel1'):
            mut_start = gen_start + int(mut.group('indel1'))
            if mut.group('indel2'):
                mut_end = gen_start + int(mut.group('indel2'))
        elif mut.group('sub'):
            mut_start = gen_start + int(mut.group('sub'))
        elif mut.group('any_miss'):
            mut_start = gen_start + int(mut.group('any_miss'))
        else:
            # in the case of frameshift or large deletion, there is no specific position
            mut_start = None
        if not mut_end:
            # in the case of point mutations
            mut_end = mut_start
    ##########################################################################################
    #                     GENE 1: 20__________________________40                             #
    #            COV : 0_ _ _ _ _ 22 23_ _ _30 31_ _ _ 37 38_ _ _ _ _45                      #
    #                                                                                        #
    #  start cov <= start gene - window  AND     OR    start cov <= end gene - window  AND   #
    #  end cov  2 >= start gene                         end   cov >= end gene                #
    #                                            OR                                          #
    #                            start cov >= start gen - window AND                         #
    #                            end cov   <= end gen                                        #
    ##########################################################################################
    window_neg = window if strand == '-' else 0
    window_pos = window if strand == '+' else 0
    gen_dps = []
    mut_dps = []
    for row in genome_cov:
        chr, cov_start, cov_end, dp = row.split('\t')
        cov_start, cov_end, dp = int(cov_start), int(cov_end), int(dp)
        if not mutation:
            if (cov_start <= (gen_start - window_pos) <= cov_end or
                cov_start <= (gen_end + window_neg) <= cov_end or
                cov_start >= (gen_start - window_pos)) and cov_end <= gen_end + window_neg:
                gen_dps.append(dp)
        # if mutation has a position or if mutation is True, when False, only compute complete region
        if mutation and mut_start:
            if (cov_start <= mut_start <= cov_end) or \
                    (cov_start <= mut_end <= cov_end):
                mut_dps.append(dp)
        if cov_start > gen_end:
            break

    return mut_dps, gen_dps

def load_db(res_db_path, acc):
    """
    read resistance DB file and load into structured dict
    :param res_db_path: DB file path
    :return: structured dict with locus_mutation as keys
    """
    rows = get_db_path(acc, res_db_path)
    # remove empty strings from list
    rows = [x for x in rows if x]
    # initiate resistance dict with with key and set for all unique loci
    resdb = {'all_genes': {}, 'res_hgvs': {}, 'res_other': {}}
    # regex for all mutation types: protein, del/ins, substitution, frameshift, large del, large missense
    mut_pat = re.compile(r'(^[pn]\.-?[a-z]{3}(?P<prot>(\d){1,6})[a-z\*]{1,3})|'
                         r'(^[cn]\.-?(?P<indel1>(\d){1,6})(_-?(?P<indel2>(\d){1,6}))*(del|ins)[atcg]*)|'
                         r'(^[cn]\.-?(?P<sub>(\d){1,6})[atcg]*>[atcg]*)|(?P<large_del>large_deletion)|'
                         r'(?P<shift>frameshift)|any_missense_codon_(?P<any_miss>(\d)*)', re.I)
    # add the mutations to a mutation_drug dict
    for item in rows:
        gene, mut, drug, lit = [s.strip() for s in item.split('\t')]
        # replace r. notation with n. for compatibility with SNPEff (HGVS)
        mut = re.sub('^r.', 'n.', mut).lower()
        mut = mut.lower()
        # check GFF if locus/gene exists for the particular organism
        locus, gene, start, end, strand, succes = GFF_check(gene)
        if not succes:
            # stop execution if unknown gene in resistance DB
            raise Exception('Gene in resistance DB unknown:\n'+ gene)
        # make unique ids by combining locus and hgvs mutation
        unq_entry = '{}_{}'.format(locus, mut)
        # there are multiple drugs for the same mutations, check if mutation is \
        # already in the dictionary, if true: append info, else make new entry
        match = mut_pat.match(mut)
        if match:
            subdict = 'res_hgvs'
            if match.group('shift') or match.group('sub') or match.group('any_miss'):
                subdict = 'res_other'
            # find the coverage (depth) of mutation positions + whole gene region
            cov_mut, cov_overall = compare_genome_cov(int(start), int(end), match, strand)
            if unq_entry in resdb[subdict]:
                resdb[subdict][unq_entry]['drugs'].append(drug)
                if lit:
                    resdb[subdict][unq_entry]['lit'].append([lit])
            else:
                resdb[subdict][unq_entry] = {'drugs': [drug], 'lit': [lit]}
        # mutation must match the standard annotation, else raise error
        else:
            raise Exception('Unknown mutation type in database: {}, '
                            '\nplease modify script or correct DB'.format(resdb['other'][locus]))
        # add locus to all genes to keep a separate record of resistance genes and corresponding depth
        if locus not in resdb['all_genes']:
            cov_mut, cov_overall = compare_genome_cov(int(start), int(end))
            resdb['all_genes'][locus] = {'gene_name': gene, 'mutations': [], 'mutation_depths': [], 'coverage': cov_overall}
        else:
            resdb['all_genes'][locus]['mutation_depths'].append({'hgvs': mut, 'dp': cov_mut})

    return resdb

def get_db_path(acc, resDB):
    basename = '{}_resistance_db.tsv'
    dbpath = os.path.join(resDB, basename.format(acc))
    if os.path.exists(dbpath):
        return read_file(dbpath)
    else:
        print('no resistance db for this species, exiting..')
        exit(0)

def read_file(path):
    try:
        with open(path, 'r') as infile:
            lines = infile.readlines()
        return lines
    except Exception as e:
            exc_info = sys.exc_info()
            print(e)
            traceback.print_exception(*exc_info)
            sys.exit(1)

def write_results(db):
    """
    write the generated json with mutation info to pipeline json
    :param db: dict with converted mutation info from tbdb with coverage from bam
    """
    with open('resDB_cov.json', 'w') as dbout:
        json.dump(db, dbout)

def write_BED(bed_string):
    """
    write the TBDB gene regions to a BED format
    :param bed_string: string containing the regions in BED format
    :return:
    """
    with open('TBDB_regions.bed', 'w') as out:
        out.write(bed_string)

def TBDB_toBED(res_db_path, acc, window=200):
    """
    read resistance DB file and load into structured dict
    :param res_db_path: DB file path
    :param acc: accession code for determining which resistance DB to use
    :param window: window to use (up/downstream for creating BED file)
    :return: structured dict with locus_mutation as keys
    """
    rows = get_db_path(acc, res_db_path)
    # remove empty strings from list
    genes = [i.split('\t')[0].strip() for i in rows if i]
    genes = list(set(genes))
    # make al list for the BED file
    bed_list = []
    # loop over all the mutations to find corresponding position in GFF
    for gene in genes:
        # check GFF if locus/gene exists for the particular organism
        chr, start, end, strand, succes = GFF_check(gene, return_chr=True)
        if window:
            start -= window
            end += window
        # locus, gene_name, start, end, strand, succes = GFF_check(gene, return_chr=True)
        bed_line = '{}\t{}\t{}'.format(chr, start, end)
        if bed_line not in bed_list:
            bed_list.append(bed_line)
    bed_string = '\n'.join(bed_list)

    write_BED(bed_string)

def run(args):
    if args.makeBED:
        TBDB_toBED(args.resDB, args.acc, args.window)
    else:
        db = load_db(args.resDB, args.acc)
        write_results(db)


if __name__ == '__main__':
    args = parse_args()
    gff_lines = read_file(args.gff)
    #
    # if args.bedcov:
    #     r = args.bedcov.read()
    #     print(r)
    #     if not sys.stdin.isatty():
    #         print('yes')
    #     for line in args.bedcov:
    #         print(line)
    #     # for line in sys.stdin:
    #     #     print(line)
    #     # if sys.stdin.isatty():
    #     #     genome_cov = read_file(args.bedcov)
    # else:
    #     genome_cov = None
    genome_cov = None
    run(args)
