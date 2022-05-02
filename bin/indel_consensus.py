#!/usr/bin/env python

import pandas as pd
import argparse
import os
import re

desc = 'module for merging indels of two VCF files, '
def parse_args():
    """
    Argument parser
    """
    parser = argparse.ArgumentParser(description=desc,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--vcf1", "-1", type=str, required=True,
                        help="First VCF file, keep all the lines, except when overlapping with second VCF "),
    parser.add_argument("--vcf2", "-2", type=str, required=True,
                        help="Second VCF file (reference for INDELS), keep all the lines"),
    parser.add_argument('--window', '-w', type=str, required=False,
                        help='The window to determine if two regions are overlapping', default=15)
    parser.add_argument('--stdout', required=False, action='store_true',
                        help='Output to stdout, default is: write to output file')
    # parse all arguments
    args = parser.parse_args()

    return args

def sort_vcf(body_vcf, contig):
    """
    sort the composite VCF on genomic position after merging the two input VCF
    :param body_vcf: body of VCF (VCF lines minus the header)
    :return: sorted body as list
    """
    # split vcf strings to convert to pandas dataframe
    body_vcf = [l.split('\t') for l in body_vcf]
    # convert to dataframe for sorting VCF lines on genomic position
    # change type of genomic pos column for numeric sorting
    body_df = pd.DataFrame(body_vcf)
    bodies = []
    for c in contig:
        chr_code = re.search(r'<ID=((\w)+(\.(\d))?),', c, re.I).group(1)
        body_contig = body_df[body_df[0] == chr_code]
        body_contig[1] = pd.to_numeric(body_df[1])
        body_contig.sort_values([1], inplace=True)
        body_contig[1] = body_df[1].apply(str)
        bodies.append(body_contig)

    body_df = pd.concat(bodies)
    # convert pandas df to list for merging with header
    body_list = ['\t'.join(x) for x in body_df.values.tolist()]

    return body_list

def compare(vcf1, vcf2, window):
    """
    compare two pairs of VCF lines, when overlapping, keep variants of second VCF
    :param vcf1: first VCF, throw away indels when overlapping with vcf2
    :param vcf2: reference VCF for Indels (keep these lines when overlapping with vcf1)
    :param window: range to determine if there is real overlap, default=15
    :return: new (sorted) consensus VCF
    """
    # function to determine end position of variant
    def end_pos(pos, ref, alt, info):
        if alt == '<DEL>':
            return int(re.search(r';END=((\d)+);', info).group(1))
        else:
            return int(pos) + (len(ref) - len(alt))

    vcf1_lines = open_vcf(vcf1)
    vcf2_lines = open_vcf(vcf2)
    # check if no variants are found
    if len(vcf2_lines[1]) == 0:
        output_vcf = '\n'.join(vcf1_lines[0] + vcf1_lines[0])
        return output_vcf

    new_header, contig = sort_header_lines(vcf1_lines[0], vcf2_lines[0])
    # add the reference lines to the VCF, we will keep all of those
    body_vcf = list(vcf2_lines[1])
    # append all the VCF lines of the first file to the new VCF
    for line_vcf1 in vcf1_lines[1]:
        app = True
        # chr_1, pos_1, id_1, ref_1, alt_1, _, filt1, info1, *n = line_vcf1.split('\t')
        chr_1, pos_1, id_1, ref_1, alt_1, _, filt1, info1 = line_vcf1.split('\t')[0], line_vcf1.split('\t')[1], line_vcf1.split('\t')[2], line_vcf1.split('\t')[3], line_vcf1.split('\t')[4], line_vcf1.split('\t')[5], line_vcf1.split('\t')[6], line_vcf1.split('\t')[7]
        for line_vcf2 in vcf2_lines[1]:
            chr_2, pos_2, id_2, ref_2, alt_2, _, filt2, info2 = line_vcf2.split('\t')[0], line_vcf2.split('\t')[1], line_vcf2.split('\t')[2], line_vcf2.split('\t')[3], line_vcf2.split('\t')[4], line_vcf2.split('\t')[5], line_vcf2.split('\t')[6], line_vcf2.split('\t')[7]
            # fetch start and end positions
            end_pos1 = end_pos(pos_1, ref_1, alt_1, info1)
            end_pos2 = end_pos(pos_2, ref_2, alt_2, info2)
            # make cutoff for Delly variants, if called deletion is too large, skip this variant
            if (int(end_pos2) - int(pos_2)) > 100000:
                continue
            ###############################################################################
            #                     GAP 1: 20________________________40                     #
            #       GAP 2: 0___________________25  GAP 2: 32__________________50          #
            #                                                                             #
            #  start 2 <= start 1 + window  AND    OR     start 2 < end 1 - window  AND   #
            #  end   2 > start 1                          end   2 >= end 1                #
            #                                      OR                                     #
            #                            start cov >= start gen - window AND              #
            #                            end cov   <= end gen                             #
            ###############################################################################
            if (int(pos_2) <= (int(pos_1) + window)) and (int(end_pos2) > int(pos_1)) or \
                    (int(pos_2) < (int(end_pos1) - window)) and (int(end_pos2) >= int(end_pos1)) or \
                    (int(pos_2) >= (int(pos_1) + window)) and int(end_pos2) <= (int(end_pos1)- window):
                # determine if variant is a deletion with ALT column = <DEL>
                # or: bp length in ALT column is smaller than bp length in REF column
                if (alt_1 == '<DEL>' or (len(alt_1) < len(ref_1))) and \
                        (alt_2 == '<DEL>' or (len(alt_2) < len(ref_2))):
                    # if lines are overlapping, do not append line of first VCF
                    app = False
                    break
        # if app(end) is false, do not append
        if app:
            body_vcf.append(line_vcf1)

    # sort VCF on genomic position
    body_list = sort_vcf(body_vcf, contig)
    # merge head and body
    new_vcf = '\n'.join(new_header + body_list)

    return new_vcf

def sort_header_lines(head1, head2):
    """
    sort the header lines of both input VCF files, throw away duplicates
    :param head1: header of VCF1
    :param head2: header of VFC2
    :return: sorted, merged unique VCF header
    """
    alt_lines = []
    filter_lines = []
    format_lines = []
    info_lines = []
    vcf_version = []
    contig_line = []
    header_line = []

    def remove_dups(x):
        return list(dict.fromkeys(x))

    def remove_one(x):
        return [x[0]]

    # subcategorize header lines in separate lists
    for head in [head1, head2]:
        for line in head:
            if line.startswith('##ALT'):
                alt_lines.append(line)
            elif line.startswith('##FILTER'):
                filter_lines.append(line)
            elif line.startswith('##FORMAT'):
                format_lines.append(line)
            elif line.startswith('##INFO'):
                info_lines.append(line)
            elif line.startswith('##fileformat'):
                vcf_version.append(line)
            elif line.startswith('##contig'):
                contig_line.append(line)
            elif line.startswith('#CHROM'):
                header_line.append(line)

    # remove duplicates and ensure unique header lines for contig line and column header
    lists = [alt_lines, filter_lines, format_lines, info_lines, vcf_version, contig_line]
    alt, filter, format, info, version, contig = list(map(remove_dups, lists))
    head = remove_one(header_line)
    new_header = version + alt + filter + format + info + contig + head

    return new_header, contig

def open_vcf(path, mode='r'):
    """
    open VCF file, seperate header and body lines
    :param path: path of VCF file
    :param mode: file open mode, default: read
    :return:
    """
    if not os.path.exists(path):
        raise Exception('Path of vcf file does not exist')
    with open(path, mode) as inp:
        lines = inp.readlines()
        # filter header and empty lines
        body_lines = [l.rstrip() for l in lines if not l.startswith('#') and l]
        head_lines = [l.rstrip() for l in lines if l.startswith('#')]
    return [head_lines, body_lines]

def write_vcf(new_vcf, vcf1, vcf2, stdout):
    """
    output VCF writer, write to file if stdout = False
    :param new_vcf: new consensus VCF
    :param vcf1: file path of VCF1
    :param vcf2: file path of VCF2
    :param stdout: if True, print to stdout
    """
    if stdout:
        print(new_vcf)
    else:
        vcf1_basename = os.path.basename(vcf1).rstrip('.vcf')
        vcf2_basename = os.path.basename(vcf2).rstrip('.vcf')
        # join VCF file names and construct new file path
        out_name = '{}_{}_merged.vcf'.format(vcf1_basename, vcf2_basename)
        with open(out_name, 'w') as out:
            out.write(new_vcf)

def run(args):
    # compare two VCF, keep reference lines when overlapping
    new_vcf = compare(args.vcf1, args.vcf2, args.window)
    # write out new consensus VCF
    write_vcf(new_vcf, args.vcf1, args.vcf2, args.stdout)

if __name__ == "__main__":
    args = parse_args()
    run(args)
