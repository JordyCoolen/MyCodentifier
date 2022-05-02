#!/usr/bin/env python
import re
import os
import sys
import glob2
import argparse
import simplejson as json
from datetime import datetime

INFO = 'Module for reporting resistant SNPs by mapping a VCF to a resistance database'

date = datetime.now().strftime('%Y%m%d%H%M%S')
out_file = "comparison_myco_tb_{}.txt".format(date)

samples_diff = []

def parse_args():
    parser = argparse.ArgumentParser(description=INFO, \
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--mycojson', '-m', nargs='?', type=argparse.FileType('r'),
                    help="JSON results file from the MycoProfiler pipeline")
    parser.add_argument('--mycodir', nargs='?',
                        help='directory of json results from myco pipeline')
    parser.add_argument('--tbjson', '-t', nargs='?', type=argparse.FileType('r'),
                    help="JSON results file from the TBProfiler pipeline")
    parser.add_argument('--tbdir', nargs='?',
                        help='directory of json results from TB pipeline')
    parser.add_argument('--sample', '-s', nargs='?',
                        help="sample name to include in output")
    parser.add_argument('--outdir', '-o', nargs='?', default='./',
                        help='output directory, default: current')
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    # parse all arguments
    args = parser.parse_args()
    return args

def update_out_path(outdir):
    global out_file
    out_file = os.path.join(outdir, out_file)

def get_sample_name(sample_path):
    sample_name = re.split(r'_(\d)*.*.json$', os.path.basename(sample_path[0]))[0]
    return sample_name

def write_header(not_paired):
    with open(out_file, 'a') as out:
        out.write(not_paired)

def write_footer():
    total = len(samples_diff)
    footer = "Total differences: {}\nSamples: {}".format(total, '\n'.join(samples_diff))
    with open(out_file, 'a') as out:
        out.write(footer)

def get_tb_mutations(json_rec):
    json_doc = json.load(json_rec)
    mutations = json_doc['dr_variants']
    tb_muts = []
    for mut in mutations:
        print(mut)
        # replace r. notation with c. for compatibility
        # n supposed to be c?
        hg_mut = re.sub('^r.', 'n.', mut['change']).lower()
        tb_muts.append([mut['locus_tag'].lower(), hg_mut.lower()])

    return tb_muts

def get_myco_mutations(json_rec):
    json_doc = json.load(json_rec)
    loci = json_doc['variants']['resistance']
    myco_muts = []
    for locus in loci:
        for mut in loci[locus]['mutations']:
            # replace r. notation with c. for compatibility
            #mutc = re.sub('^n.', 'c.', mut['hgvs_c'])
            mutc = mut['hgvs_c']
            myco_muts.append([locus.lower(), mutc.lower(), mut['hgvs_p'].lower()])

    return myco_muts

def compare_results(tb_muts, myco_muts, sample):
    # resolve which notation we need to compare for each mutation
    for i, myco_mut in enumerate(myco_muts):
        myco_muts[i] = [myco_muts[i][0], myco_muts[i][1]]
        success = False
        for tb_mut in tb_muts:
            if myco_mut[0] == tb_mut[0]:
                if myco_mut[2] == tb_mut[1]:
                    myco_muts[i] = [myco_mut[0], myco_mut[2]]
                    success = True
                    break
                elif myco_mut[1] == tb_mut[1]:
                    myco_muts[i] = [myco_mut[0], myco_mut[1]]
                    success = True
                    break
        if not success:
            pass

    myco_muts = [(x,y) for x,y in myco_muts]
    tb_muts = [(x,y) for x,y in tb_muts]

    myco_muts = set(tuple(myco_muts))
    tb_muts = set(tuple(tb_muts))

    myco_tb_diff = list(myco_muts.difference(tb_muts))
    tb_myco_diff = list(tb_muts.difference(myco_muts))

    myco_tb_diff = '' if not myco_tb_diff else myco_tb_diff
    myco_tb_diff = [', '.join(item) for item in myco_tb_diff]

    output_txt = []

    if myco_tb_diff:
        myco_not_tb = 'Mutations detected by MycoProfiler, but not by TBprofiler:\n{}\n'\
            .format('\n'.join(myco_tb_diff))
    else:
        myco_not_tb = 'Mutations detected by MycoProfiler, but not by TBprofiler:\nNone'

    output_txt.append(myco_not_tb)

    tb_myco_diff = [', '.join(item) for item in tb_myco_diff]
    tb_myco_diff = '' if not tb_myco_diff else tb_myco_diff

    if tb_myco_diff:
        tb_not_myco = 'Mutations detected by TBprofiler, but not by MycoProfiler:\n{}\n'\
            .format('\n'.join(tb_myco_diff))
    else:
        tb_not_myco ='Mutations detected by TBprofiler, but not by MycoProfiler:\nNone'

    output_txt.append(tb_not_myco)

    if not myco_tb_diff or not tb_myco_diff:
        output_txt = ['No difference found']
    else:
        global samples_diff
        samples_diff.append(sample)

    return output_txt

def get_parallel_lists(mycodir, tbdir):
    mycofiles = glob2.glob(os.path.abspath(mycodir) +'/*/mutations/*.json')
    tbfiles = glob2.glob(os.path.abspath(tbdir) +'/*/results/*.json')

    paired_json = [[my, tb] for tb in tbfiles for my in mycofiles if
                   re.split(r'_(\d)*.*.json$', os.path.basename(my))[0] == os.path.basename(tb).split('.results.json')[0]]

    paired_tb = [item[1] for item in paired_json]
    paired_my = [item[0] for item in paired_json]

    not_paired_tb = [tb for tb in tbfiles if tb not in paired_tb]
    not_paired_my = [my for my in mycofiles if my not in paired_my]

    not_paired = ''

    if not_paired_tb:
        not_paired += 'TBProfiler did analyze: {}, but MycoProiler did not\n'\
                       .format('\n'.join(not_paired_tb))
    if not_paired_my:
        not_paired += 'mycoprofiler did analyze {}, but tb-profiler did not\n'\
                      .format('\n'.join(not_paired_my))

    if not_paired:
        print(not_paired)

    sample_names = list(map(get_sample_name, paired_json))

    return paired_json, sample_names, not_paired


def write_results(output, samplename, tbfile, mycofile):
    output_txt = ['tbprofiler: {}\nmycoprofiler: {}\n'.format(tbfile, mycofile)] +\
                 ['sample: {}\n'.format(samplename)] + output + \
                 ['______________________________________________________________\n']
    output_str = '\n'.join(output_txt)

    with open(out_file, 'a') as out:
        out.write(output_str)

def run(args):
    update_out_path(args.outdir)
    if args.mycodir and args.tbdir:
        myco_json, sample_names, not_paired = get_parallel_lists(args.mycodir, args.tbdir)
        write_header(not_paired)
        for item, name in zip(myco_json, sample_names):
            tb_muts = get_tb_mutations(open(item[1]))
            myco_muts = get_myco_mutations(open(item[0]))
            output_txt = compare_results(tb_muts, myco_muts, name)
            write_results(output_txt, name, item[1], item[0])
        write_footer()
    else:
        if not args.sample:
            args.sample = get_sample_name(args.mycojson)
        tb_muts = get_tb_mutations(args.tbjson)
        myco_muts = get_myco_mutations(args.mycojson)
        output_txt = compare_results(tb_muts, myco_muts)
        write_results(output_txt, args.sample, args.tbjson, args.mycojson)


if __name__ == "__main__":
    args = parse_args()
    run(args)