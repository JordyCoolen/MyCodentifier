#!/usr/bin/env python
import argparse
import requests
import traceback
import sys
import os
import re
import ftplib
from time import sleep

desc = "Tool for automatic download of Fasta Genbank and GFF files"
__version__ = 0.20

def parse_args():
    """
        Argument parser
    """
    parser = argparse.ArgumentParser(description=desc, \
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--taxid", type=str, required=False,
                        help="taxid in order to retrieve accession number from conversion table"),
    parser.add_argument("--table", type=str, required=False,
                        help="path to conversion table to convert taxid to accession"),
    # parser.add_argument("--standir", type=str, required=False,
    #                      help="standard output directory to check if genome files already exist")
    parser.add_argument('--accession', type=str, required=False,
                        help="accession code to retrieve file without conversion table"),
    parser.add_argument("-f", "--format", type=str, required=True,
                            choices=['genbank','fasta', 'gff', 'all'])
    parser.add_argument("--sample", type=str, required=False,
                            help="name of sequence sample"),
    parser.add_argument("-o", "--outdir", type=str, required=False,
                        help="full path of output folder", default=os.path.abspath("./"))

    parser.add_argument("-v", "--version", action="version",
                    version="%(prog)s {version}".format(version=__version__))

    # parse all arguments
    args = parser.parse_args()

    return args

def convert_id(taxid, table_f):
    """
    create dictionary from conversion table to retrieve accession code with taxid
    """
    with open(table_f, 'r') as table:
        lines = table.read().splitlines()
        rows = {x.split('\t')[1]: x.split('\t')[0] for x in lines}
    if taxid not in rows:
        raise Exception("given taxid not in sequence genbank conversion table")
    return rows[taxid]

def fetch_meta_info(acc, form='gb'):
    """
    retrieve genbank file from NCBI with requests library, use gbk file to retrieve meta data for fetching genome files
    """
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&rettype={}&id={}&retmode=text"
    form_url = url.format(form, acc)
    timer = 0
    try:
        r = requests.get(form_url)
        while r.status_code != 200:
            r = requests.get(url)
            sleep(10)
            timer += 10
            # stop at 10 tries, raise exception
            if timer > 600:
                raise Exception('Too many tries (wait time 10 min) for request: {}'.format(form_url))
        meta_info = r.text
        # get Assembly ID of genome from gbk file, use this to resolve FTP path of corresponding genomes
        match = re.search(r'Assembly:(\s)*((\w)*_(\d)*\.(\d){1,2})\n', meta_info)
        if not match:
            raise Exception('assembly id not in genbank file: {}'.format(form_url))
        else:
            return match.group(2)
    except requests.exceptions.RequestException as e:
        print(e)
        exit(1)

def fetch_records(form, acc, outdir):
    """

    :param form: type of record: GFF, Fasta or Genbank
    :param acc: accession code
    :param outdir: output directory for genome files
    :return:
    """
    assembly_id = fetch_meta_info(acc)

    digits = re.search(r'GCF_((\d)*)\.(\d)', assembly_id).group(1)
    ftp_parent = '/genomes/all/GCF/{}/{}/{}/'.format(digits[0:3], digits[3:6], digits[6:9])

    ftp_connect = ftplib.FTP(host='ftp.ncbi.nlm.nih.gov', user='anonymous', passwd='heleen.severin@radboudumc.nl')

    # get assembly version subfolders, check most recent one
    dir_list = list(ftp_connect.mlsd(ftp_parent))
    assem_folders = [x for x in dir_list if x[0] not in ['.', '..']]

    # obtain newest assembly version
    versions = dict()
    for x in assem_folders:
            versions[x[0]] = x[0][-1]

    recent_assem = max(versions, key=versions.get)

    # this is not used because of not correctly modifying dates at NCBI
    # mod_dates = [x[1]['modify'] for x in assem_folders]
    # recent_assem = assem_folders[mod_dates.index(max(mod_dates))][0]

    target_folder = os.path.join(ftp_parent, recent_assem)
    print('fetching genome file from: ftp.ncbi.nlm.nih.gov{}'.format(target_folder))
    target_file = os.path.join(target_folder, '{}{}'.format(recent_assem, form))

    write_out(target_file, ftp_connect, acc, form, outdir)


def write_out(addr, ftp_handle, acc, form, outdir):
    """
    write request object to record with appropriate name in chosen folder
    """
    out_name = "{}{}".format(acc, form)
    out_path = os.path.join(outdir, out_name)

    # make output path if not exists
    if not os.path.exists(outdir):
        os.mkdir(outdir, 0o755)
    # write out new genome file
    with open(out_path, 'wb') as out_file:
        print('writing out: {}'.format(out_path))
        ftp_handle.retrbinary('RETR %s' % addr, out_file.write)


def check_files_existance(acc, format, basepath):
    """
    Check if genome files are already present in standard folder
    :param acc: accession code
    :param format: extention of genome file
    :param basepath: standard path for genome files
    :return:
    """
    all_present = True
    file_paths = []
    exts = {'gff': ['_genomic.gff.gz'], 'fasta': ['_genomic.fna.gz'],
            'genbank': ['_genomic.gbff.gz'],
            'all': ['_genomic.gff.gz', '_genomic.fna.gz', '_genomic.gbff.gz']}

    if os.path.exists(basepath):
        for ext in exts[format]:
            # recreate standard file name format to perform check
            file_name = '{}{}'.format(acc, ext)
            file_path = os.path.join(basepath, file_name)
            if not os.path.exists(file_path):
                print('not present in standard dir: {}'.format(file_path))
                all_present = False
            file_paths.append(file_path)
        # if all the needed files are present in the given directory, cp to current
        # if not all present, program will not exit and download all the needed files
        if all_present:
            print("All reference genome files are present, unzip to current dir....")
            # for f in file_paths:
            #     copy(f, './')
            sys.exit(0)

def run(args):
    """
    main run function to control flow based on arguments
    """
    try:
        outdir = args.outdir
        if args.accession:
            acc = args.accession
        # if taxid is given as input, conversion to accession code is needed
        else:
            if not args.taxid:
                raise Exception("provide a taxid if input is not accession")
            if not args.table:
                raise Exception("provide a conversion table if taxid is input")
            acc = convert_id(args.taxid, args.table)

        # check if files already exist in the given outputdir
        if args.outdir:
            check_files_existance(acc, args.format, args.outdir)
        # download file format based on 'format' argument
        if args.format == 'fasta':
            fetch_records('_genomic.fna.gz', acc, outdir)
        elif args.format == 'genbank':
            fetch_records('_genomic.gbff.gz', acc, outdir)
        elif args.format == 'gff':
            fetch_records('_genomic.gff.gz', acc, outdir)
        else:
            fetch_records('_genomic.fna.gz', acc, outdir)
            fetch_records('_genomic.gff.gz', acc, outdir)
            fetch_records('_genomic.gbff.gz', acc, outdir)

    except Exception as e:
        exc_info = sys.exc_info()
        print(e)
        traceback.print_exception(*exc_info)
        sys.exit(1)


if __name__ == "__main__":
    # load arguments, run main code
    args = parse_args()
    run(args)
