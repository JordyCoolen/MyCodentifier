#!/usr/bin/env python

"""
    SNP-IT (to identify Mycobacterium tuberculosis complex (MTC) lineage and sublineage
    https://github.com/philipwfowler/snpit
"""
INFO = "SNP-IT"
__version__ = 1.0

import logging, os
from utils.CommandRunner import exe
import utils.Functions as Functions
from string import Template
import storage.storage as storage
from os import path
import argparse
import sys

logging = logging.getLogger('test')

def parse_args():
    """
        Argument parser
    """
    parser = argparse.ArgumentParser(description=INFO,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--input', '-i', nargs='?', type=argparse.FileType('r'),
                    default=sys.stdin, help="SNPit tab delimited output")
    parser.add_argument("-o", "--outputDir", type=str, required=False,
                        help="full path of output folder", default=os.path.abspath("./"))
    parser.add_argument("--sample", type=str, required=True,
                        help="name of sequence sample"),
    parser.add_argument("--JSON", type=str, required=False,
                        help="full path to JSON file", default=None)
    parser.add_argument("--snpit", type=str, required=False,
                        help="relative path to snipit-run.py", default='snpit/bin/snpit-run.py')
    parser.add_argument("-v", "--version", action="version",
                    version="%(prog)s {version}".format(version=__version__))
    
    args = parser.parse_args()
    return args

def parse_snpit_output(args):
    """
    Parses stdout and add that to JSON file
        
    :param stdout: Output of SNP-IT (adjusted to output tabdelim)
    :param JSON: JSON dictionary file
    :return JSON: returns JSON object
    """
    # split stdout
    snpit_stdin = args.input.read()
    lines = snpit_stdin.split("\n")
    split = lines[1].split("\t") # take second line
    
    # parse stdout
    sample_name = split[0]
    species = split[1]
    lineage = split[2]
    sublineage = split[3]
    name = split[4]
    percentage = split[5]

    # list of results
    results = [('database', 'SNPIT'),
               ('version', __version__),
               ('sample_name', sample_name),
               ('species', species),
               ('lineage', lineage),
               ('sublineage', sublineage),
               ('name', name),
               ('percentage', percentage)]
    
    # if JSON is present use exiting, else create new unique name
    JSON = storage.JSON()
    if not args.JSON:
        JSON.name(args.sample)
    JSON.open(args.JSON)
    JSON.add_results('SNP-IT', results)
    JSON.pretty_print()
    JSON.write(args.outputDir)

    logging.info(results)


if __name__ == "__main__":
    # get all arguments
    args = parse_args()
    # error handling if pipeline fails with BUZZ feedback
    parse_snpit_output(args)