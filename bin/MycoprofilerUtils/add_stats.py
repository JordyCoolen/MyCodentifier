#!/usr/bin/env python

import argparse
import storage.storage as storage
import logging
import os

######
INFO = "Extract data and add calculations to JSON file"
__version__ = 0.10
######

def parse_args():
    """
        Argument parser
    """
    parser = argparse.ArgumentParser(description=INFO, \
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--JSON", type=str, required=True,
                        help="full path to JSON file", default=None)
    parser.add_argument("-o", "--outputDir", type=str, required=False,
                        help="full path of output folder", default=os.path.abspath("./"))
    parser.add_argument("-v", "--version", action="version",
                        version="%(prog)s {version}".format(version=__version__))

    # parse all arguments
    args = parser.parse_args()

    return args

def add_calculate_percentage_reads(JSON, key, readtype, inputreads):
    """
        Add ((value * 2) / all input reads) * 100 to the data of centrifuge

        :param JSON: JSON object file
        :param key: (str) name of the key in the dictionary/json
        :param readtype: (str) this is either numReads or numUniqueReads
        :param inputreads: (int) total number of input reads
    """
    # create key percentage_numReads
    JSON.data[key]['data'][f"%{readtype}"] = {"0": 0}
    # loop over the numReads and then calculate percentage using inputreads
    # result will be added to the dictionary
    for k, v in JSON.data[key]['data'][readtype].items():
        JSON.data[key]['data'][f"%{readtype}"].update({k: round(((v*2)/inputreads * 100), 2)})

    return JSON

def get_key_of_besthit(JSON, key, readtype):
    """
        Obtain the besthit via the taxID.
        This is needed to extract the correct readtype for
        calculating the percentage.

        :param JSON: JSON object file
        :param key: (str) name of the key in the dictionary/json
        :param readtype: (str) this is either numReads or numUniqueReads
    """
    taxid = JSON.data[key]['taxID']
    for k,v in JSON.data[key]['data']['taxID'].items():
        if v == taxid:
            return JSON.data[key]['data'][f"%{readtype}"][k]

def run(args):
    # if JSON is present use exiting, else create new unique name
    JSON = storage.JSON()
    JSON.open(args.JSON)

    # fastp calculations
    before = JSON.data['fastp']['summary']['before_filtering']['total_reads']
    after = JSON.data['fastp']['summary']['after_filtering']['total_reads']
    JSON.update('fastp', [("%reads", round(after / before * 100, 2))])

    JSON = add_calculate_percentage_reads(JSON, 'centrifuge', 'numReads', after)
    # get percentage of besthit
    cen_best = get_key_of_besthit(JSON, 'centrifuge', 'numReads')
    JSON.update('centrifuge', [("%numReads", cen_best)])

    JSON = add_calculate_percentage_reads(JSON, 'WGS_typing', 'numReads', after)
    WGS_best = get_key_of_besthit(JSON, 'WGS_typing', 'numReads')
    JSON.update('WGS_typing', [("%numReads", WGS_best)])

    JSON = add_calculate_percentage_reads(JSON, 'contaminants', 'numReads', after)
    con_best = get_key_of_besthit(JSON, 'contaminants', 'numReads')
    JSON.update('contaminants', [("%numReads", con_best)])

    JSON = add_calculate_percentage_reads(JSON, 'centrifuge', 'numUniqueReads', after)
    # get percentage of besthit
    cen_best_un = get_key_of_besthit(JSON, 'centrifuge', 'numUniqueReads')
    JSON.update('centrifuge', [("%numUniqueReads", cen_best_un)])

    JSON = add_calculate_percentage_reads(JSON, 'WGS_typing', 'numUniqueReads', after)
    WGS_best_un = get_key_of_besthit(JSON, 'WGS_typing', 'numUniqueReads')
    JSON.update('WGS_typing', [("%numUniqueReads", WGS_best_un)])

    JSON = add_calculate_percentage_reads(JSON, 'contaminants', 'numUniqueReads', after)
    con_best_un = get_key_of_besthit(JSON, 'contaminants', 'numUniqueReads')
    JSON.update('contaminants', [("%numUniqueReads", con_best_un)])

    # write and save JSNO file to disk
    JSON.write(args.outputDir)

if __name__ == "__main__":
    # load arguments set global workdir
    args = parse_args()
    logging.info(args)
    run(args)
    logging.info('Finished adding stats to JSON')