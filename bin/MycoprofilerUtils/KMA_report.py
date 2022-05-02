#!/usr/bin/env python

######
INFO = "Code to extract single species from KMA report file"
__version__ = 1.0
######

import argparse
import os
import logging
import utils.pandas as upd
import storage.storage as storage
import sys

def parse_args():
    """
        Argument parser
    """
    parser = argparse.ArgumentParser(description=INFO, \
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--inputFile", type=str, required=True,
                        help="input full path of input report file")
    parser.add_argument("--JSON", type=str, required=False,
                        help="full path to JSON file", default=None)
    parser.add_argument("--name", type=str, required=False,
                        help="name to store in JSON as tool name", default='KMA')
    parser.add_argument("--database", type=str, required=False,
                        help="name of database to put in JSON", default='KMA')
    parser.add_argument("--sample", type=str, required=True,
                        help="name of sequence sample"),
    parser.add_argument("-o", "--outputDir", type=str, required=False,
                        help="full path of output folder", default=os.path.abspath("./"))
    parser.add_argument("-n", "--numberofhits", type=int, required=False,
                        help="Number of top hits to store in data", default=10)
    parser.add_argument("--metric", type=str, required=False,
                        choices=['q_value','score'],
                        default='q_value',
                        help="metric to use for getting best hit")
    parser.add_argument("-v", "--version", action="version",
                    version="%(prog)s {version}".format(version=__version__))

    # parse all arguments
    args = parser.parse_args()

    return args

def KMA_filter(df, identity, coverage):
    '''
        :param df: dataframe file from pandas
        :param identity: minimal identity of match
        :param coverage: minimal coverage of match
        
        :result df: dataframe filtered on minimal identity and coverage
        
    '''
    
    # filter on identity
    df = df[(df['Template_Identity'] >= identity)]
    df = df[(df['Query_Identity'] >= identity)]
    # filter on coverage
    df = df[(df['Template_Coverage'] >= coverage)]
    df = df[(df['Query_Coverage'] >= coverage)]

    return df

def extract_best_match(df, column):
    """
        Extract row with best matching species
        and taxid of best match

        :param df: pandas dataframe in format output centrifuge
        :param column: which column to use as besthit metric
        :return name: name of besthit
    """

    if df.empty:
        name = "Failed"
        metric = "NA"
        depth = "NA"
        identity = "NA"
        coverage = "NA"
        return (name, metric, depth, identity, coverage)

    BestHit = df[df[column]==df[column].max()]
    name = str(BestHit['#Template'].values[0])
    metric = int(BestHit[column].values[0])
    depth = int(BestHit["Depth"].values[0])
    identity = int(BestHit["Template_Identity"].values[0])
    coverage = int(BestHit["Template_Coverage"].values[0])

    return(name, metric, depth, identity, coverage)

def run(args, name):
    """
        Main run function

        :param args: arguments
        :param name: samplename
    """
    print('Start')
    df = upd.read_report(args.inputFile)
    df = KMA_filter(df, 70, 70)
    besthit, metric, depth, identity, coverage = extract_best_match(df, args.metric)
    
    # list of tuples of results
    results = [("database", args.database),
               ("version", __version__),
               ("besthit", besthit),
               ("metric_score", metric),
               ("metric_type", args.metric),
               ("depth", depth),
               ("identity", identity),
               ("coverage", coverage),
               ("data", df.to_dict())]

    # if JSON is present use exiting, else create new unique name
    JSON = storage.JSON()
    if args.JSON != None:
        JSON.open(args.JSON)
    else:
        JSON.name(args.sample)

    JSON.add_results(args.name, results)
    JSON.pretty_print()
    JSON.write(args.outputDir)

    logging.info(results)
    
if __name__ == "__main__":
    # load arguments set global workdir
    args = parse_args()
    logging.info(args)
    name = args.inputFile.replace('.report','')
    name = name.replace('.','_')
    run(args, name)
    logging.info('Finished extract')
    sys.exit(0)