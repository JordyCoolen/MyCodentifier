#!/usr/bin/env python

######
INFO = "Code to extract single species from centrifuge report file"
__version__ = 0.10
######

import argparse
import os
import utils.pandas as upd
import storage.storage as storage
import logging
logging = logging.getLogger('test')

import matplotlib as mpl
#mpl.use('TkAgg')
mpl.use('agg')
import matplotlib.pyplot as plt

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
                        help="name to store in JSON as tool name", default='centrifuge')
    parser.add_argument("--sample", type=str, required=True,
                        help="name of sequence sample"),
    parser.add_argument("-o", "--outputDir", type=str, required=False,
                        help="full path of output folder", default=os.path.abspath("./"))
    parser.add_argument("-n", "--numberofhits", type=int, required=False,
                        help="Number of top hits to store in data", default=10)
    parser.add_argument("--metric", type=str, required=False,
                        choices=['numUniqueReads','numReads','abundance'],
                        default='abundance',
                        help="metric to use for getting best hit")
    parser.add_argument('--conv_table', type=str, required=False, default=None,
                        help='full path to accession/taxid conversion table')
    parser.add_argument('--omit_accession', action='store_true', required=False, default=False,
                        help='This will omit retrieving of accesionnumbers using the\
                        conversion table')
    parser.add_argument("-v", "--version", action="version",
                    version="%(prog)s {version}".format(version=__version__))
    
    # parse all arguments
    args = parser.parse_args()
    
    if args.omit_accession == False and args.conv_table == None:
        raise Exception('conv_table parameter not set')
    
    return args

def extract_hits_centrifuge(df, number, column):
    """
        Extracts top hits from centrifuge.report file

        :param df: pandas dataframe file created with read_centrifuge_report
        :param number: number of top hits to extract
        :param column: the column to use as metric for ranking
        
        :return df: new dataframe of top hits
    """
    # sort on column
    df = df.sort_values(column, ascending=False)
    # obtain top
    df = df.iloc[0:number]
    return(df)

def bar_plot(df, column, name, outputDir):
    """
        Matplotlib to plot barplot of dataframe

        :param df: pandas dataframe (with x='name' and y=column)
        :param column: column name to plot
        :param name: name of plot
        :outputDir: directory to place output

        :return res: full path to outputFile
    """

    res = os.path.join(outputDir, "{}_bar.png".format(name))

    df.plot.bar(y=column, x='name')
    plt.savefig(res, bbox_inches = 'tight')
    return(res)

def get_accession(taxid, table_f):
    """
    create dictionary from conversion table to retrieve accession code with taxid
    """
    with open(table_f, 'r') as table:
        lines = table.read().splitlines()
        rows = {x.split('\t')[1].strip(): x.split('\t')[0].strip() for x in lines}
    if taxid not in rows:
        raise Exception("given taxid not in sequence genbank conversion table: {}"\
                        .format(taxid))
    return rows[taxid]

def extract_best_match(df, column, conv_table, omit_accession):
    """
        Extract row with best matching species
        and taxid of best match

        :param df: pandas dataframe in format output centrifuge
        :param column: which column to use as besthit metric
        :return name: name of besthit
        :return taxID: taxID number of besthit
        :return metric: column used to get besthit
        :return abundance: column with abundance as calculated by centrifuge
        :return numreads: column with the number of reads on hit
        :return numuniquereads: column with number of unique reads on hit
    """

    BestHit = df[df[column]==df[column].max()]
    name = str(BestHit['name'].values[0])
    taxID = int(BestHit['taxID'].values[0])
    abundance = int(BestHit['abundance'].values[0])
    numreads = int(BestHit['numReads'].values[0])
    numuniquereads = int(BestHit['numUniqueReads'].values[0])
    if omit_accession:
        acc = 'NA'
    else:
        acc = get_accession(str(taxID), conv_table)
    metric = round(float(BestHit[column].values[0]), 3)
    return(name, taxID, acc, metric, abundance, numreads, numuniquereads)

def run(args, name):
    """
        Main run function

        :param args: arguments
        :param name: samplename
    """

    metric = args.metric
    conv_table = args.conv_table

    df = upd.read_report(args.inputFile)
    df = extract_hits_centrifuge(df, args.numberofhits, metric)
    besthit, taxID, acc, metric_value, abundance, numreads, numuniquereads = extract_best_match(df, metric, conv_table, args.omit_accession)

    barplot = bar_plot(df, metric, name, args.outputDir)

    # list of tuples of results
    results = [("version", __version__),
               ("besthit", besthit),
               ("taxID", taxID),
               ("accession", acc),
               ("metric", metric_value),
               ("metric_type", metric),
               ("abundance", abundance),
               ("numReads", numreads),
               ("numUniqueReads", numuniquereads),
               ("data", df.to_dict()),
               ("barplot", barplot)]

    # if JSON is present use exiting, else create new unique name
    JSON = storage.JSON()
    if not args.JSON:
        JSON.name(args.sample)
    JSON.open(args.JSON)
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
