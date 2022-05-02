#!/usr/bin/env python

from os import path
import glob2
import re
import sys
import argparse
import simplejson as json
import pandas as pd
from datetime import datetime

INFO = "Scipt to convert result json to a csv, can be used to extract only identification data or " \
       "experimental include the QC and resistance data"

__version__ = 0.3

drug_dict = {'ethionamide',
             'rifampicin',
             'amikacin',
             'delamanid',
             'isoniazid',
             'kanamycin',
             'ofloxacin',
             'capreomycin',
             'levofloxacin',
             'aminoglycosides',
             'ciprofloxacin',
             'fluoroquinolones',
             'cycloserine',
             'streptomycin',
             'para-aminosalicylic_acid',
             'ethambutol',
             'moxifloxacin',
             'pyrazinamide',
             'linezolid',
             'clofazimine',
             'bedaquiline'
             }

date = datetime.now().strftime('%Y%m%d%H%M%S')
out_file = "myco_collate_{}_{}.csv"

def parse_args():
    """
        Argument parser
    """
    parser = argparse.ArgumentParser(description=INFO, \
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--json", type=str, required=False,
                        help="full path to JSON file", default=None)
    parser.add_argument("--sample", type=str, required=False, default=None,
                        help="name of sequence sample"),
    parser.add_argument("--mycoresults", type=str, required=False,
                        help='directory of mycoprofiler results for batch analysis (containing jsons in the subfolders)')
    parser.add_argument("-o", "--outdir", type=str, required=False,
                        help="full path of output folder", default=path.abspath("./"))
    parser.add_argument("--onlyidentification", required=False, default=False, action='store_true',
                        help="only rapports identifications to csv")
    parser.add_argument("-v", "--version", action="version",
                    version="%(prog)s {version}".format(version=__version__))

    # parse all arguments
    # args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()

    return args

def get_sample_name(sample_path):
    sample_name = re.split(r'_(\d)*.*.json$', path.basename(sample_path))[0]
    return sample_name

def generate_csv(data, samplename):

    # check if expected keys are present
    #check_contents(data)

    # get detected resistances
    RES = data['variants']['resistance']

    # create pandas dataframe of resistance mutation
    RES_df = get_resistance_mutations(RES)

    # create susceptibility table
    RStable, mut_table = susceptibility(RES_df, data, samplename)

    # # obtain pandas dataframe of target gene mutations
    # RGENE = data['variants']['resistance_gene']
    # RGENE_df = get_resistance_gene(RGENE)
    #
    # # get detected resistances
    # RES = data['variants']['all_drugs']

    return RStable, mut_table

def sample_info(sample, data):
    '''
        Extract data from json and convert to dataframe.

        :params sample: (str) name of sample
        :params data: full JSON data file

        :return df: dataframe containing sample info

    '''
    # transform to dataframe
    df = pd.DataFrame(data={'Sample': sample,
                      'id': [data['id']],
                      'runID': [data['fastq']['runID']],
                      'barcode': [data['fastq']['barcode']]})
    return df

def get_fastp(data):
    '''
        Extract data from json and convert to dataframe.

        :params data: full JSON data file

        :return df: dataframe containing the fastp percentage

    '''
    # calculate percentage and add to dataframe
    df = pd.DataFrame(data ={'raw_reads': [data['fastp']['summary']['before_filtering']['total_reads']],
                             'filtered_reads': [data['fastp']['summary']['after_filtering']['total_reads']],
                             '%reads': [data['fastp']['%reads']]})
    return df

def get_centrifuge_data(data, key):
    '''
        Extract data from json and convert to dataframe.

        :params data: full JSON data file
        :params key: key to extract centrifuge data from

        :return df: dataframe containing the contaminants data

        keys: centrifuge, WGS_typing, contaminants,

    '''
    # transform to dataframe
    df = pd.DataFrame(data={key: [data[key]['besthit']],
                            'metric_contaminants': [data[key]['metric']],
                            '%numReads': [data[key]['%numReads']],
                            '%numUniqueReads': [data[key]['%numUniqueReads']]})
    return df

def extract_besthit(data, key, metric):
    df = pd.DataFrame.from_dict(data[key]['data'])
    BestHit = df[df[metric] == df[metric].max()]
    name = str(BestHit['name'].values[0])
    taxID = int(BestHit['taxID'].values[0])
    abundance = float(BestHit['abundance'].values[0])
    numreads = float(BestHit['%numReads'].values[0])
    numuniquereads = float(BestHit['%numUniqueReads'].values[0])
    # transform to dataframe
    res = pd.DataFrame(data={f'{key}_{metric}': [name],
                            f'taxID_{metric}': [taxID],
                            f'metric': [metric],
                            f'abundance_{metric}': [abundance],
                            f'%numReads_{metric}': [numreads],
                            f'%numUniqueReads_{metric}': [numuniquereads]})
    return(res)

def get_hsp65(data):
    '''
        Extract data from json and convert to dataframe.

        :params data: full JSON data file

        :return df: dataframe containing the hsp65 data

    '''
    # transform to dataframe
    df = pd.DataFrame(data={'hsp65': [data['KMA']['besthit']],
                            'depth': [data['KMA']['depth']],
                            'identity': [data['KMA']['identity']],
                            'coverage': [data['KMA']['coverage']]})
    return df

def get_SNPIT(data):
    '''
        Extract data from json and convert to dataframe.

        :params data: full JSON data file

        :return df: dataframe containing the SNPIT data
                    if SNP-IT is empty add NA
    '''

    # transform to dataframe
    try:
        df = pd.DataFrame(data={'SNPIT': [data['SNP-IT']['name']],
                            'percentage': [data['SNP-IT']['percentage']]})
    except KeyError:
        df = pd.DataFrame.from_dict(data={'SNPIT': ['NA'], 'percentage': ['NA']})
    return df

def get_QC(data):
    '''
        Extract data from json and convert to dataframe.

        :params data: full JSON data file

        :return df: dataframe containing the QC

    '''
    # transform to dataframe
    df = pd.DataFrame([data['BAM_QC']])
    return df

def get_info(data, sample):
    '''
        Perform all functions to extract data from json
        and convert to dataframe.

        :params data: full JSON data file
        :params sample: (str) sample name

        :return df: dataframe containing
                    all identification results
    '''
    sampleinfo = sample_info(sample, data)
    fastp = get_fastp(data)
    iden1 = extract_besthit(data, "centrifuge", 'abundance')
    iden2 = extract_besthit(data, "centrifuge", 'numReads')
    iden3 = extract_besthit(data, "centrifuge", 'numUniqueReads')
    cont1 = extract_besthit(data, "contaminants", 'abundance')
    cont2 = extract_besthit(data, "contaminants", 'numReads')
    cont3 = extract_besthit(data, "contaminants", 'numUniqueReads')
    hsp65 = get_hsp65(data)
    snpit = get_SNPIT(data)
    wgstyp1 = extract_besthit(data, "WGS_typing", 'abundance')
    wgstyp2 = extract_besthit(data, "WGS_typing", 'numReads')
    wgstyp3 = extract_besthit(data, "WGS_typing", 'numUniqueReads')
    qc = get_QC(data)
    df = pd.concat([sampleinfo, fastp, iden1, iden2, iden3, hsp65,
                        wgstyp1, wgstyp2, wgstyp3, snpit, cont1, cont2, cont3, qc], axis=1, sort=False)
    return df

def susceptibility(df, data, samplename):
    '''
        Get all drugs in resistance database
        and compare with found resistance.
        Make susceptibility table.

        :params df: dataframe obtained from get_resistance_mutations
        :params data: full JSON data file

        :return df: dataframe containing
                    drugs susceptibility #mutations
    '''

    #print(['sample'] + data['variants']['all_drugs'])
    # create pandas dataframe drug summary
    columns = list(df.columns.values)

    SUM = df['drug'].value_counts().reset_index()
    SUM.columns = ['drug', '#mutations']
    header_line = ['sample'] + data['variants']['all_drugs']
    table = []
    # compare found resistance to all_drugs list
    for d in data['variants']['all_drugs']:
        if d in SUM['drug'].to_list():
            i = SUM['drug'].to_list().index(d)
            table.append([d, 'R', SUM['#mutations'][i]])
        else:
            table.append([d, 'S', 0])

    report_df = pd.DataFrame(table).sort_values([1, 0])
    report_df.columns = ['drugs', 'Susceptibility', '#mutations']

    # header line
    SRtable_line = [samplename]
    muttable_line = [samplename]

    for drug in data['variants']['all_drugs']:
        if drug in SUM['drug'].to_list():
            SRtable_line.append('R')
        else:
            SRtable_line.append('S')

    for drug in data['variants']['all_drugs']:
        drug_col = ''
        for index,row in df.iterrows():
            if row['drug'] == drug:
                if row['aa']:
                    drug_col += row['aa']
                else:
                    drug_col += row['nt']
        muttable_line.append(drug_col)


    # add line as row in dataframe, not column!
    RStable = pd.DataFrame(columns=header_line)
    RStable.loc[0] = SRtable_line

    muttable = pd.DataFrame(columns=header_line)
    muttable.loc[0] = muttable_line
    muttable.columns = header_line

    return RStable, muttable

def get_resistance_mutations(data):
    '''
        parse JSON data to pandas dataframe

        examples of data:
        data['variants']['resistance']


        :params data:   data['variants']['resistance']
                        all mutations in resistance database

        :return df: [locus, gene, nt, aa, eff, drug] in pandas dataframe
    '''

    # resistance mutation table
    table = []

    # loop to obtain mutation per locus per drug
    for locus in data:
        gene = data[locus]["gene_name"]

        for i, _ in enumerate(data[locus]['mutations']):
            nt = data[locus]['mutations'][i]['hgvs_c']
            aa = data[locus]['mutations'][i]['hgvs_p']
            eff = data[locus]['mutations'][i]['eff']
            dp = data[locus]['mutations'][i]['dp']
            for d in data[locus]['mutations'][i]['drugs']:
                row = [locus, gene, dp, nt, aa, eff, d]
                table.append(row)

    # generate pandas table
    df = pd.DataFrame(table, columns=['locus', 'gene', 'depth', 'nt', 'aa', 'eff', 'drug'])
    df = df.sort_values(['drug'])

    return df

def get_json_list(mycodir):
    jsonfiles = glob2.glob(path.join(path.abspath(mycodir), '*.json'))
    return jsonfiles

def write_df_out(df, outdir, typecollate):
    out_path = path.join(outdir, out_file.format(typecollate, date))
    with open(out_path, 'w') as out:
        print("writing out: {}".format(out_path))
        df.to_csv(out, header=True, index=False, sep=',')

def run(args):
    '''
        General run code
    '''
    
    # for generating a csv on batch json
    if args.mycoresults:
        json_files = get_json_list(args.mycoresults)

        samplefile = json_files[0]
        data = json.load(open(samplefile, 'r'))
        args.sample = get_sample_name(samplefile)

        # extract identification information from json
        INFO = get_info(data, args.sample)

        # toggle for only identification results or not
        if not args.onlyidentification:
            RS_table, mut_table = generate_csv(data, args.sample)
            result = pd.concat([INFO, RS_table], axis=1, sort=False)
        else:
            result = INFO

        for samplefile in json_files[1:]:
            data = json.load(open(samplefile, 'r'))
            args.sample = get_sample_name(samplefile)

            # extract identification information from json
            INFO =  get_info(data, args.sample)

            # toggle for only identification results or not
            if not args.onlyidentification:
                RS_table, mut_table = generate_csv(data, args.sample)
                intermediate = pd.concat([INFO, RS_table], axis=1, sort=False)
            else:
                intermediate = INFO

            # merge all results of samples
            result = result.append(intermediate, ignore_index=True).reindex(intermediate.columns, axis=1)

            # mutation output
            if not args.onlyidentification:
                collated_mut = collated_mut.append(mut_table, ignore_index=True).reindex(mut_table.columns, axis=1)

        # toggle output files
        if args.onlyidentification:
            write_df_out(result, args.outdir, 'Identification')
        else:
            write_df_out(result, args.outdir, 'RS')
            write_df_out(collated_mut, args.outdir, 'mut')

    # for generating a csv on a single json
    elif args.json:
        if not args.sample:
            args.sample = get_sample_name(args.json)

        # read JSON result file
        data = json.load(open(args.json, 'r'))

        # extract_besthit(data, "numReads")
        # extract_besthit(data, "abundance")
        # extract_besthit(data, "numUniqueReads")
        # sys.exit()

        # extract identification information from json
        INFO = get_info(data, args.sample)

        # toggle for only identification results or not
        if not args.onlyidentification:
            RS_table, mut_table = generate_csv(data, args.sample)
            result = pd.concat([INFO, RS_table], axis=1, sort=False)
        else:
            result = INFO

        # toggle output files
        if args.onlyidentification:
            write_df_out(result, args.outdir, f"Identification_{args.sample}")
        else:
            write_df_out(result, args.outdir, f"RS_{args.sample}")
            write_df_out(mut_table, args.outdir, f"mut{args.sample}")
    else:
        print('please provide directory of myoprofiler results or single json')

if __name__ == "__main__":
    # load arguments set global workdir
    args = parse_args()
    run(args)
    print("Finished")
