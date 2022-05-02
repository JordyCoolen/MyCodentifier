#!/usr/bin/env python

# create html template
# .json input file
# read .json
# put all results in pdf format
# output location of pdf

######
INFO = "Convert JSON results to PDF report"
__version__ = 0.11
######

import os
import argparse
import storage.storage as storage
import pandas as pd

def parse_args():
    """
        Argument parser
    """
    parser = argparse.ArgumentParser(description=INFO, \
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--JSON", type=str, required=False,
                        help="full path to JSON file", default=None)
    parser.add_argument("--sample", type=str, required=True,
                        help="name of sequence sample"),
    parser.add_argument("-o", "--outputDir", type=str, required=False,
                        help="full path of output folder", default=os.path.abspath("./"))
    parser.add_argument("-v", "--version", action="version",
                    version="%(prog)s {version}".format(version=__version__))

    # parse all arguments
    args = parser.parse_args()

    return args

def check_contents(data):
    """
        code to check if expected keys are present in JSON
        :param JSON: JSON object containing all the results
        :return:
    """

    exp_keys = {'fastq': [],
                'centrifuge': [],
                'contaminants': [],
                'KMA': [],
                'WGS_typing': [],
                'variants': ['resistance',
                             'resistance_gene',
                             'all_variants',
                             'unknown',
                             'all_drugs'],
                'BAM_QC': ['avg_dp',
                           'low_cov_frac',
                           'dp_thres',
                           'cov_thres',
                           'qc_passed']}
    try:
        for k in exp_keys:
            data[k]
            for s in exp_keys[k]:
                data[k][s]
    except KeyError as e:
        key = e.args[0]
        print('key {} not present in json, skip report generator'.format(key))
        exit(0)


def fill_html(JSON, sample, outputDir):
    '''
        Code to fill in the placeholders in the html
        and generate a html and pdf
        
        :params JSON: JSON object containing all the results
        :params outputDir: directory to store the results
        
        :out pdf: pdf report of the results
        :out html: html report of the results
    '''
    
    import matplotlib
    matplotlib.use('Agg')
    from weasyprint import HTML
    from jinja2 import Environment, FileSystemLoader
    
    print('Start Filling')
    
    localdir = os.path.dirname(os.path.realpath(__file__))
    
    # create and render html file with tables
    env = Environment(loader=FileSystemLoader(localdir))
    template = env.get_template('report/final_report_template.html')
    
    # location of logo
    logo = os.path.join(localdir, "report/logo.png")
    logo = logo.replace(' ','%20')

    # load contents from JSON
    data = JSON.data
    # check if expected keys are present
    check_contents(data)
    
    # centrifuge
    CF_df = pd.DataFrame.from_dict(data['centrifuge'], orient='index')
    CF_df = CF_df.drop(['metric', 'data', 'barplot'])
    
    # centrifuge data
    CFd_df = pd.DataFrame.from_dict(data['centrifuge']['data'], orient='index')
    CFd_df = CFd_df.drop(['taxRank'])
    CFd_df = CFd_df.transpose()

    #BAM QC
    BAM_QC = pd.DataFrame.from_dict(data['BAM_QC'], orient='index')

    # contaminants
    CON_df = pd.DataFrame.from_dict(data['contaminants']['data'], orient='index')
    CON_df = CON_df.drop(['taxRank'])
    CON_df = CON_df.transpose()
    
    # Identification methods
    if "SNP-IT" in data.keys():
        SNPIT = pd.DataFrame.from_dict(data['SNP-IT'], orient='index')
        SNPIT = SNPIT.drop(['sample_name'])
    else:
        SNPIT = pd.DataFrame(['NO DATA'])

    # HSP65 identification
    HSP = pd.DataFrame.from_dict(data['KMA'], orient='index')
    HSP = HSP.drop(['data'])

    # WGS typing
    WGS_IDEN = pd.DataFrame.from_dict(data['WGS_typing']['data'], orient='index')
    WGS_IDEN = WGS_IDEN.drop(['taxRank'])
    WGS_IDEN = WGS_IDEN.drop(['taxID'])
    WGS_IDEN = WGS_IDEN.transpose()

    # get detected resistances
    RES = data['variants']['resistance']
    
    # create pandas dataframe of resistance mutation
    RES_df = get_resistance_mutations(RES)
    
    # obtain pandas dataframe of target gene mutations
    RGENE = data['variants']['resistance_gene']
    RGENE_df = get_resistance_gene(RGENE)
    
    # create susceptibility table
    SUS = susceptibility(RES_df, data)
    
    # fill html
    template_vars = {
                    # pretty things
                    "logo": logo,
                    "version": __version__,
                    
                    # general info
                    "strainname": data['id'],
                    "runID": data['fastq']['runID'],
                    "barcode": data['fastq']['barcode'],
                    
                    # data tables
                    "CENTRIFUGE":CF_df.to_html(index=True, header=False),
                    "CENTRIFUGEDATA":CFd_df.to_html(index=False, header=True),
                    "BAMQC":BAM_QC.to_html(header=False),
                    "CONTAMINANTS":CON_df.to_html(index=False, header=True),           
                    "SNPIT":SNPIT.to_html(index=True, header=False),
                    "HSP": HSP.to_html(index=True, header=False),
                    "WGS_IDEN": WGS_IDEN.to_html(index=False, header=True),
                    "SUSCEPTIBILITY":SUS.to_html(index=False),
                    "RESISTANCE":RES_df.to_html(index=False, header=True),
                    "RGENE":RGENE_df.to_html(index=False, header=True),
                    }
    
    # output pdf
    outfile = os.path.join(outputDir,'{}.pdf'.format(sample))
    
    # render html and write
    html_out = template.render(template_vars)
    with open(os.path.join(outputDir, '{}.html'.format(sample)), 'w') as html_file:
        html_file.write(html_out)
           
    # save html as pdf to disc
    HTML(string=html_out, base_url=__file__).write_pdf(outfile, stylesheets=[os.path.join(localdir, 'report/style.css')])

def susceptibility(df, data):
    '''
        Get all drugs in resistance database
        and compare with found resistance.
        Make susceptibility table.
        
        :params df: dataframe obtained from get_resistance_mutations
        :params data: full JSON data file
        
        :return df: dataframe containing
                    drugs susceptibility #mutations
    '''
    
    # create pandas dataframe drug summary
    SUM = df['drug'].value_counts().reset_index()
    SUM.columns = ['drug','#mutations']
    
    table = []
    # compare found resistance to all_drugs list
    for d in data['variants']['all_drugs']:
        if d in SUM['drug'].to_list():
            i = SUM['drug'].to_list().index(d)
            table.append([d, 'R', SUM['#mutations'][i]])
        else:
            table.append([d, 'S', 0])
    
    df = pd.DataFrame(table).sort_values([1,0])
    df.columns = ['drugs','Susceptibility','#mutations']
    
    return df
    
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
            af = data[locus]['mutations'][i]['af']
            for d in data[locus]['mutations'][i]['drugs']:
                row = [locus, gene, dp, af, nt, aa, eff, d]
                table.append(row)
    
    # generate pandas table
    df = pd.DataFrame(table, columns=['locus', 'gene', 'depth', 'freq', 'nt', 'aa', 'eff', 'drug'])
    df = df.sort_values(['drug'])
    
    return df

def get_resistance_gene(data):
    '''
        parse JSON data to pandas dataframe
        
        examples of data:
        data['variants']['resistance_gene'] 

        
        :params data:   data['variants']['resistance_gene']
                        all unknown mutations located on genes that
                        have known mutations causing resistance.
        
        :return df: [locus, gene, nt, aa, eff, dp] in pandas dataframe
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
            af = data[locus]['mutations'][i]['af']
            row = [locus, gene, dp, af, nt, aa, eff]
            table.append(row)
    
    # generate pandas table
    df = pd.DataFrame(table, columns=['locus', 'gene', 'depth', 'freq', 'nt', 'aa', 'eff'])
    df = df.sort_values(['locus'])
    
    return df

def run(args):
    '''
        General run code
    '''
    
    # read JSON result file
    JSON = storage.JSON()
    JSON.open(args.JSON)
    
    # fill html file
    fill_html(JSON, args.sample, args.outputDir)
    
if __name__ == "__main__":
    # load arguments set global workdir
    args = parse_args()
    run(args)
    print("Finished")
    
