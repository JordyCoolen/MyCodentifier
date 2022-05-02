#!/usr/bin/env python

import argparse
import os
import re

desc = "module for generating a SNPEff db for specific organism on the fly, based on accession code"

table_translation = \
    {1: 'Standard',
     2: 'Vertebrate_Mitochondrial',
     3: 'Yeast_Mitochondrial',
     4: 'Mold_Mitochondrial',
     5: 'Invertebrate_Mitochondrial',
     6: 'Ciliate_Nuclear',
     9: 'Echinoderm_Mitochondria',
     10: 'Euplotid_Nuclear',
     11: 'Bacterial_and_Plant_Plastid',
     12: 'Alternative_Yeast_Nuclear',
     13: 'Ascidian_Mitochondrial',
     14: 'Alternative_Flatworm_Mitochondrial',
     16: 'Chlorophycean_Mitochondrial',
     21: 'Trematode_Mitochondrial',
     22: 'Scenedesmus_obliquus_Mitochondrial',
     23: 'Thraustochytrium_Mitochondrial',
     24: None,
     25: None,
     26: None,
     27: None,
     28: None,
     29: None,
     30: None,
     31: None,
     33: None
     }

def parse_args():
    """
        Argument parser
    """
    parser = argparse.ArgumentParser(description=desc,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--acc", "-t", type=str, required=True,
                        help="accession code for the naming of organism in SNPEff index"),
    parser.add_argument("--genbank", "-g", type=str, required=True,
                        help="genbank file to retrieve the NCBI codon index")
    parser.add_argument("--template", "-p", type=str, required=True,
                        help="template file path for SNPEff DB")
    parser.add_argument("--outDir", type=str, required=True, default=None,
                        help="output directory for runtime file")
    parser.add_argument('--db', type=str, required=True,
                        help='path of snpEff index to the reference files')

    # parse all arguments
    args = parser.parse_args()

    if args.outDir is None:
        args.outDir = args.template

    return args

def write_db(temp_path, acc, outdir, org_class):
    """
    write the required lines in the SNPEff config file (genome id (referring to index) and codon table)
    :param temp_path: the template file to use for the config file
    :param acc: organism accession code, used as genome id for SNPEff
    :param outdir: directory for runtime.config output
    :param org_class: the SNPEff codon table corresponding to the NCBI table
    :return:
    """
    temp_lines = open_file(temp_path)
    genome = '{0}.genome : {0}'.format(str(acc))
    codon = '{0}.Chromosome.codonTable : {1}'.format(str(acc), org_class)
    temp_lines += '\n{}\n{}\n'.format(genome, codon)
    # check if outDir folder exist
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    new_config = os.path.join(outdir, 'runtime.config')
    with open(new_config, 'w') as new_file:
        new_file.write(temp_lines)

def get_codon_table(genbk_path):
    """
    get the SNPEff codon table corresponding to the NCBI table
    :param genbk_path: path of the genbank file to retrieve the NCBI codon table code
    :return: name of SNPEff codon table
    """
    template_lines = open_file(genbk_path)
    table_match = re.search(r'/transl_table=((\d){1,3})', template_lines)
    if table_match:
        table_idx = int(table_match.group(1))
        org_class = table_translation[table_idx]
    else:
        org_class = 'Standard'
    # table_pos = template_lines.find('/transl_table=')
    # if table_pos == '-1':
    #     raise Exception('No translation table specified in genbank')
    # codon_table = template_lines[table_pos+14]
    # org_class = table_translation.get(int(codon_table))

        # print('Warning: No suitable codon table found for corresponding NCBI codon table:'
        #       ' {} Using standard table instead'.format(codon_table))
    return org_class

def check_db_existance(acc_code, db_path):
    """
    check if index for specific organism already exists (check for accession code)
    :param acc_code: accession code
    :param db_path: path of the SNPEff index
    :return:
    """
    basepath = os.path.join(db_path, acc_code)
    if os.path.exists(basepath):
        if os.path.exists(os.path.join(basepath, 'snpEffectPredictor.bin')) and \
                os.path.exists(os.path.join(basepath, 'genes.gbk')):
            return True
    return False

def run(args):
    """
    create SNPEff config file if index does not exist then return 0, else return 1
    """
    try:
        org_class = get_codon_table(args.genbank)
        write_db(args.template, args.acc, args.outDir, org_class)
        passed = check_db_existance(args.acc, args.db)
        if passed:
            print('1')
        else:
            print('0')
    except Exception as e:
        print(e)

def open_file(path, mode='r'):
    if not os.path.exists(path):
        raise Exception('Path {} does not exist'.format(path))
    with open(path, mode) as inp:
        lines = inp.read()
    return lines


if __name__ == "__main__":
    args = parse_args()
    run(args)