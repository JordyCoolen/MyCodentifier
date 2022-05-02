#!/usr/bin/env python

import sys
import argparse
import traceback
import re
import os
import copy
import simplejson as json

INFO = 'Module for reporting resistant SNPs by mapping a VCF to a resistance database'

# EXAMPLE JSON
# {"resistance":
# "Rv0667": {
#     "gene_name": "rpoB",
#     "mutations": [
#       {
#         "hgvs_c": "c.1349c>t",
#         "hgvs_p": "p.ser450leu",
#         "eff": "missense_variant",
#         "drugs": [
#           "rifampicin"
#         ],
#         "lit": [
#           ""
#         ],
#         "dp": "43"
#       }
#     ]
# }
# "all_variants":
# "Rv0002": {
#     "gene_name": "dnaN",
#     "mutations": [ ...

# }

def parse_args():
    parser = argparse.ArgumentParser(description=INFO,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--SNPsift', '-s', nargs='?', type=argparse.FileType('r'),
                    default=sys.stdin, help="VCF file processed by 'SNPsift extractFields ANN[*]'")
    parser.add_argument('--gff', '-g', nargs='?',
                    help='GFF file to check VCF and resistance DB')
    # parser.add_argument('--resDB', '-d', nargs='?',
    #                 help='resistance db in format: gene_or_locus\thgvs_notation\tdrug\n')
    parser.add_argument('--resDB_dir', '-d', nargs='?',
                    help='dir of resistance DBs with format: gene_or_locus\thgvs_notation\tdrug\n')
    parser.add_argument('--json', '-j', nargs='?', type=argparse.FileType('r'),
                    help='Json file with meta info about the run')
    parser.add_argument('--acc', '-a', nargs='?',
                    help='Accession code to select corresponding resistance DB for organism')
    parser.add_argument('--use_old', '-o', action='store_true', default=False,
                        help='use the old locus tags for the new versions of the GFF files')
    parser.add_argument('--check_range', '-r', action='store_true', default=False,
                        help='check a range of position within a window to match TBDB and reported mutation (fix'
                             'for "shifted" SNPEff annotation')

    # parse all arguments
    args = parser.parse_args()
    return args

def run(check_range):
    try:
        # read the snpsift input excluding the typed command from stdin
        rows = snpsift.readlines()[1:]
        # parse VCF for resistance mutations
        variants = VCF_parse(rows, check_range=args.check_range)
        # write to json file
        write_results(variants)
    except Exception as e:
        exc_info = sys.exc_info()
        print(e)
        traceback.print_exception(*exc_info)
        sys.exit(1)

def get_db_path(acc, resDB):
    basename = '{}_resistance_db.tsv'
    dbpath = os.path.join(resDB, basename.format(acc))
    if os.path.exists(dbpath):
        return read_file(dbpath)
    else:
        print('no resistance db for this species, exiting..')
        exit(0)

def load_json(f):
    try:
        dict_ob = json.load(f)
        return dict_ob
    except json.JSONDecodeError:
        raise Exception('Provide a json file as run info file')

def res_db_map(locus, mut_c, mut_p, eff, cds_var, check_range=False):
    """
    do the actual mapping with the resistance DB
    :param locus: locus tag
    :param mut_c: DNA mutation notation
    :param mut_p: Protein mutation notation
    :param eff: mutation effect(e.g. synonymous)
    :return: resistance entry, in_gene entry, default: False, False
    """
    res_entry = False
    in_gene = False
    other_res_genes = [l.split('_')[0] for l in res_db['res_other']]
    # make unique mutation ids by concatenating hgvs mutation to locus
    unq_mut_c = '{}_{}'.format(locus, mut_c.lower())
    unq_mut_p = '{}_{}'.format(locus, mut_p.lower())
    del_ins = re.match(r'.*(del|ins|dup).*', mut_c)

    def get_res_entry(unq_mut, sub_rec):
        drugs = res_db[sub_rec][unq_mut]['drugs']
        lit = res_db[sub_rec][unq_mut]['lit']
        return {'hgvs_c': mut_c, 'hgvs_p': mut_p, 'eff': eff, 'drugs': drugs, 'lit': lit}

    # compare both dna mutation and protein mutation
    if locus in res_db['all_genes']:
        if unq_mut_c in res_db['res_hgvs']:
            res_entry = get_res_entry(unq_mut_c, 'res_hgvs') # compare dna mutation with DB
        elif unq_mut_p in res_db['res_hgvs']:
            res_entry = get_res_entry(unq_mut_p, 'res_hgvs') # compare protein mutation with DB
        elif del_ins and check_range and varmap_to_tbdb(mut_c, locus, del_ins.group(1), res_db['res_hgvs']):
            res_entry = varmap_to_tbdb(mut_c, locus, del_ins.group(1), res_db['res_hgvs']) # compare del/ins mutation in range

        # check if locus is associated with large deletion, frameshift or any missense mutations
        # mutation must be in CDS (in_gene: result of regex for CDS HGVS mutation)
        # why? elif locus in other_res_genes and cds_var:
        elif locus in other_res_genes and cds_var:
            # get the posion digits from the hgvs notation
            pos = re.findall(r'[0-9]+', mut_p)[0]
            # use mut_c or mut_p for pos or dependent on variant?
            #pos = re.findall(r'[0-9]+', mut_c)[0]
            # check if deletion of this gene is present in resistance DB
            if '{}_large_deletion'.format(locus) in res_db['res_other'] and eff in ['inframe_deletion',
                                                                                    'disruptive_inframe_insertion']:
                    res_entry = get_res_entry('{}_large_deletion'.format(locus), 'res_other')
            # check if frameshift of this gene is present in resistance DB
            elif '{}_frameshift'.format(locus) in res_db['res_other'] and eff in ['frameshift_variant']:
                    res_entry = get_res_entry('{}_frameshift'.format(locus), 'res_other')
            # check if missense of this gene is present in resistance DB
            elif '{}_any_missense_codon_{}'.format(locus, pos) in res_db['res_other']:
                res_entry = get_res_entry('{}_any_missense_codon_{}'.format(locus, pos), 'res_other')
        else:
            # variant is on a resistance related gene
            in_gene = {'hgvs_c': mut_c, 'hgvs_p': mut_p, 'eff': eff}
    return res_entry, in_gene

def append_to_dict(dic, subdic, locus, gene, entry):
    try:
        dic[subdic][locus]['mutations'].append(entry)
    except KeyError:
        dic[subdic][locus] = {'gene_name': gene, 'mutations': [entry]}
    return dic

def varmap_to_tbdb(mut, locus, muttype, resmuts):
    """
    Map del/ins mutations to DB with a range of offsets -4 up to +4.
    Pos might be offset because of consequence call interpretation
    :param mut: hgvs nucleotide
    :param locus: gene id
    :param muttype: type of mutation: deletion, insertion, duplicate(insertion)
    :param resmuts: dictionary with all resistance mutations
    :return:
    """

    digits = re.findall(r'\d{1,5}', mut)
    nuc_match = re.match(r'.*(dup|ins|del)([a-zA-Z]*)', mut)
    for i in range(-4, 5):
        search_digits = [str(int(digits[0])+i), str(int(digits[0])+i)]
        if len(digits) != 1:
            search_digits[i] = str(int(digits[1])+i)
        elif muttype in ['ins', 'dup']:
            search_digits[1] = str(int(digits[0])+i+1)

        # search_digits = {d: d+i for d in digits}
        if muttype == 'del':
            search_mut = 'c.{}_{}del'.format(search_digits[0], search_digits[1])
        else:
            muttype = 'ins' if muttype == 'dup' else muttype
            nucs = nuc_match.group(2).lower()
            search_mut = 'c.{}_{}{}{}'.format(search_digits[0], search_digits[1], muttype, nucs)
            # find nucleotides, append to search mut
        unq_mut = '{}_{}'.format(locus, search_mut)
        if unq_mut in resmuts:
            return get_res_entry(unq_mut, 'res_hgvs')
    return ''
    # take start - gene - vcf pos (+1 pos index correction)(+1?)

def variant_map(variants, var_info, var_dp, var_af, cds_var, check_range=False):
    """
    map the mutation to the resistance database and subcategorize mutation
    :param variants: dict with subdicts to dump mutations in different categories
    :param var_info: tuple containing: gene_id, gene_name, mutation_c, mutation_p, eff
    :return: updated variant dict
    """
    # elements of variant-gene entry
    gene_id, gene_name, mutation_c, mutation_p, eff = var_info

    # check with gff if the gene name/locus is known
    locus, gene, success = GFF_check(gene_id)
    if success:
        # do resistance DB check
        res_entry, in_gene_entry = res_db_map(locus, mutation_c, mutation_p, eff, cds_var, check_range)
        if res_entry:
            # append entry to corresponding locus key
            res_entry['dp'] = var_dp
            res_entry['af'] = var_af
            variants['resistance'][locus]['mutations'].append(res_entry)
        elif in_gene_entry:
            in_gene_entry['dp'] = var_dp
            in_gene_entry['af'] = var_af
            variants['resistance_gene'][locus]['mutations'].append(in_gene_entry)
        else:
            # add entry to the all_variants key
            variants = append_to_dict(variants, 'all_variants', locus, gene_name , \
                                      {'hgvs_c': mutation_c, 'hgvs_p': mutation_p, 'eff': eff, 'dp': var_dp, 'af': var_af})
    else:
        # gff check not successful, append to not_found dict
        variants = append_to_dict(variants, 'unknown', gene_id, gene_name , \
                                  {'hgvs_c': mutation_c, 'hgvs_p': mutation_p, 'eff': eff, 'dp': var_dp, 'af': var_af})
    return variants

def VCF_parse(vcf, check_range=False):
    """
    parse SNPEFF annotated VCF file and extract relevant info to filter on
    :param vcf: SNPEFF annotated VCF file
    :return: updated variant dictionary
    """
    # copy all resistance genes to subdirectories
    res_muts = copy.deepcopy(res_db['all_genes'])
    cds_mut = copy.deepcopy(res_db['all_genes'])
    all_drugs = copy.deepcopy(res_db['all_drugs'])
    all_vars = {}
    not_found = {}
    # combine sub dictionaries in one dictionary
    variants = {'resistance': res_muts, 'resistance_gene': cds_mut, 'all_variants': all_vars,
                'all_drugs': all_drugs, 'unknown': not_found}
    for var in vcf:
        # check if nonempty line below the header
        if not var.startswith('#') and var:
            # get ann field and depth field from input
            var_ann = var.split('\t')[:-2]
            var_dp = var.split('\t')[-2].strip()
            # get allele frequency field from input, if type is del, dp and af will be '', if type is ins, af will be 0
            if var_dp and int(var_dp) > 0:
                var_al1 = var.split('\t')[-1].strip().split(',')[0]
                var_al2 = var.split('\t')[-1].strip().split(',')[1]
                var_af = str(format(int(var_al2)/int(var_dp), '.2f'))
            else:
                var_dp = '0'
                var_af = '0'

            # for alle frequency, two??? extra tabs from snpsift
            # 0 0,59
            # 0 0,88
            # 0 0,89
            # 0 0,63
            # 0,89 0
            # 0, 108

            # 1/0 does not occur only 0/1 and 1/1

            # for every gene associated with the variant
            for gene in var_ann:
                elements = gene.split('|')
                eff = elements[1]
                gene_name = elements[3]
                gene_id = elements[4]
                mutation_c = elements[9]
                mutation_p = elements[10]
                distance = elements[14] if elements[14] else 0
                # tuple containing elements needed for mapping
                var_info = gene_id, gene_name, mutation_c, mutation_p, eff
                # pattern matching the notation for in-gene mutations, point or indels
                # downstream gene has a - sign, upstream has a * sign, so will not match
                cds_mut = re.match(r'c\.\d{1,6}[a-zA-Z_]>[a-zA-Z_]|'
                                   r'c\.(\d{1,6}_)?\d{1,6}(del|ins|dup)[ATCGatcg]?', mutation_c)
                # if the mutation is not in a gene and the distance is more than 150 bp,
                # or if the mutation startswith n (intragenic, intergenic) or empty, continue to next var
                if (not cds_mut and int(distance) > 200) or \
                   (not mutation_c and not mutation_p):
                    continue
                # above filter is passed, call resistance mapping function
                variants = variant_map(variants, var_info, var_dp, var_af, cds_mut, check_range)
                # if mutation is in-gene (matches cds_mut), don't loop upstream/downstream, continue to next mutation
                if cds_mut:
                    break
                # if eff del or insertion, loop to next snpeff entry to get absolute del position (avoid interpretation mismatch)
                # compute new relative position with info from gff file. E.g n.4326997delC matches c.460del from TBDB
                # so start (from GFF) - abs pos +1 (offset 1)

                # also, map "c.1080delA" to c.1080_1080del
    # return top dictionary with all subdicts filled with vcf info
    return variants

def GFF_check(id):
    """
    check gene names and locus tags with entries in the reference GFF
    :param id: gene name or locus tag
    :return: locus_tag, gene_name, success (boolean)
    """
    # read the lines and skip headers
    for line in gff_lines:
        if not line.startswith("#"):
            if line.split('\t')[2].strip() == 'gene':
                info = line.split('\t')[8]
                gene_id = re.search(r'ID=([0-9a-zA-Z_-]*)(;)?', info).group(1)
                if old_locus and 'old_locus_tag' in info:
                    locus = re.search(r'old_locus_tag=([0-9a-zA-Z_]*)(;)?', info).group(1)
                else:
                    locus = re.search(r'locus_tag=([0-9a-za-zA-Z_]*)(;)?', info).group(1)
                if 'gene=' in info:
                    gene_name = re.search(r'gene=([0-9a-zA-Z_]*)(;)?', info).group(1)
                else:
                    gene_name = ''
                if id == locus or id == gene_name or id == gene_id:
                    return locus, gene_name, True
    return '', '', False


def load_db(res_db_path, acc):
    """
    read resistance DB file and load into structured dict
    :param res_db_path: DB file path
    :return: structured dict with locus_mutation as keys
    """
    rows = get_db_path(acc, res_db_path)
    # remove empty strings from list
    rows = [x for x in rows if x]
    # sublist all the drugs
    all_drugs = []
    # initiate resistance dict with with key and set for all unique loci
    resdb = {'all_genes': {}, 'res_hgvs': {}, 'res_other': {}}
    # regex for all mutation types: protein, del/ins, substitution
    mut_pat = re.compile(r'(^[pn]\.-?[a-z]{3}(\d){1,6}[a-z\*]{1,3})'
                         r'|(^[cn]\.-?(\d){1,6}(_-?(\d){1,6})*(del|ins)[atcg]*)'
                         r'|(^[cn]\.-?(\d){1,5}[atcg]*>[atcg]*)', re.I)
    # add the mutations to a mutation_drug dict
    for item in rows:
        gene_db, mut, drug, lit = [s.strip() for s in item.split('\t')]
        # add the drug
        all_drugs.append(drug.lower())
        # replace r. notation with n. for compatibility with SNPEff (HGVS)
        mut = re.sub('^r.', 'n.', mut).lower()
        mut = mut.lower()
        # check GFF if locus/gene exists for the particular organism
        locus, gene, succes = GFF_check(gene_db)
        if not succes:
            # stop execution if unknown gene in resistance DB
            raise Exception('Gene in resistance DB unknown:\n'+ gene_db)
        # make unique ids by combining locus and hgvs mutation
        unq_entry = '{}_{}'.format(locus, mut)
        # there are multiple drugs for the same mutations, check if mutation is \
        # already in the dictionary, if true: append info, else make new entry
        if mut in ['large_deletion', 'frameshift'] or re.match(r'any_missense_codon_(\d)*', mut):
            if unq_entry in resdb['res_other']:
                resdb['res_other'][unq_entry]['drugs'].append(drug)
                if lit:
                    resdb['res_other'][unq_entry]['lit'].append([lit])
            else:
                resdb['res_other'][unq_entry] = {'drugs': [drug], 'lit': [lit]}
        elif mut_pat.match(mut):
            if unq_entry in resdb['res_hgvs']:
                resdb['res_hgvs'][unq_entry]['drugs'].append(drug)
                if lit:
                    resdb['res_hgvs'][unq_entry]['lit'].append([lit])
            else:
                resdb['res_hgvs'][unq_entry] = {'drugs': [drug], 'lit': [lit]}
        # mutation must match the standard annotation, else raise error
        else:
            raise Exception('Unknown mutation type in database: {}, '
                            '\nplease modify script or correct DB'.format(resdb['other'][locus]))
        # add locus to all genes to keep a separate record of resistance genes
        if locus not in resdb['all_genes']:
            resdb['all_genes'][locus] = {'gene_name': gene, 'mutations': []}
    # add a unique list of all the drugs
    resdb['all_drugs'] = list(set(all_drugs))
    return resdb


def write_results(variant_info):
    """
    write the generated json with mutation info to pipeline json
    :param variant_info: json with mutation info
    """
    meta = load_json(meta_json)
    meta['variants'] = variant_info
    with open(meta_json.name, 'w') as dbout:
        json.dump(meta, dbout)

def read_file(path):
    try:
        with open(path, 'r') as infile:
            lines = infile.readlines()
            return lines
    except Exception as e:
            exc_info = sys.exc_info()
            print(e)
            traceback.print_exception(*exc_info)
            sys.exit(1)


if __name__ == '__main__':
    args = parse_args()
    # boolean for fetching old locus tag from GFF
    old_locus = args.use_old
    # SNPsift annotated VCF
    snpsift = args.SNPsift
    # read the gff file
    gff_lines = read_file(args.gff)
    # load the resistance database
    res_db = load_db(args.resDB_dir, args.acc)
    meta_json = args.json
    # start run of program
    run(args.check_range)

