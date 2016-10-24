#!/usr/bin/env python
from __future__ import division, absolute_import

__author__ = 'chris'

description = """

"""

import sys
import os
import operator
import argparse

from pythomics.templates import CustomParser

parser = CustomParser(description=description)
group = parser.add_argument_group('Protein Inference File')
group.add_argument('--inference', help="The protein inference file (your peptide file with gene/protein annotations). For multiple files, separate by spaces (must be in same order as mods).", nargs='+', type=argparse.FileType('r'), required=True)
group.add_argument('--gene', help="The Gene column name", type=str, default='Gene')
group.add_argument('--protein', help="The Protein column name", type=str, default='Protein')
group.add_argument('--peptide', help="The Peptide column name", type=str, default='Peptide')
group.add_argument('--quant', help="The name of quantification columns (such as Heavy/Light). Separate multiple columns by spaces", nargs='+', default=['Heavy/Light'])
mods = parser.add_argument_group('Modification File')
mods.add_argument('--mods', help="The modifications file (the file with sites, peptides). For multiple files, separate by spaces (must be in same order as inference).", nargs='+', type=argparse.FileType('r'), required=True)
mods.add_argument('--site-protein', help="The mod file protein column name", type=str, default='Protein')

parser.add_argument('--no-log2', help='Do not log2 normalize quantification values.', action='store_true')
parser.add_argument('--no-median', help='Do not normalize quantification values by the median of the experiment.', action='store_true')
parser.add_argument('--wp', help="The whole proteome inference file, if it exists. For multiple replicates, separate by spaces.", nargs='+', type=argparse.FileType('r'))
parser.add_argument('--non-mod-norm', help='Normalize the data by the non-modified peptides.', action='store_true')

parser.add_argument('--site-file', help='The output path for the file with sumamries at the site level.', default=sys.stdout, type=argparse.FileType('wb'))
parser.add_argument('--peptide-file', help='The output path for the file with sumamries at the site and peptide level.', default=sys.stdout, type=argparse.FileType('wb'))


def main():
    args = parser.parse_args()
    inference_files = args.inference
    mod_files = args.mods
    wp_files = args.wp if args.wp else []
    quant_cols = args.quant
    gene_col = args.gene
    prot_col = args.protein
    pep_col = args.peptide
    log_norm = not args.no_log2
    med_norm = not args.no_median
    site_col = args.site_protein
    non_mod_norm = args.non_mod_norm

    import pandas as pd
    import numpy as np
    # do wp things
    wp_medians = []
    wp_reps = []
    sample_names = []
    for wp_index, wp_rep in enumerate(wp_files):
        rep = pd.read_table(wp_rep.name)
        sample_names.append(os.path.split(wp_rep.name)[1])
        wp = pd.DataFrame(columns=quant_cols+['WP_Mean'])
        wp_reps.append(wp)
        for quant_col in quant_cols:
            qrep = rep[pd.isnull(rep[quant_col])==False]
            rep[quant_col] = np.log2(qrep[quant_col]) if log_norm else qrep[quant_col]
            #normalize each to their median
            if med_norm:
                median = rep[quant_col].median()
                rep[quant_col] = rep[quant_col]-median
                wp_medians.append(median)
            wp[(wp_index, quant_col)] = rep.groupby(gene_col)[quant_col].median()
    wp_summary = pd.DataFrame(columns=quant_cols+['Mean'], index=set([i for wp in wp_reps for i in wp.index]))
    for quant_col in quant_cols:
        wp_summary[quant_col] = pd.concat([wp[(wp_index, quant_col)] for wp_index, wp in enumerate(wp_reps)], axis=1).median(axis=1) if wp_reps else 0
    wp_medians = float(sum(wp_medians))/len(wp_medians) if wp_medians else 0

    # ptm stuff
    import copy
    from collections import OrderedDict

    def process_file(mods=None, all_data=None, col_names=None, file_quant_cols=None, wp_median=None, df=None, df_collapsed=None):
        for index, row in mods.iterrows():
            site = row['Site'].lower()
            protein = row[site_col]
            peptides = list(set(row['Peptide'].split(';')))
            for peptide in peptides:
                genes = all_data[(all_data[pep_col] == peptide)==True][gene_col].str.split(';').fillna('NA')
                gene = ';'.join(list(OrderedDict.fromkeys([j for i in genes for j in i])))
                matches = all_data[all_data[pep_col] == peptide]
                df_col_index = (gene, protein, site)
                for quant_col in file_quant_cols:
                    med = ((np.log2(matches[quant_col]) if log_norm else matches[quant_col])-wp_median).median()
                    for col_name in col_names[quant_col]:
                        df_index = (gene, protein, site, peptide)
                        if df_index in df and col_name in df[df_index]:
                            continue
                        try:
                            df[df_index][col_name] = med
                        except KeyError:
                            df[df_index] = {col_name: med}
                        try:
                            df_collapsed[df_col_index][col_name].add(med)
                        except KeyError:
                            if df_col_index in df_collapsed:
                                df_collapsed[df_col_index][col_name] = set([med])
                            else:
                                df_collapsed[df_col_index] = {col_name: set([med])}

    df = {}
    df_collapsed = {}
    for ptm_file, mod_file in zip(inference_files, mod_files):
        p1 = pd.read_table(ptm_file.name)
        file_name = os.path.split(ptm_file.name)[1]
        p1_mods = pd.read_table(mod_file.name)
        colnames = {}
        for quant_col in quant_cols:
            colname = '{}-{}'.format(file_name,quant_col)
            try:
                colnames[quant_col].append(colname)
            except KeyError:
                colnames[quant_col] = [colname]
        process_file(mods=p1_mods, all_data=p1, col_names=colnames, file_quant_cols=quant_cols, wp_median=wp_medians, df=df, df_collapsed=df_collapsed)
        if non_mod_norm:
            a = p1[p1['Modifications'].str.contains('Phospho')==False].groupby((gene_col, prot_col, pep_col)).median()
            a = a.groupby(level=[0,1]).median()
            a = a.groupby(level=[0]).median()
            wp_summary = a.loc[:, quant_cols]

    # do a median on the collapsed items
    for i,v in df_collapsed.iteritems():
        for col in v.keys():
            collapsed_items = [k for k in v[col] if not pd.isnull(k)] if isinstance(v[col], set) else v[col]
            v[col] = np.median(collapsed_items) if collapsed_items else np.nan

    df = pd.DataFrame(df).T
    df_collapsed = pd.DataFrame(df_collapsed).T
    for quant_col in quant_cols:
        df['{}_Medians'.format(quant_col)] = df.loc[:, colnames[quant_col]].median(axis=1)
        df_collapsed['{}_Medians'.format(quant_col)] = df_collapsed.loc[:, colnames[quant_col]].median(axis=1)
        if not wp_summary.empty:
            df['WP_{}'.format(quant_col)] = wp_summary[quant_col]
    df.index.names = [gene_col, prot_col, 'Site', pep_col]
    df_collapsed.index.names = [gene_col, prot_col, 'Site']

    for gene, row in wp_summary.iterrows():
        try:
            df.loc[gene]
        except KeyError:
            continue
        for quant_col in quant_cols:
            df.loc[gene, 'WP_{}'.format(quant_col)] = row[quant_col]
            df_collapsed.loc[gene, 'WP_{}'.format(quant_col)] = row[quant_col]

    df.to_csv(args.peptide_file, sep='\t')
    df_collapsed.to_csv(args.site_file, sep='\t')

if __name__ == "__main__":
    sys.exit(main())
