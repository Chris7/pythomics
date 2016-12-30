#!/usr/bin/env python

__author__ = 'chris'

description = """
This script will take a delimited file and collapse features together, optionally
creating summary statistics for them.

For instance, gene ids can be selected and their FPKM/iBAQ values combined.
Also, features can be can be grouped into longer sequences
with the --substring flag (ex: peptides LNGERPEPTIDE and ERPEPT will be merged
into LNGERPEPTIDE).
"""

import argparse, sys, re, csv, copy, decimal
from pythomics.templates import CustomParser
import pythomics.proteomics.config as config
import pythomics.proteomics.digest as digest
import pythomics.parsers.fasta as fasta
from pythomics.utils import ColumnFunctions

parser = CustomParser(description = description)
parser.add_delimited_file(cols=['--group-on'])
parser.add_out()
parser.add_argument('--substring', help='If set, merge features by partial matches (such as collapsing peptides into larger peptides)', action='store_true')
parser.add_column_function('--summary-col', col_help="The function to apply to grouped entries in modification columns.")
parser.add_argument('--summary-col-delimiter', help="If the summary column has a delimiter, such as a ; for multiple proteins.")
parser.add_argument('--strict', help='For numeric operations, fail if types are incorrect (converting NA to a float for instance).', action='store_true')
parser.add_argument('--merge', help='Merge together identical entries.', action='store_true')
# parser.add_argument('--merge-columns', help="If set, columns of merged peptides will be combined.", action='store_true')
# parser.add_argument('--merge-delimiter', help='The delimiter for column merges.', type=str, default=';')
parser.add_argument('--case-sensitive', help="Treat peptides as case-sensitive (ie separate modified peptides)", action='store_true')

def main():
    args = parser.parse_args()
    peptide_colname = False
    try:
        peptide_column = int(args.group_on)
        peptide_column = peptide_column-1 if peptide_column > 0 else peptide_column
    except ValueError:
        peptide_colname = True
        peptide_column = args.group_on
    tsv_file = args.tsv
    header_lines = args.header
    delimiters = ''.join(list(set([',','\t',args.delimiter])))
    peptide_join = args.substring
    col_func = ColumnFunctions(args)
    try:
        mod_col = int(args.summary_col)-1 if args.summary_col else False
    except ValueError:
        mod_col = None
    mod_col_func = getattr(col_func, args.summary_col_func, col_func.concat)
    summary_col_delim = args.summary_col_delimiter
    merge = args.merge
    # col_delimiter = args.merge_delimiter
    # merge_columns = args.merge_columns
    case_sens = args.case_sensitive
    peptide_history = {}
    headers = []
    with tsv_file as f:
        dialect = csv.Sniffer().sniff(f.read(1024*16), delimiters=delimiters)
        f.seek(0)
        reader = csv.reader(f, dialect)
        for line_num, entry in enumerate(reader):
            if line_num < header_lines:
                headers.append(entry)
                if peptide_colname:
                    for i, v in enumerate(entry):
                        if v.lower() == peptide_column.lower():
                            peptide_column = i
                            break
                if mod_col is None:
                    for i, v in enumerate(entry):
                        if v.lower() == args.summary_col.lower():
                            mod_col = i
                            break
            else:
                peptide = entry[peptide_column]
                if not case_sens:
                    peptide = peptide.upper()
                try:
                    peptide_history[peptide].append(entry)
                except KeyError:
                    peptide_history[peptide] = [entry]
    #order the peptides by length and make a single string of them for searching
    peptide_keys = peptide_history.keys()
    peptide_keys.sort(key=lambda x: len(x), reverse=True)
    peptide_string = '\n'.join(peptide_keys)
    #now we go from the back of the list (smallest, and merge them into the largest entry)
    if peptide_join:
        lpos = len(peptide_string)
    else:
        lpos = len(peptide_keys)
    for peptide in reversed(peptide_keys):
        if peptide_join:
            lpos-=len(peptide)+1
            first_match = peptide_string.find(peptide, 0, lpos)
        else:
            try:
                lpos -= 1
                first_match = peptide_keys[:lpos].index(peptide)
            except ValueError:
                first_match = 0
        if first_match > 0:
            if peptide_join:
                matching_peptide = peptide_keys[peptide_string.count('\n', 0, first_match)]
            else:
                matching_peptide = peptide_keys[first_match]
            if mod_col:
                if isinstance(peptide_history[peptide], list):
                    peptide_history[matching_peptide] += peptide_history[peptide]
                else:
                    peptide_history[matching_peptide].append(peptide_history[peptide])
            del peptide_history[peptide]
    with args.out as o:
        writer = csv.writer(o, dialect=dialect)
        for i in headers:
            writer.writerow(i)
        for index, (peptide, entries) in enumerate(peptide_history.iteritems()):
            if not index%1000:
                sys.stderr.write('%d of %d complete.\n' % (index, len(peptide_history)))
            if mod_col:
                entry = []
                for peptide_index, peptide_info in enumerate(entries):
                    for i, v in enumerate(peptide_info):
                        if peptide_index == 0:
                            entry.append([v])
                        else:
                            entry[i].append(v)

                entry_string = []
                for i, v in enumerate(entry):
                    if i == mod_col:
                        if summary_col_delim:
                            l = [j for i in v for j in i.split(summary_col_delim)]
                        else:
                            l = v
                        entry_string.append(mod_col_func(set(l) if merge else l))
                    else:
                        entry_string.append(v[0])
            else:
                entry_string = entries[0]

            writer.writerow(entry_string)

if __name__ == "__main__":
    sys.exit(main())