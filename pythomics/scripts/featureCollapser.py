#!/usr/bin/env python

__author__ = 'chris'

description = """
This script will take a delimited file and collapse features together, such
as scan numbers. It can also be used to group peptides into longer sequences
with the --substring flag (ex: peptides LNGERPEPTIDE and ERPEPT will be merged
into LNGERPEPTIDE).
"""

import argparse, sys, re, csv, copy, decimal
from pythomics.templates import CustomParser
import pythomics.proteomics.config as config
import pythomics.proteomics.digest as digest
import pythomics.parsers.fasta as fasta

parser = CustomParser(description = description)
parser.add_delimited_file()
parser.add_out()
parser.add_argument('--substring', help='If set, merge features by partial matches (such as collapsing peptides into larger peptides)', action='store_true', default=False)
parser.add_argument('--merge-columns', help="If set, columns of merged peptides will be combined.", action='store_true', default=False)
parser.add_argument('--merge-delimiter', help='The delimiter for column merges.', type=str, default=';')
parser.add_argument('--case-sensitive', help="Treat peptides as case-sensitive (ie separate modified peptides)", action='store_true', default=False)

def main():
    args = parser.parse_args()
    peptide_column = args.col-1
    tsv_file = args.tsv
    header_lines = args.header
    delimiter = args.delimiter
    peptide_join = args.substring
    col_delimiter = args.merge_delimiter
    merge_columns = args.merge_columns
    case_sens = args.case_sensitive
    peptide_history = {}
    headers = []
    with tsv_file as f:
        reader = csv.reader(f, delimiter=delimiter)
        for line_num, entry in enumerate(reader):
            if line_num < header_lines:
                headers.append(entry)
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
            if merge_columns:
                if isinstance(peptide_history[peptide], list):
                    peptide_history[matching_peptide] += peptide_history[peptide]
                else:
                    peptide_history[matching_peptide].append(peptide_history[peptide])
            del peptide_history[peptide]
    with args.out as o:
        writer = csv.writer(o, delimiter=delimiter)
        for i in headers:
            writer.writerow(i)
        for index, (peptide, entries) in enumerate(peptide_history.iteritems()):
            if not index%1000:
                sys.stderr.write('%d of %d complete.\n' % (index, len(peptide_history)))
            if merge_columns:
                entry = []
                for peptide_index, peptide_info in enumerate(entries):
                    for i, v in enumerate(peptide_info):
                        if peptide_index == 0:
                            entry.append([v])
                        else:
                            entry[i].append(v)
                entry_string = [col_delimiter.join(v) if i != peptide_column else v[0] for i, v in enumerate(entry)]
            else:
                entry_string = entries[0]

            writer.writerow(entry_string)

if __name__ == "__main__":
    sys.exit(main())