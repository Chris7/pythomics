#!/usr/bin/env python

__author__ = 'chris'

description = """
This script will lookup features from one delimited file in another delimited file, and
perform various operations on the found entries in the alternative file
"""

import argparse, sys, re, csv, copy, decimal
from pythomics.templates import CustomParser
import pythomics.proteomics.config as config
import pythomics.proteomics.digest as digest
import pythomics.parsers.fasta as fasta

parser = CustomParser(description = description)
parser.add_delimited_file(files=['-a'], delimiter=['--adelim'], cols=['--acol'], header=['--aheader'], help="This is the file to lookup values from.")
parser.add_delimited_file(files=['-b'], delimiter=['--bdelim'], cols=['--bcol'], header=['--bheader'], help="This is the file to lookup values in.")
parser.add_argument('--blookup', help='The column to take entries from in file b.', type=str, default=1)
parser.add_argument('--strict', help='For numeric operations, fail if types are incorrect (converting NA to a float for instance).', action='store_true')
parser.add_out()
parser.add_argument('--function', help='The function to apply to found entries.', choices=['concat', 'mean', 'sum', 'median', 'var', 'std'], type=str, default='concat')
parser.add_argument('--colname', help='The column name to give the new appended value. Defaults to function chosen', type=str, default='')
parser.add_argument('--aregex', help='An optional regex pattern for matching columns in file a.', type=str, default='')
parser.add_argument('--bregex', help='An optional regex pattern for matching columns in file b.', type=str, default='')

class CorrelateFunctions(object):
    strict = False

    def __init__(self, parser_args):
        if parser_args.strict:
            self.strict = True

    def process_list(self, l):
        if self.strict:
            try:
                l = [float(i) for i in l if i != '']
            except ValueError:
                sys.stderr.write('Invalid entry found in list %s.\n'%l)
                raise ValueError
        else:
            nl = []
            for i in l:
                try:
                    nl.append(float(i))
                except ValueError:
                    pass
            l = nl
        return l

    def concat(self, l):
        if not len(l):
            return 'NA'
        return ';'.join([str(i) for i in l])

    def mean(self, l):
        if not len(l):
            return 'NA'
        l = self.process_list(l)
        return sum(l)/float(len(l))

    def median(self, l):
        if not len(l):
            return 'NA'
        elif len(l) == 1:
            return l[0]
        l = self.process_list(l)
        l_sorted = sorted(l)
        mid = len(l_sorted)/2
        if len(l_sorted)%2:
            # odd
            return l_sorted[mid]
        else:
            return self.mean(l_sorted[mid-1:mid+1])

    def var(self, l):
        if not len(l):
            return 'NA'
        elif len(l) == 1:
            return 0
        l = self.process_list(l)
        sample_mean = self.mean(l)
        dev = 0
        for x in l:
            dev += (x - sample_mean)*(x - sample_mean)
        return dev/(len(l)-1)

    def std(self, l):
        if not len(l):
            return 'NA'
        elif len(l) == 1:
            return 0
        l = self.process_list(l)
        import math
        return math.sqrt(self.var(l))

    def sum(self, l):
        if not len(l):
            return 'NA'
        l = self.process_list(l)
        return sum(l)

def main():
    args = parser.parse_args()
    a_colname, b_colname, bl_colname = False, False, False
    try:
        a_column = int(args.acol)
        a_column = a_column-1 if a_column > 0 else a_column
    except ValueError:
        a_colname = True
        a_column = args.acol
    try:
        b_column = int(args.bcol)
        b_column = b_column-1 if b_column > 0 else b_column
    except ValueError:
        b_colname = True
        b_column = args.bcol
    a_file = args.a
    b_file = args.b
    a_header_lines = args.aheader
    b_header_lines = args.bheader
    a_delimiter = args.adelim
    b_delimiter = args.bdelim
    correlator = CorrelateFunctions(args)
    try:
        b_lookup = int(args.blookup)
        b_lookup = b_lookup-1 if b_lookup > 0 else b_lookup
    except ValueError:
        bl_colname = True
        b_lookup = args.blookup
    func = getattr(correlator, args.function, correlator.concat)
    b_vals = {}
    colname = args.colname if args.colname else func
    with b_file as f:
        reader = csv.reader(f, delimiter=b_delimiter)
        for line_num, entry in enumerate(reader):
            if line_num < b_header_lines:
                if b_colname:
                    for i, v in enumerate(entry):
                        if v.lower() == b_column.lower():
                            b_column = i
                            break
                if not str(b_column).isdigit():
                    sys.stderr.write('Cannot find --bcol named %s.\n'%b_column)
                    sys.stderr.write('Columns available are %s.\n'%entry)
                    return 1
                if bl_colname:
                    for i, v in enumerate(entry):
                        if v.lower() == b_lookup.lower():
                            b_lookup = i
                            break
                if not str(b_lookup).isdigit():
                    sys.stderr.write('Cannot find --blookup named %s.\n'%b_lookup)
                    sys.stderr.write('Columns available are %s.\n'%entry)
                    return 1
                continue
            else:
                try:
                    b_vals[entry[b_column]].append(entry[b_lookup])
                except KeyError:
                    b_vals[entry[b_column]] = [entry[b_lookup]]
                except IndexError:
                    pass
    with args.out as o:
        writer = csv.writer(o, delimiter=a_delimiter)
        with a_file as f:
            reader = csv.reader(f, delimiter=a_delimiter)
            for line_num, entry in enumerate(reader):
                if line_num < a_header_lines:
                    if a_colname:
                        for i, v in enumerate(entry):
                            if v.lower() == a_column.lower():
                                a_column = i
                                break
                    if not str(a_column).isdigit():
                        sys.stderr.write('Cannot find --acol named %s.\n'%a_column)
                        sys.stderr.write('Columns available are %s.\n'%entry)
                        return 1
                    entry.append(colname)
                else:
                    try:
                        lookup_val = entry[a_column]
                        entry.append(func(b_vals.get(lookup_val, [])))
                    except IndexError:
                        entry.append('')
                writer.writerow(entry)

if __name__ == "__main__":
    sys.exit(main())