#!/usr/bin/env python

__author__ = 'chris'

description = """
This script will annotate a tab delimited text file with which protein headers
correspond to peptides present.
"""

import argparse, sys, itertools, re, csv
import pythomics.proteomics.config as config
import pythomics.proteomics.digest as digest
import pythomics.parsers.fasta as fasta

parser = argparse.ArgumentParser(description = description)
parser.add_argument('-f', '--fasta', nargs='?', help="The fasta file to cleave.", type=argparse.FileType('r'), default=sys.stdin)
parser.add_argument('-o', '--out', nargs='?', help="The file to write digested products to.", type=argparse.FileType('w'), default=sys.stdout)
parser.add_argument('-t', '--tsv', help="The tab separated file.", type=argparse.FileType('r'), required=True)
parser.add_argument('-d', '--delimiter', help="The delimiter for fields.", type=str, default='\t')
parser.add_argument('-c', '--col', help="The column with peptides (default: 1).", type=int, default=1)
parser.add_argument('--header', help="The number of headers lines (default: 1).", type=int, default=1)
parser.add_argument('-r', '--regex', help="A perl regular expression determining which parts of the header to capture.", type=str)

def main():
    args = parser.parse_args()
    fasta_file = fasta.FastaIterator(args.fasta)
    peptide_column = args.col-1
    tsv_file = args.tsv
    header_lines = args.header
    delimiter = args.delimiter
    sys.stderr.write("Reading in Fasta file.\n")
    fasta_headers, protein_sequences = zip(*[(header,sequence) for header,sequence in fasta_file])
    #replace headers with parsed ones
    if args.regex:
        regex = re.compile(args.regex)
        fasta_headers = [regex.search(header) for header in fasta_headers]
        protein_sequences = [protein_sequences[i] for i,v in enumerate(fasta_headers) if v]
        sys.stderr.write("%d header sequences did not match regex and have been discarded." % (len(fasta_headers)-len(protein_sequences)))
        fasta_headers = [' '.join(i.groups()) for i in fasta_headers if i]
    protein_sequences = '\n'.join(protein_sequences)
    with args.out as o:
        writer = csv.writer(o, delimiter=delimiter)
        with tsv_file as f:
            reader = csv.reader(f, delimiter=delimiter)
            for line_num, entry in enumerate(reader):
                if line_num < header_lines:
                    entry.append('Proteins')
                else:
                    peptide = entry[peptide_column].upper()
                    sites = [match.start() for match in re.finditer(peptide, protein_sequences)]
                    matches = ';'.join([fasta_headers[protein_sequences.count('\n', 0, match_position)] for match_position in sites])
                    entry.append(matches)
                writer.writerow(entry)

if __name__ == "__main__":
    sys.exit(main())