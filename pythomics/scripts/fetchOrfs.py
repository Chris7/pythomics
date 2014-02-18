#!/usr/bin/env python

description = """
This script will accept a given nucleotide fasta file and output
found ORFs. ORFs are annotated by which stop codon they are a part
of. As in, ORF 3 is annotated as the 3rd sequence if the translated
sequence is divided by stop codons. This is prevent ambiguity with
differing minimum lengths of ORFs.
"""

import argparse
import sys
import pythomics.parsers.fasta as fasta

parser = argparse.ArgumentParser(description = description, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-f', '--fasta', nargs='?', help="The fasta file to cleave.", type=argparse.FileType('r'), default=sys.stdin)
parser.add_argument('-o', '--out', nargs='?', help="The file to write digested products to.", type=argparse.FileType('w'), default=sys.stdout)
parser.add_argument('--min', help="Minimum ORF length in amino acids.", type=int, default=50)
parser.add_argument('--both-strands', help="Search both strands for ORFs.", action='store_true', default=False)

def main():
    args = parser.parse_args()
    file_name = args.fasta
    orf_min = args.min
    fasta_file = fasta.FastaIterator(file_name)
    negative_strand = args.both_strands
    with args.out as o:
        for header, sequence in fasta_file:
            for i in xrange(3):
                strand='+'
                translation = fasta._translate(sequence[i:])
                translation = translation.split('*')
                for protein_index,protein_sequence in enumerate(translation):
                    if len(protein_sequence) >= orf_min and protein_sequence[0] == 'M':
                        o.write('>%s F:%s%d Orf:%d\n%s\n' % (header,strand,i+1,protein_index+1,protein_sequence))
                if negative_strand:
                    strand = '-'
                    translation = fasta._translate(fasta._reverse_complement(sequence)[i:])
                    for protein_index,protein_sequence in enumerate(translation):
                        if len(protein_sequence) >= orf_min and protein_sequence[0] == 'M':
                            o.write('>%s F:%s%d Orf:%d\n%s\n' % (header,strand,i+1,protein_index+1,protein_sequence))

if __name__ == "__main__":
    sys.exit(main())
