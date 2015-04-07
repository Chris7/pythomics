#!/usr/bin/env python

description = """
This script will accept a given nucleotide fasta file and output
found ORFs. ORFs are annotated by which stop codon they are a part
of. As in, ORF 3 is annotated as the 3rd sequence if the translated
sequence is divided by stop codons. This is prevent ambiguity with
differing minimum lengths of ORFs.
"""

from pythomics.templates import CustomParser
import sys, argparse
import pythomics.parsers.fasta as fasta

parser = CustomParser(description = description, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_fasta()
parser.add_out()
parser.add_argument('--min', help="Minimum ORF length in amino acids.", type=int, default=50)
parser.add_argument('--both-strands', help="Search both strands for ORFs.", action='store_true')
parser.add_argument('--no-met-start', help="Output ORFs starting with amino acids other than MET", action='store_true')
parser.add_argument('--from-met', help="Truncate leading amino acids up to MET", action='store_true')
parser.add_argument('--from-met-keep', help="Truncate leading amino acids up to MET, but keep the untruncated version as well.", action='store_true')

def main():
    args = parser.parse_args()
    file_name = args.fasta
    orf_min = args.min
    fasta_file = fasta.FastaIterator(file_name)
    negative_strand = args.both_strands
    no_met = args.no_met_start
    from_met = args.from_met
    from_met_keep = args.from_met_keep
    with args.out as o:
        for header, sequence in fasta_file:
            for i in xrange(3):
                strand='+'
                translation = fasta._translate(sequence[i:])
                translation = translation.split('*')
                for protein_index,protein_sequence in enumerate(translation):
                    protein_sequences = [protein_sequence]
                    if from_met:
                        if not from_met_keep:
                            protein_sequences = []
                        pos = protein_sequence.find('M')
                        if pos == -1:
                            if not from_met_keep:
                                continue
                        protein_sequences.append(protein_sequence[pos:])
                    for protein_sequence in protein_sequences:
                        if len(protein_sequence) >= orf_min and (no_met or protein_sequence[0] == 'M'):
                            o.write('>%s F:%s%d Orf:%d\n%s\n' % (header,strand,i+1,protein_index+1,protein_sequence))
                if negative_strand:
                    strand = '-'
                    translation = fasta._translate(fasta._reverse_complement(sequence)[i:])
                    for protein_index,protein_sequence in enumerate(translation):
                        protein_sequences = [protein_sequence]
                        if from_met:
                            if not from_met_keep:
                                protein_sequences = []
                            pos = protein_sequence.find('M')
                            if pos == -1:
                                if not from_met_keep:
                                    continue
                            protein_sequences.append(protein_sequence[pos:])
                        for protein_sequence in protein_sequences:
                            if len(protein_sequence) >= orf_min and (no_met or protein_sequence[0] == 'M'):
                                o.write('>%s F:%s%d Orf:%d\n%s\n' % (header,strand,i+1,protein_index+1,protein_sequence))

if __name__ == "__main__":
    sys.exit(main())
