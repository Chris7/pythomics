#!/usr/bin/env python

import argparse, sys, os
import pythomics.proteomics.config as config
import pythomics.proteomics.digest as digest
import pythomics.parsers.fasta as fasta

parser = argparse.ArgumentParser()
parser.add_argument('--enzyme', help="The enzyme to cleave with.", choices=config.ENZYMES.keys(), type=str, default='trypsin')
parser.add_argument('-f', '--file', nargs='?', help="The fasta file to cleave.", type=argparse.FileType('r'), default=sys.stdin)
parser.add_argument('-o', '--out', nargs='?', help="The file to write digested products to.", type=argparse.FileType('w'), default=sys.stdout)
parser.add_argument('-t', '--type', help="The type of fasta file (default protein).", choices=['prot','nt'], type=str, default='prot')
parser.add_argument('--frame', help="If using a nucleotide file, translate in how many frames?", choices=[1,3,6], type=int)
parser.add_argument('--min', help="Minimum cleavage length", type=int, default=7)
parser.add_argument('--max', help="Maximum cleavage length", type=int, default=30)

def main():
    args = parser.parse_args()
    file_name = args.file
    enzyme_choice = args.enzyme
    digest_type = args.type
    digest_frame = args.frame
    digest_negative = False
    if digest_frame == 6:
        digest_negative = True
        digest_frame = 3
    digest_min = args.min
    digest_max = args.max
    if digest_type == 'prot' and digest_frame:
        print "Protein digestions cannot have a frame."
        return 1
    if digest_type == 'nt' and not digest_frame:
        print "Nucleotide digestions must specify the frame."
        return 1
    fasta_file = fasta.FastaIterator(file_name)
    enzyme = digest.Enzyme( enzyme=enzyme_choice )
    with args.out as o:
        if digest_type == 'nt':
            for header, sequence in fasta_file:
                for i in xrange(digest_frame):
                    strand='+'
                    for protein_index,protein_sequence in enumerate(fasta._translate(sequence[i:]).split('*')):
                        peptides = enzyme.cleave(protein_sequence, min=digest_min, max=digest_max)
                        for peptide_index,peptide in enumerate(peptides):
                            o.write('>%s F:%s%d Orf:%d Pep:%d \n%s\n' % (header,strand,i+1,protein_index+1,peptide_index+1,peptide))
                    if digest_negative:
                        strand = '-'
                        for protein_index,protein_sequence in enumerate(fasta._translate(fasta._reverse_complement(sequence)[i:]).split('*')):
                            peptides = enzyme.cleave(protein_sequence, min=digest_min, max=digest_max)
                            for peptide_index,peptide in enumerate(peptides):
                                o.write('>%s F:%s%d Orf:%d Pep:%d \n%s\n' % (header,strand,i+1,protein_index+1,peptide_index+1,peptide))
        else:
            for header, sequence in fasta_file:
                peptides = enzyme.cleave(sequence, min=digest_min, max=digest_max)
                for peptide_index,peptide in enumerate(peptides):
                    o.write('>%s Pep:%d \n%s\n' % (header,peptide_index+1,peptide))
        
    
if __name__ == "__main__":
    sys.exit(main())