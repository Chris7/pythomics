#!/usr/bin/env python

import argparse, sys, itertools
import pythomics.proteomics.config as config
import pythomics.proteomics.digest as digest
import pythomics.parsers.fasta as fasta

parser = argparse.ArgumentParser()
parser.add_argument('--enzyme', help="The enzyme to cleave with.", choices=config.ENZYMES.keys(), type=str, default='trypsin')
parser.add_argument('-f', '--file', nargs='?', help="The fasta file to cleave.", type=argparse.FileType('r'), default=sys.stdin)
parser.add_argument('-o', '--out', nargs='?', help="The file to write digested products to.", type=argparse.FileType('w'), default=sys.stdout)
parser.add_argument('-t', '--type', help="The type of fasta file (default protein).", choices=['prot','nt'], type=str, default='prot')
parser.add_argument('--frame', help="If using a nucleotide file, translate in how many frames?", choices=[1,3,6], type=int)
parser.add_argument('--genome', help="Are we translating a genome? This will keep chromosome positions in the header.", choices=[True,False], type=bool)
parser.add_argument('--min', help="Minimum cleavage length", type=int, default=7)
parser.add_argument('--max', help="Maximum cleavage length", type=int, default=30)
parser.add_argument('--unique', help="Only return unique peptides per cleavage", type=bool, default=False, choices=[True,False])

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
    genome = args.genome
    unique_digest = args.unique
    #if we're splitting a genome
    if genome:
        import re
        regex = re.compile(r'([\*])')
        digest_type = 'nt'
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
                if genome:
                    slen = len(sequence)
                for i in xrange(digest_frame):
                    strand='+'
                    translation = fasta._translate(sequence[i:])
                    if genome:
                        position = i+1
                        translation = [j for j in regex.split(translation)]
                        translation = [''.join(j) for j in itertools.izip_longest(translation[0::2],translation[1::2],fillvalue='')]
                    else:
                        translation = translation.split('*')
                    for protein_index,protein_sequence in enumerate(translation):
                        if genome:
                            peptides = enzyme.cleave(protein_sequence, min=0, max=99999, unique=unique_digest)
                        else:
                            peptides = enzyme.cleave(protein_sequence, min=digest_min, max=digest_max, unique=unique_digest)
                        for peptide_index,peptide in enumerate(peptides):
                            if genome:
                                if len(peptide)>=digest_min:
                                    if peptide.endswith('*'):
                                        o.write('>%s F:%s%d Start:%d End:%d \n%s\n' % (header,strand,i+1,position,position+len(peptide)*3-1,peptide[:-1]))
                                    else:
                                        o.write('>%s F:%s%d Start:%d End:%d \n%s\n' % (header,strand,i+1,position,position+len(peptide)*3-1,peptide))
                                position+=len(peptide)*3
                            else:
                                o.write('>%s F:%s%d Orf:%d Pep:%d \n%s\n' % (header,strand,i+1,protein_index+1,peptide_index+1,peptide))
                    if digest_negative:
                        strand = '-'
                        translation = fasta._translate(fasta._reverse_complement(sequence)[i:])
                        if genome:
                            position = slen-i
                            translation = [j for j in regex.split(translation)]
                            translation = [''.join(j) for j in itertools.izip_longest(translation[0::2],translation[1::2],fillvalue='')]
                        else:
                            translation = translation.split('*')
                        for protein_index,protein_sequence in enumerate(translation):
                            if genome:
                                peptides = enzyme.cleave(protein_sequence, min=0, max=999999, unique=unique_digest)
                            else:
                                peptides = enzyme.cleave(protein_sequence, min=digest_min, max=digest_max, unique=unique_digest)
                            for peptide_index,peptide in enumerate(peptides):
                                if genome:
                                    if len(peptide)>=digest_min:
                                        if peptide.endswith('*'):
                                            o.write('>%s F:%s%d Start:%d End:%d \n%s\n' % (header,strand,i+1,position-len(peptide)*3+1,position,peptide[:-1]))
                                        else:
                                            o.write('>%s F:%s%d Start:%d End:%d \n%s\n' % (header,strand,i+1,position-len(peptide)*3+1,position,peptide))
                                    position-=(len(peptide)*3)
                                else:
                                    o.write('>%s F:%s%d Orf:%d Pep:%d \n%s\n' % (header,strand,i+1,protein_index+1,peptide_index+1,peptide))
        else:
            for header, sequence in fasta_file:
                peptides = enzyme.cleave(sequence, min=digest_min, max=digest_max)
                for peptide_index,peptide in enumerate(peptides):
                    o.write('>%s Pep:%d \n%s\n' % (header,peptide_index+1,peptide))
        
    
if __name__ == "__main__":
    sys.exit(main())