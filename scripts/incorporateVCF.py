#!/usr/bin/env python

__author__ = 'Chris Mitchell'

import sys
from pythomics.templates import CustomParser
import pythomics.parsers.fasta as fasta
import pythomics.genomics.parsers as gp

description = """
This script will incorporate the variants in a given VCF file into a specified
fasta file.
"""

parser = CustomParser(description = description)
parser.add_fasta(help="The fasta file to incorporate changes into.")
parser.add_out(help="The file to write resulting fasta file to.")
parser.add_vcf()

def main():
    args = parser.parse_args()
    file_name = args.fasta
    vcf = args.vcf
    snps = args.no_snps
    dels = args.dels
    ins = args.ins
    homs = args.no_homozygous
    hets = args.heterozygous
    individual = args.individual-1
    fasta_file = fasta.FastaIterator(file_name)
    vcf_file = gp.VCFIterator( vcf )
    #store our vcf file first
    entries = {}
    to_append = 'chr' if args.append_chromosome else ''
    for info in vcf_file:
        checked = False
        valid_variant = False
        if homs:
            if info.is_homozygous()[individual]:
                if ((snps and not info.has_snp(individual=individual)) and
                    (dels and not info.has_deletion(individual=individual)) and
                    (ins and not info.has_insertion(individual=individual))):
                        checked = True
                        continue
                valid_variant = True
                try:
                    entries['%s%s' % (to_append,info.chrom)][int(info.pos)-1] = info
                except KeyError:
                    entries['%s%s' % (to_append,info.chrom)] = {int(info.pos)-1: info}
        if hets:
            if info.is_heterozygous()[individual]:
                if ((not valid_variant and not checked) and
                    (snps and not info.has_snp(individual=individual)) and
                    (dels and not info.has_deletion(individual=individual)) and
                    (ins and not info.has_insertion(individual=individual))):
                        continue
                try:
                    entries['%s%s' % (to_append,info.chrom)][int(info.pos)-1] = info
                except KeyError:
                    entries['%s%s' % (to_append,info.chrom)] = {int(info.pos)-1: info}
    with args.out as o:
        for header, sequence in fasta_file:
            d = entries.get(header, None)
            if d:
                bases = d.keys()
                bases.sort(reverse=True)
                sequence = list(sequence)
                #we go from the back of the sequence so we don't have to bother
                #with offsets if we are inserting/deleting bases as well
                for i in bases:
                    var_info = d[i]
                    ref = var_info.ref
                    alt = var_info.get_alt(individual=individual)[0]
                    # sys.stderr.write('swapping %s with %s\n' % (ref,alt))
                    sequence[i:i+len(ref)] = list(alt)
                sequence = ''.join(sequence)
            o.write('>%s\n%s\n' % (header, sequence))



if __name__ == "__main__":
    sys.exit(main())