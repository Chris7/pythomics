#!/usr/bin/env python

__author__ = 'Chris Mitchell'

import argparse, sys
import pythomics.parsers.fasta as fasta
import pythomics.genomics.parsers as gp

description = """
This script will incorporate the variants in a given VCF file into a specified
fasta file.
"""

parser = argparse.ArgumentParser(description = description)
parser.add_argument('--vcf', help="The VCF file to use.", type=argparse.FileType('r'), default=sys.stdin)
parser.add_argument('-f', '--file', nargs='?', help="The fasta file to incorporate changes into.", type=argparse.FileType('r'), default=sys.stdin)
parser.add_argument('-o', '--out', nargs='?', help="The file to write resulting fasta file to.", type=argparse.FileType('w'), default=sys.stdout)
parser.add_argument('--no-homozygous', help="Don't include homozygous variants (default to include)", action='store_false')
parser.add_argument('--heterozygous', help="Use heterozygous variants", action='store_true')
parser.add_argument('--no-snps', help="Don't use SNPs (default to true).", action='store_false')
parser.add_argument('--dels', help="Use Deletions.", action='store_true')
parser.add_argument('--ins', help="Use Insertions.", action='store_true')
parser.add_argument('--individual', help="The Individual to use.", type=int, default=1)
parser.add_argument('--append-chromosome', help="Should 'chr' be appended to the chromosome column?.", action='store_true')

def main():
    args = parser.parse_args()
    file_name = args.file
    vcf = args.vcf
    snps = not args.no_snps
    dels = args.dels
    ins = args.ins
    homs = not args.no_homozygous
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