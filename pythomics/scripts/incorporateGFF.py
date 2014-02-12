#!/usr/bin/env python

__author__ = 'Chris Mitchell'

import argparse
import sys
import operator
import pythomics.parsers.fasta as fasta
import pythomics.genomics.parsers as gp

description = """
This script will incorporate the a given GFF file into a specified
fasta file. It can also incorporate variants given in a VCF file
while generating this fasta file (eventually).
"""

parser = argparse.ArgumentParser(description = description)
parser.add_argument('-f', '--fasta', nargs='?', help="The fasta file to reference.", type=argparse.FileType('r'), default=sys.stdin)
parser.add_argument('-o', '--out', nargs='?', help="The file to write resulting fasta file to.", type=argparse.FileType('w'), default=sys.stdout)
parser.add_argument_group('GFF file related options')
parser.add_argument('--gff', help="The GFF file to use.", type=argparse.FileType('r'))
parser.add_argument('--group-on', help="The key to group entries together by (such as transcript_id)", type=str, default='ID')
parser.add_argument('--cufflinks', help="If the gff file is in the standard cufflinks output", action='store_true')
# Not In Yet
# parser.add_argument_group('VCF file related options')
# parser.add_argument('--vcf', help="The VCF file to use.", type=argparse.FileType('r'))
# parser.add_argument('--no-homozygous', help="Don't include homozygous variants (default to include)", action='store_false')
# parser.add_argument('--heterozygous', help="Use heterozygous variants", action='store_true')
# parser.add_argument('--no-snps', help="Don't use SNPs (default to true).", action='store_false')
# parser.add_argument('--dels', help="Use Deletions.", action='store_true')
# parser.add_argument('--ins', help="Use Insertions.", action='store_true')
# parser.add_argument('--individual', help="The Individual to use.", type=int, default=1)
# parser.add_argument('--append-chromosome', help="Should 'chr' be appended to the chromosome column?.", action='store_true')

def main():
    args = parser.parse_args()
    snps = not args.no_snps
    dels = args.dels
    ins = args.ins
    homs = not args.no_homozygous
    hets = args.heterozygous
    individual = args.individual-1
    fasta_file = fasta.FastaIterator(args.fasta)
    if args.vcf:
        vcf_file = gp.VCFIterator(args.vcf)
    id_tag = args.group_on
    if args.cufflinks:
        gff = gp.GFFReader(args.gff, preset='cufflinks')
    else:
        gff = gp.GFFReader(args.gff, tag_map={'ID': id_tag, 'Parent': 'Parent'})
    with args.out as o:
        for feature_name, feature in gff.feature_map.iteritems():
            if args.cufflinks:
                gff_objects = [(gff_object, gff_object.start) for gff_object in feature.parts() if gff_object.feature_type != 'transcript']
            else:
                gff_objects = [(gff_object, gff_object.start) for gff_object in feature.parts()]
            gff_objects.sort(key=operator.itemgetter(1))
            seq = ''.join([fasta_file.get_sequence(gff_object.seqid, gff_object.start, gff_object.end) for gff_object,_ in gff_objects])
            if gff_object.strand == '-':
                seq = fasta._reverse_complement(seq)
            o.write('>%s\n%s\n' % (feature_name, seq))

if __name__ == "__main__":
    sys.exit(main())