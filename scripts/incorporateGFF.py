#!/usr/bin/env python

__author__ = 'Chris Mitchell'

from pythomics.templates import CustomParser
import sys
import argparse
import operator
import pythomics.parsers.fasta as fasta
import pythomics.genomics.parsers as gp

description = """
This script will incorporate the a given GFF file into a specified
fasta file. It can also incorporate variants given in a VCF file
while generating this fasta file.
"""

parser = CustomParser(description = description)
parser.add_fasta(help="The fasta file to reference.")
parser.add_out(help="The file to write resulting fasta file to.")
gff_group = parser.add_argument_group('GFF file related options')
gff_group.add_argument('--gff', help="The GFF file to use.", type=argparse.FileType('r'), required=True)
gff_group.add_argument('--group-on', help="The key to group entries together by (such as transcript_id)", type=str, default='ID')
gff_group.add_argument('--feature', help="The feature to use for fetching coordinates (such as CDS, does not apply with cufflinks flag).", type=str, default='')
gff_group.add_argument('--cufflinks', help="If the gff file is in the standard cufflinks output", action='store_true')
vcf_group = parser.add_argument_group('VCF file related options')
vcf_group.add_vcf()
vcf_group.add_argument('--variants-only', help="Only output transcripts with variants.", action='store_true')
splice_group = parser.add_argument_group('Splice Junction Options (if a variant falls over a exon-exon junction. Default is to ignore.)')
splice_group.add_argument('--splice-partial', help="Partially splice variants (only include exonic portions of variant)", action='store_true')


def main():
    args = parser.parse_args()
    snps = args.no_snps
    dels = args.dels
    ins = args.ins
    homs = args.no_homozygous
    hets = args.heterozygous
    individual = 0 if args.individual is None else args.individual
    fasta_file = fasta.FastaIterator(args.fasta)
    splice_variants = args.splice_partial
    id_tag = args.group_on
    vars_only = args.variants_only
    chosen_feature = args.feature
    vcf = None
    if args.vcf:
        vcf = gp.VCFReader(args.vcf, append_chromosome=args.append_chromosome, sample=individual)
    if args.cufflinks:
        gff = gp.GFFReader(args.gff, preset='cufflinks')
    else:
        gff = gp.GFFReader(args.gff, tag_map={'ID': id_tag, 'Parent': 'Parent'})
    with args.out as o:
        for feature_name, feature in gff.feature_map.iteritems():
            header = feature_name
            if args.cufflinks:
                gff_objects = [(gff_object, gff_object.start) for gff_object in feature.parts()
                               if gff_object.feature_type != 'transcript']
            else:
                if chosen_feature:
                    gff_objects = [(gff_object, gff_object.start) for gff_object in feature.parts()
                                   if gff_object.feature_type == chosen_feature]
                else:
                    gff_objects = [(gff_object, gff_object.start) for gff_object in feature.parts()]
            if not gff_objects:
                continue
            if vcf:
                #for vcf, we want to sort from the end to the start so we can incorporate variants without having to
                #worry about an offset
                gff_objects.sort(key=operator.itemgetter(1), reverse=True)
                seq = []
                variant_info = []
                for gff_object, _ in gff_objects:
                    tseq = list(fasta_file.get_sequence(gff_object.seqid, gff_object.start, gff_object.end))
                    overlapping_variants = [(int(entry.pos), entry) for entry in
                                            vcf.contains(gff_object.seqid, gff_object.start, gff_object.end)]
                    #sort our variants from end to start as well
                    overlapping_variants.sort(key=operator.itemgetter(0), reverse=True)
                    to_remove = []
                    for position, vcf_entry in overlapping_variants:
                        checked = False
                        valid_variant = False
                        if homs:
                            if vcf_entry.is_homozygous(individual=individual):
                                if ((snps and not vcf_entry.has_snp(individual=individual)) and
                                    (dels and not vcf_entry.has_deletion(individual=individual)) and
                                    (ins and not vcf_entry.has_insertion(individual=individual))):
                                        continue
                                valid_variant = True
                        if hets:
                            if vcf_entry.is_heterozygous(individual=individual):
                                if ((not valid_variant and not checked) and
                                    (snps and not vcf_entry.has_snp(individual=individual)) and
                                    (dels and not vcf_entry.has_deletion(individual=individual)) and
                                    (ins and not vcf_entry.has_insertion(individual=individual))):
                                        continue
                                valid_variant = True
                        if not valid_variant:
                            to_remove.append(vcf_entry)
                            continue
                        position -= gff_object.start
                        ref = vcf_entry.ref
                        lref = len(ref)
                        if splice_variants and position < 0:
                            if args.splice_partial:
                                alt = max(vcf_entry.get_alt(individual=individual))
                                alt = ''.join(list(alt)[abs(position):])
                                lref += position
                                position=0
                        elif position > 0:
                            alt = max(vcf_entry.get_alt(individual=individual))
                        else:
                            continue
                        variant_info.append('%s %s %s->%s' % (vcf_entry.chrom, vcf_entry.pos, ref, alt))
                        tseq[position:position+lref] = list(alt)
                    vcf.remove_variants(to_remove)
                    seq.append(''.join(tseq))
                if variant_info:
                    header += '\t%s' % ';'.join(variant_info)
                seq.reverse()
                seq = ''.join(seq)
            else:
                gff_objects.sort(key=operator.itemgetter(1))
                seq = ''.join([fasta_file.get_sequence(gff_object.seqid, gff_object.start, gff_object.end)
                               for gff_object, _ in gff_objects])
            if gff_object.strand == '-':
                seq = fasta._reverse_complement(seq)
            if seq and (not vars_only or (vars_only and variant_info)):
                o.write('>%s\n%s\n' % (header, seq))

if __name__ == "__main__":
    sys.exit(main())
