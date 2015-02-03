#!/usr/bin/env python

description = """
This produces a bam file corresponding to junctional regions in a given gtf file
"""

import sys
import pysam
from pythomics.genomics.parsers import GFFReader
from pythomics.templates import CustomParser

parser = CustomParser(description = description)
parser.add_bam()
parser.add_bam_out()
parser.add_gff()

def main():
    args = parser.parse_args()
    samfile = pysam.Samfile(args.bam, 'rb')
    junctionreads = pysam.Samfile(args.out_bam, 'wb', template=samfile)
    id_tag = args.group_on
    chosen_feature = args.feature
    if args.cufflinks:
        gff = GFFReader(args.gff, preset='cufflinks')
    else:
        gff = GFFReader(args.gff, tag_map={'ID': id_tag, 'Parent': 'Parent'})
    written = set([])
    for feature_name, feature in gff.get_features():
        try:
            children = feature.children
        except AttributeError:
            continue
        if len(children) > 1:
            starts = dict([(j.start, j) for i,v in children.iteritems() for j in v.parts()])
            if len(starts) > 1:
                parts = [(v.seqid, v.start, v.end) for i,v in starts.iteritems()]
                parts.sort(key=lambda x: x[1])
                for ri, read in enumerate(parts[:-1]):
                    read2 = parts[ri+1]
                    reads = set([])
                    reads2 = set([])
                    read_dict = {}
                    try:
                        for i in samfile.fetch(read[0], int(read[2])-1, read[2]):
                            if not i.overlap(int(read[2])-1, int(read[2])) or i.qname in written:
                                continue
                            reads.add(i.qname)
                            read_dict[i.qname] = i
                            # if not i.mate_is_unmapped:
                            #     mate = samfile.mate(i)
                            #     reads.add(mate.qname)
                            #     read_dict[mate.qname] = mate
                        for i in samfile.fetch(read2[0], read2[1], int(read2[1])+1):
                            if not i.overlap(int(read2[2])-1, int(read2[2])) or i.qname in written:
                                continue
                            reads2.add(i.qname)
                            read_dict[i.qname] = i
                            # if not i.mate_is_unmapped:
                            #     mate = samfile.mate(i)
                            #     reads2.add(mate.qname)
                            #     read_dict[mate.qname] = mate
                        for i in reads&reads2:
                            written.add(i)
                            junctionreads.write(read_dict[i])
                    except ValueError:
                        continue
    pysam.sort(args.out_bam, '%s_sort'%args.out_bam)
    pysam.index('%s_sort.bam'%args.out_bam)
                    
if __name__ == "__main__":
    sys.exit(main())