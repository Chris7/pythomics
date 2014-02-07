__author__ = 'Chris Mitchell'

import sys
import pythomics.templates as templates
import pythomics.genomics.structures as structure

class VCFIterator(templates.GenericIterator):
    
    def __init__(self, filename):
        """An iterator over the VCF file format
        
        This reads VCF files and has been tested on VCF 4.0.
        The returned items are VCFEntry
        
        """
        super(VCFIterator, self).__init__(filename)
        #process our meta information
        row = self.filename.next().strip()
        self.vcf_file = structure.VCFFile(filename)
        self.inum=0
        while row[:2] == '##':
            if row.startswith('##INFO='):
                assert(self.vcf_file.add_info(row)==True)
            elif row.startswith('##FILTER='):
                assert(self.vcf_file.add_filter(row))
            elif row.startswith('##FORMAT='):
                assert(self.vcf_file.add_format(row))
            elif row.startswith('##CONTIG='):
                assert(self.vcf_file.add_contig(row))
            elif row.startswith('##ALT='):
                assert(self.vcf_file.add_alt(row))
            row = self.filename.next().strip()
        #got to the end of meta information, we now have the header
        assert(self.vcf_file.add_header(row))
            
    def __iter__(self):
        return self
    
    def next(self):
        row = self.filename.next()
        while not row:
            row = self.filename.next()
        self.inum+=1
        if self.inum%100000 == 0:
            sys.stderr.write('Processed %d VCF entries\n' % self.inum)
        return self.vcf_file.parse_entry( row.strip() )

class GFFIterator(templates.GenericIterator):
    def __init__(self, filename, info_delimiter=';', key_delimiter='=', quotechar=''):
        super(GFFIterator, self).__init__(filename)
        self.info_delimiter=info_delimiter
        self.key_delimiter=key_delimiter
        self.quotechar=quotechar

    def __iter__(self):
        return self

    def next(self):
        row = self.filename.next()
        while not row or row[0] == '#':#skip blanks and comments
            row = self.filename.next()
        ob = structure.GFFObject(info_delimiter=self.info_delimiter,
                                 key_delimiter=self.key_delimiter,
                                 quotechar=self.quotechar)
        ob.parse_entry(row)
        return ob

class BedIterator(templates.GenericIterator):
    pass

class BamIterator(templates.GenericIterator):
    pass

class SamIterator(templates.GenericIterator):
    pass

class PileupIterator(templates.GenericIterator):
    pass

class GFFReader(templates.GenericIterator):
    def __init__(self, filename, fast_filter=None, info_delimiter=';', key_delimiter='=', quotechar=''):
        """This will read an entire GFF/GTF/GFF3 file and establish parent-child relations and
         allow for querying based on attributes such as transcripts, gene ids, and allow writing
         of useful file types such as fasta files of coding sequences, etc.

        """
        super(GFFReader, self).__init__(filename)
        self.positions = {}
        self.fast_attributes = {}
        self.filters = set(fast_filter) if fast_filter else []
        for gff_object in GFFIterator(filename, info_delimiter=info_delimiter,
                                      key_delimiter=key_delimiter, quotechar=quotechar):
            chrom, start, end = gff_object.seqid, int(gff_object.start), int(gff_object.end)
            try:
                self.positions[chrom][(start, end)].append(gff_object)
            except KeyError:
                try:
                    self.positions[chrom][(start, end)] = [gff_object]
                except KeyError:
                    self.positions[chrom] = {(start, end): [gff_object]}
            fast_lookup = [fast_access for fast_access in self.filters for attribute in gff_object.attributes if fast_access in attribute]
            for fast_access in fast_lookup:
                try:
                    self.fast_attributes[fast_access].append(gff_object)
                except KeyError:
                    self.fast_attributes[fast_access] = [gff_object]

    def get_attribute(self, attribute, value=None):
        """This returns a list of gff objects with the given attribute and if supplied, those
        attributes with the specified value

        """
        if attribute in self.filters:
            return self.fast_attributes[attribute] if not value else\
                [i for i in self.fast_attributes[attribute] if i.attributes.get(attribute, False)]

    def get_entry(self, seqid, start, end, overlap=True):
        d = self.positions.get(seqid,[])
        if overlap:
            return [gff_object for gff_start, gff_end in d
                    for gff_object in d[(gff_start, gff_end)]
                    if not (end <= gff_start or start >= gff_end)]
        else:
            return [gff_object for gff_start, gff_end in d
                    for gff_object in d[(gff_start, gff_end)]
                    if (gff_start >= start and gff_end >= end)]