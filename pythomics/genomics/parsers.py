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
    def __init__(self, filename, filter=[], info_delimiter=';', key_delimiter='=', quotechar=''):
        """This will read an entire GFF/GTF/GFF3 file and establish parent-child relations and
         allow for querying based on attributes such as transcripts, gene ids, and allow writing
         of useful file types such as fasta files of coding sequences, etc.
         
        """
        super(GFFReader, self).__init__(filename)
        pass