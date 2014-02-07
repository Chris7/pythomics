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
    pass

class GTFIterator(templates.GenericIterator):
    pass

class BedIterator(templates.GenericIterator):
    pass

class BamIterator(templates.GenericIterator):
    pass

class SamIterator(templates.GenericIterator):
    pass

class PileupIterator(templates.GenericIterator):
    pass