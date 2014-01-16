import pythomics.templates as templates
import pythomics.genomics.structures as structure
import csv

class VCFIterator(templates.GenericIterator):
    def __init__(self, filename):
        super(VCFIterator, self).__init__(filename)
        #process our meta information
        row = self.filename.readline().strip()
        self.vcf_file = structure.VCFFile(filename)
        while row[:2] == '##':
            if row.startswith('##INFO='):
                assert(self.vcf_file.add_info(row)==True)
            elif row.startswith('##FILTER='):
                assert(self.vcf_file.add_filter(row))
            elif row.startswith('##FORMAT='):
                assert(self.vcf_file.add_format(row))
            elif row.startswith('##CONTIG='):
                assert(self.vcf_file.add_contig(row))
            row = self.filename.readline().strip()
        #got to the end of meta information, we now have the header
        assert(self.vcf_file.add_header(row))
            
    def __iter__(self):
        return self
    
    def next(self):
        row = self.filename.readline()
        while not row:
            row = self.filename.readline()
        return self.vcf_file.add_entry( row.strip() )

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

