__author__ = 'Chris Mitchell'

import os, unittest, hashlib
import pythomics.genomics.parsers as parser

class Test_VCF_Iterator(unittest.TestCase):

    def setUp(self):
        self.handle = os.path.join(os.path.split(__file__)[0], 'valid-4.0.vcf')

    def test_vcf_iterator(self):
        out = ""
        f = parser.VCFIterator( self.handle )
        assert(isinstance( f, parser.VCFIterator))
        out = '\n'.join([str(row) for row in f])
        digest = hashlib.sha224(out).hexdigest()
        self.assertEqual( 'debfe5ab13d9e1fe2abe860aa06ca243c7d61d82fff3c98c3569f711', digest, "VCF Iterator Failure")

    def test_vcf_zygosity(self):
        pass

    def test_vcf_variants(self):
        pass

    def test_vcf_filters(self):
        pass

    def test_vcf_alleles(self):
        pass


suite = unittest.TestLoader().loadTestsFromTestCase(Test_VCF_Iterator)
unittest.TextTestRunner(verbosity=2).run(suite)