import os, unittest, hashlib, sys
import pythomics.parsers.fasta as parser
 
class Test_Fasta_Iterator(unittest.TestCase):
     
    def setUp(self):
        self.handle = '%s.%s'%(os.path.splitext(__file__)[0], 'fasta')
        self.index = '%s.%s'%(os.path.splitext(__file__)[0], 'fasta.fai')
         
    def test_fasta_iterator(self):
        out = ""
        f = parser.FastaIterator(self.handle, delimiter='>')
        assert(isinstance(f, parser.FastaIterator))
        for header, sequence in f:
            out += ">%s\n%s\n"%(header,sequence)
        digest = hashlib.sha224(out).hexdigest()
        self.assertEqual( 'a4b6987095e97824cbcb36674f9757c4ccfad161eeb9fd8a993e851a', digest, "Fasta Iterator Failure")
        
    def test_fasta_get_sequence(self):
        f = parser.FastaIterator(self.handle, index=self.index)
        out = f.get_sequence('c1', 5, 30)
        digest = hashlib.sha224(out).hexdigest()
        self.assertEqual( 'ddb5a96ada0f651bffeb8ef856c76faf610ca669a68be904b0acb8b8', digest, "Fasta get_sequence #1 Failure")
        f.fasta_file.close()

    def test_fasta_index_build(self):
        f = parser.FastaIterator(self.handle)
        f.build_fasta_index()
        out = '\n'.join([row.strip() for row in open(self.index, 'rb')])
        digest = hashlib.sha224(out).hexdigest()
        self.assertEqual( 'e071a4ec04e59d55231dc667e06b81b17d96fad0d40fe2ac883e9fe3', digest, "Fasta Index Build Failure")
        
        
suite = unittest.TestLoader().loadTestsFromTestCase(Test_Fasta_Iterator)
unittest.TextTestRunner(verbosity=2).run(suite)