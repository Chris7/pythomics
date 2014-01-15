import os, unittest, hashlib
import parsers.fasta as parser
 
class Test_Fasta_Iterator(unittest.TestCase):
     
    def setUp(self):
        self.handle = '%s.%s'%(os.path.splitext(__file__)[0], 'fasta')
        self.index =  '%s.%s'%(os.path.splitext(__file__)[0], 'fasta.fai')
         
    def test_fasta_iterator(self):
        out = ""
        f = parser.FastaIterator( self.handle, delimiter= '>' )
        assert(isinstance( f, parser.FastaIterator)) 
        for header, sequence in f:
            out += ">%s\n%s\n"%(header,sequence)
        digest = hashlib.sha224(out).hexdigest()
        self.assertEqual( '4208fd8f0e8e8fa4f27e002431d624421da220872d75172db76257a1', digest, "Fasta Iterator Failure")
        
    def test_fasta_index_build(self):
        f = parser.FastaIterator( self.handle )
        f.build_fasta_index()
        out = '\n'.join([row for row in open(self.index)])
        digest = hashlib.sha224(out).hexdigest()
        self.assertEqual( 'e0c787da7b19f79b7c1c261fbd4fd55bd0d7e59adb72bf3c3b9e4f5a', digest, "Fasta Index Build Failure") 
        
    def test_fasta_get_sequence(self):
        f = parser.FastaIterator( self.handle, index = self.index )
        out = f.get_sequence('c1', 5, 30)
        digest = hashlib.sha224(out).hexdigest()
        self.assertEqual( 'ddb5a96ada0f651bffeb8ef856c76faf610ca669a68be904b0acb8b8', digest, "Fasta get_sequence #1 Failure") 
        out = f.get_sequence('c2', 70, 300)
        digest = hashlib.sha224(out).hexdigest()
        self.assertEqual( '755df1058e7dc0f22ecd3ca8d73fc784e78dfdf3d8f984591af75782', digest, "Fasta get_sequence #2 Failure")
        out = f.get_sequence('c2', 59, 80)
        digest = hashlib.sha224(out).hexdigest()
        self.assertEqual( 'b93fb5c1cf028b56b59bd20446bfaaa0a390f6d9a971c3c6cce2e202', digest, "Fasta get_sequence #3 Failure")
        
        
suite = unittest.TestLoader().loadTestsFromTestCase(Test_Fasta_Iterator)
unittest.TextTestRunner(verbosity=2).run(suite)