import os, unittest, hashlib, subprocess
 
class Test_Script_Fasta_Digest(unittest.TestCase):

    def setUp(self):
        self.script_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'scripts')
        self.nucleotide_file = os.path.join(os.path.dirname(__file__), 'nucleotide.fasta')
        self.iterator_file = os.path.join(os.path.dirname(__file__), 'test_fasta_iterator.fasta')
         
    def test_fasta_genome_digest(self):
        job = subprocess.Popen([os.path.join(self.script_dir, 'fastadigest.py'), '-f', self.nucleotide_file, '--genome', '--min', '0',
                                '--frame', '6', '--max', '99999', '--enzyme', 'none'], stdout=subprocess.PIPE)
        digest = hashlib.sha224(job.stdout.read()).hexdigest()
        self.assertEqual( '21ded78615bbfb3ef03969fee9a3d0475f8dac2e77ac7026039b693d', digest, "Fasta Genome Digest Test 1 Failure")
        job = subprocess.Popen([os.path.join(self.script_dir, 'fastadigest.py'), '-f', self.nucleotide_file, '--genome', '--min', '0',
                                '--frame', '6', '--max', '99999', '--enzyme', 'trypsin'], stdout=subprocess.PIPE)
        digest = hashlib.sha224(job.stdout.read()).hexdigest()
        self.assertEqual( 'b79ebb042e5139f0a1608793968f1c26428cb79195d91678b997bbfc', digest, "Fasta Genome Digest Test 2 Failure")
        job = subprocess.Popen([os.path.join(self.script_dir, 'fastadigest.py'), '-f', self.nucleotide_file, '--genome', '--min', '6',
                                '--frame', '6', '--max', '30', '--enzyme', 'trypsin'], stdout=subprocess.PIPE)
        digest = hashlib.sha224(job.stdout.read()).hexdigest()
        self.assertEqual( '123c89194061247b7bfa0f73c9e34e81dccf74465cfa87ffd82973e8', digest, "Fasta Genome Digest Test 3 Failure")
        
    def test_fasta_protein_digest(self):
        job = subprocess.Popen([os.path.join(self.script_dir, 'fastadigest.py'), '-f', self.iterator_file, '--min', '0',
                                '--max', '99999', '--enzyme', 'trypsin'], stdout=subprocess.PIPE)
        digest = hashlib.sha224(job.stdout.read()).hexdigest()
        self.assertEqual( '4d905d579e244676605cda804b30a09428d2ea8b4e6867f0ff7b7891', digest, "Fasta Protein Digest Test 1 Failure")
        job = subprocess.Popen([os.path.join(self.script_dir, 'fastadigest.py'), '-f', self.iterator_file, '--min', '6',
                                '--max', '30', '--enzyme', 'trypsin'], stdout=subprocess.PIPE)
        out = job.stdout.read()
        digest = hashlib.sha224(out).hexdigest()
        self.assertEqual( '23bacc02c77c49287ae35d9d02baca0620b8788e9d87432a73a8600b', digest, "Fasta Protein Digest Test 2 Failure")
        
suite = unittest.TestLoader().loadTestsFromTestCase(Test_Script_Fasta_Digest)
unittest.TextTestRunner(verbosity=2).run(suite)