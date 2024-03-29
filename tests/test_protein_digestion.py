import os, unittest, hashlib
import pythomics.parsers.fasta as parser
import pythomics.proteomics.digest as digest


class Test_Protein_Digestion(unittest.TestCase):
    def setUp(self):
        base_dir = os.path.split(__file__)[0]
        data_dir = os.path.join(base_dir, "fixtures")
        self.handle = os.path.join(data_dir, "test_fasta_iterator.fasta")

    def test_protein_digestion(self):
        out = ""
        f = parser.FastaIterator(self.handle, delimiter=">")
        assert isinstance(f, parser.FastaIterator)
        enzyme = digest.Enzyme(enzyme="trypsin")
        assert isinstance(enzyme, digest.Enzyme)
        for __, sequence in f:
            out += sequence
        peptides = "".join(enzyme.cleave(out, min=7, max=30))
        hash_sum = hashlib.sha224(peptides.encode("utf-8")).hexdigest()
        self.assertEqual(
            "31c6612b85dcea10c26e35826f4e5577b674624725477eb5202b18bb",
            hash_sum,
            "Protein Digestion With Trypsin Failure",
        )
        enzyme = digest.Enzyme(enzyme="lysc")
        peptides = "".join(enzyme.cleave(out, min=0, max=9999, unique=True))
        hash_sum = hashlib.sha224(peptides.encode("utf-8")).hexdigest()
        self.assertEqual(
            "2b5e17ce606e9a296095d8b4b9cf75d44ba662d5eb3531e0a187def4",
            hash_sum,
            "Unique Protein Digestion with Lys-C Failure",
        )


if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(Test_Protein_Digestion)
    unittest.TextTestRunner(verbosity=2).run(suite)
