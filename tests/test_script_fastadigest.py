import os
import platform
import subprocess
import sys
import unittest


def fix_eol(string):
    if platform.system() == "Windows":
        return string.replace("\r\n", "\n")
    return string


class Test_Script_Fasta_Digest(unittest.TestCase):
    def setUp(self):
        base_dir = os.path.split(__file__)[0]
        data_dir = os.path.join(base_dir, "fixtures")
        self.script_dir = os.path.join(base_dir, "..", "scripts")
        self.nucleotide_file = os.path.join(data_dir, "nucleotide.fasta")
        self.iterator_file = os.path.join(data_dir, "test_fasta_iterator.fasta")

    def test_fasta_genome_digest_no_enzyme(self):
        job = subprocess.Popen(
            [
                sys.executable,
                os.path.join(self.script_dir, "fastadigest.py"),
                "-f",
                self.nucleotide_file,
                "--genome",
                "--min",
                "0",
                "--frame",
                "6",
                "--max",
                "99999",
                "--enzyme",
                "none",
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        res = fix_eol(job.stdout.read().decode())
        assert ">chr17_gl000203_random F:+1 Start:46 End:57\nIFF" in res
        assert ">chr1 F:+3 Start:10416 End:10421\nP\n" in res
        assert (
            ">chr17_gl000203_random F:-1 Start:23 End:100\nGLPYSAGWRVWYNLIRKKFKYLSKS\n"
            in res
        )

    def test_fasta_genome_digest_no_min_length(self):
        job = subprocess.Popen(
            [
                sys.executable,
                os.path.join(self.script_dir, "fastadigest.py"),
                "-f",
                self.nucleotide_file,
                "--genome",
                "--min",
                "0",
                "--frame",
                "6",
                "--max",
                "99999",
                "--enzyme",
                "trypsin",
            ],
            stdout=subprocess.PIPE,
        )
        res = fix_eol(job.stdout.read().decode())
        assert ">chr17_gl000203_random F:-3 Start:879 End:902\nNVNHIIDK\n" in res
        assert ">chr17_gl000203_random F:+3 Start:3834 End:3839\nIK\n" in res
        # Assert we have non K/R endings if its a stop codon (the nucleotide sequence goes to the end of the stop codon
        # even though the peptide sequence does not include it)
        assert (
            ">chr17_gl000203_random F:+3 Start:1737 End:1760\nNTSLEFM\n>chr17_gl000203_random F:+3 Start:1761 End:1775\nTSPPR\n"
            in res
        )

        # check for semi-tryptic digests
        assert ">chr1 F:+1 Start:10231 End:10263\nPLTLTLNPKP\n" in res

    def test_fasta_genome_digest_trypsin_with_min_length(self):
        job = subprocess.Popen(
            [
                sys.executable,
                os.path.join(self.script_dir, "fastadigest.py"),
                "-f",
                self.nucleotide_file,
                "--genome",
                "--min",
                "6",
                "--frame",
                "6",
                "--max",
                "30",
                "--enzyme",
                "trypsin",
            ],
            stdout=subprocess.PIPE,
        )
        res = fix_eol(job.stdout.read().decode()).split("\n")
        # ensure the length boundaries are respected
        for sequence in res[1::2]:
            assert 6 <= len(sequence) <= 30

    def test_fasta_protein_digest_no_min_length(self):
        job = subprocess.Popen(
            [
                sys.executable,
                os.path.join(self.script_dir, "fastadigest.py"),
                "-f",
                self.iterator_file,
                "--min",
                "0",
                "--max",
                "99999",
                "--enzyme",
                "trypsin",
            ],
            stdout=subprocess.PIPE,
        )
        res = fix_eol(job.stdout.read().decode())
        assert ">c3 Pep:48\nDAR\n" in res

    def test_fasta_protein_digest_trypsin(self):
        job = subprocess.Popen(
            [
                sys.executable,
                os.path.join(self.script_dir, "fastadigest.py"),
                "-f",
                self.iterator_file,
                "--min",
                "6",
                "--max",
                "30",
                "--enzyme",
                "trypsin",
            ],
            stdout=subprocess.PIPE,
        )
        res = fix_eol(job.stdout.read().decode())
        # Assert we find semi-tryptic products
        assert ">c1 Pep:11\nNKPGVYTK\n" in res

    def test_fasta_protein_digest_respects_length(self):
        job = subprocess.Popen(
            [
                sys.executable,
                os.path.join(self.script_dir, "fastadigest.py"),
                "-f",
                self.iterator_file,
                "--min",
                "6",
                "--max",
                "30",
                "--enzyme",
                "trypsin",
            ],
            stdout=subprocess.PIPE,
        )
        res = fix_eol(job.stdout.read().decode()).split("\n")
        for sequence in res[1::2]:
            assert 6 <= len(sequence) <= 30


if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(Test_Script_Fasta_Digest)
    unittest.TextTestRunner(verbosity=2).run(suite)
