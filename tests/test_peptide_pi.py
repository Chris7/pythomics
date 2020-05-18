import unittest
import pythomics.proteomics.peptide as peptide


class Test_Peptide_Structure(unittest.TestCase):
    def test_peptide_pi(self):
        self.assertEqual(
            "%1.4f" % peptide.Peptide("PEPTIDE").getPI(),
            "2.8812",
            "Peptide pI test failure.",
        )

    def test_peptide_mass(self):
        self.assertAlmostEqual(
            peptide.Peptide("PEPTIDE").getMass(), 799.359964, places=6
        )

    def test_peptide_charge(self):
        self.assertEqual(
            "%1.6f" % peptide.Peptide("PEPTIDE").getCharge(), "-2.998019",
        )


if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(Test_Peptide_Structure)
    unittest.TextTestRunner(verbosity=2).run(suite)
