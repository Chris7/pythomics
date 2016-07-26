import unittest
import pythomics.proteomics.peptide as peptide
 
class Test_Peptide_Structure(unittest.TestCase):
         
    def test_peptide_pi(self):
        self.assertEqual('%1.4f' % peptide.Peptide('PEPTIDE').getPI(), '2.8812', "Peptide pI test failure.")
        
    def test_peptide_mass(self):
        self.assertEqual(peptide.Peptide('PEPTIDE').getMass(), 799.3599989999999, "Peptide Mass test failure.")
        
    def test_peptide_charge(self):
        self.assertEqual('%1.6f' % peptide.Peptide('PEPTIDE').getCharge(), '-2.998019', "Peptide Mass test failure.")
        
suite = unittest.TestLoader().loadTestsFromTestCase(Test_Peptide_Structure)
unittest.TextTestRunner(verbosity=2).run(suite)