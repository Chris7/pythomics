import os
import tempfile
import unittest

from pythomics.proteomics.parsers import MZMLIterator

from mixins import FixtureMixin


class TestMZMLWriter(FixtureMixin, unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        super(TestMZMLWriter, cls).setUpClass()
        cls.mzml = os.path.join(cls.fixtures, "alkane_mix.mzml")

    def test_writes_mzml_scans(self):
        parser = MZMLIterator(open(self.mzml, "rb"))
        with tempfile.NamedTemporaryFile(delete=False, mode="wb") as tf:
            self.addCleanup(os.remove, tf.name)
            parser.writeScans(handle=tf, scans=["5", "6"])

        new_mzml = MZMLIterator(open(tf.name, "rb"))
        scan_ids = [i.id for i in new_mzml if i]
        self.assertEqual(scan_ids, ["5", "6"])
