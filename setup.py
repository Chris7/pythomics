from setuptools import setup

setup(name='pythomics',
      version='0.3.42',
      description='A multi-omic python package',
      url='https://github.com/pandeylab/pythomics',
      author='Chris Mitchell',
      author_email='cmitch48@jhmi.edu',
      install_requires = ['six'],
      license='GPL3',
      packages=['pythomics','pythomics.parsers', 'pythomics.proteomics',
                'pythomics.templates', 'pythomics.tests', 'pythomics.genomics', 'pythomics.utils'],
      scripts=['pythomics/scripts/fastadigest.py', 'pythomics/scripts/fastadigeststats.py',
               'pythomics/scripts/incorporateVCF.py', 'pythomics/scripts/fetchOrfs.py',
               'pythomics/scripts/incorporateGFF.py', 'pythomics/scripts/proteinInference.py',
               'pythomics/scripts/featureCollapser.py', 'pythomics/scripts/fastxTrimmer.py',
               'pythomics/scripts/intersectFiles.py', 'pythomics/scripts/junctionalReads.py',
               'pythomics/scripts/ptmSummary.py'],
      zip_safe=False)
