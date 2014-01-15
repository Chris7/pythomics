from setuptools import setup

setup(name='pythomics',
      version='0.1.4',
      description='A multi-omic python package',
      url='https://github.com/pandeylab/pythomics',
      author='Chris Mitchell',
      author_email='cmitch48@jhmi.edu',
      license='GPL3',
      packages=['pythomics','pythomics.parsers', 'pythomics.proteomics',
                'pythomics.templates', 'pythomics.tests'],
      scripts=['pythomics/scripts/fastadigest.py'],
      zip_safe=False)
