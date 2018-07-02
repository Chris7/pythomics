from setuptools import find_packages, setup


setup(
    name='pythomics',
    version='0.3.45',
    description='A multi-omic python package',
    url='https://github.com/pandeylab/pythomics',
    author='Chris Mitchell',
    author_email='chris.mit7@gmail.com',
    install_requires=[
        'lxml',
        'six',
    ],
    extras_require={
      'all': [
          'matplotlib',
          'pandas',
          'pysam',
          'scipy',
      ],
    },

    license='GPL3',
    packages=find_packages(),
    scripts=[
        'scripts/fastadigest.py',
        'scripts/fastadigeststats.py',
        'scripts/incorporateVCF.py',
        'scripts/fetchOrfs.py',
        'scripts/incorporateGFF.py',
        'scripts/proteinInference.py',
        'scripts/featureCollapser.py',
        'scripts/fastxTrimmer.py',
        'scripts/intersectFiles.py',
        'scripts/junctionalReads.py',
        'scripts/ptmSummary.py'
    ],
    include_package_data=True,
    zip_safe=False
)
