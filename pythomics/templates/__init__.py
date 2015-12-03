import gzip
import argparse
import sys
from distutils import spawn
from collections import deque

from six import string_types

import pythomics.proteomics.config as protein_config


class GenericIterator(object):
    gzip = False
    CHUNK_SIZE = 2**16
    UNCONSUMED = ''
    contents = []

    def __init__(self, filename, delimiter='\n', *args, **kwargs):
        self.delimiter = delimiter
        if isinstance(filename, string_types) and filename.endswith('.gz'):
            self.gzip = True
            self.filename = gzip.GzipFile(filename)
        elif isinstance(filename, string_types):
            self.filename = open(filename)
        elif isinstance(filename, (file,)):
            if filename.name.endswith('.gz'):
                self.gzip = True
                self.filename = gzip.GzipFile(filename.name)
            else:
                self.filename = filename
        else:
            raise TypeError

    def __iter__(self):
        return self

    def next(self):
        if self.gzip:
            if self.contents:
                return self.contents.popleft()
            new_contents = self.filename.read(self.CHUNK_SIZE)
            if not new_contents:
                if self.UNCONSUMED:
                    return self.UNCONSUMED
                raise StopIteration
            if new_contents and new_contents[-1] != self.delimiter:
                new_uc_index = new_contents.rfind(self.delimiter)+1
                new_unconsumed = new_contents[new_uc_index:]
                new_contents = new_contents[:new_uc_index]
            else:
                new_unconsumed = ''
            new_contents = self.UNCONSUMED+new_contents
            self.contents = new_contents.split(self.delimiter)
            self.contents = filter(None, self.contents)
            self.UNCONSUMED = new_unconsumed
            self.contents = deque(self.contents)
            if self.contents:
                return self.contents.popleft()
        else:
            return self.filename.next().strip()

    def __next__(self):
        return self.next()


class CustomParser(argparse.ArgumentParser):
    def __init__(self, *args, **kwrds):
        super(CustomParser, self).__init__(*args, **kwrds)
        self.add_argument('-p', help="Threads to run", type=int, default=1)

    def add_enzyme(self, help="The enzyme to cleave with."):
        self.add_argument('--enzyme', help=help, nargs='+',
                          choices=protein_config.ENZYMES.keys(), type=str, default='trypsin')
        self.add_argument('--enzyme-pattern', help='A regex cleavage pattern such as [KR]|{P} to cleave proteins with.', type=str)
        self.add_argument('--min', help="Minimum cleavage length", type=int, default=7)
        self.add_argument('--max', help="Maximum cleavage length", type=int, default=30)

    def add_delimited_file(self, files=None, delimiter=None, cols=None, col_default=1, header=None, help="The delimited file."):
        if files is None:
            files = ['-t', '--tsv']
        if delimiter is None:
            delimiter = ['-d', '--delimiter']
        if cols is None:
            cols = ['-c', '--col']
        if header is None:
            header = ['--header']
        self.add_argument(*files, help=help, type=argparse.FileType('r'), required=True)
        self.add_argument(*delimiter, help="The delimiter for fields.", type=str, default='\t')
        self.add_argument(*cols, help="The column of interest (default: 1). Can be a column name.", type=str, default=col_default)
        self.add_argument(*header, help="The number of headers lines (default: 1).", type=int, default=1)

    def add_fasta(self, help="The fasta file to operate on."):
        self.add_argument('-f', '--fasta', nargs='?', help=help, type=argparse.FileType('r'), required=True)

    def add_read_pair(self):
        self.add_argument('--left', help="The left (5') read pairs", nargs='?', type=argparse.FileType('r'))
        self.add_argument('--right', help="The right (3') read pairs", nargs='?', type=argparse.FileType('r'))

    def add_out(self, help='The file to write results to. Leave blank for stdout,'):
        self.add_argument('-o', '--out', nargs='?', help=help, type=argparse.FileType('w'), default=sys.stdout)

    def add_vcf(self, help="The VCF file to use."):
        self.add_argument('--vcf', help=help, type=argparse.FileType('r'))
        self.add_argument('--no-homozygous', help="Don't include homozygous variants (default to include)", action='store_false')
        self.add_argument('--heterozygous', help="Use heterozygous variants", action='store_false')
        self.add_argument('--no-snps', help="Don't use SNPs (default to true).", action='store_false')
        self.add_argument('--dels', help="Use Deletions.", action='store_true')
        self.add_argument('--ins', help="Use Insertions.", action='store_true')
        self.add_argument('--individual', help="The Individual to use.", type=str, default=None)
        self.add_argument('--append-chromosome', help="Should 'chr' be appended to the chromosome column?.", action='store_true')

    def add_bam(self, help="The SAM/BAM file to use"):
        self.add_argument('--bam', help=help, type=str)

    def add_bam_out(self):
        self.add_argument('--out-bam', type=str)

    def add_gff(self):
        gff_group = self.add_argument_group('GFF file related options')
        gff_group.add_argument('--gff', help="The GFF file to use.", type=argparse.FileType('r'), required=True)
        gff_group.add_argument('--group-on', help="The key to group entries together by (such as transcript_id)", type=str, default='ID')
        gff_group.add_argument('--feature', help="The feature to use for fetching coordinates (such as CDS, does not apply with cufflinks flag).", type=str, default='')
        gff_group.add_argument('--cufflinks', help="If the gff file is in the standard cufflinks output", action='store_true')
        return gff_group

    def add_argument_group(self, *args, **kwargs):
        group = argparse._ArgumentGroup(self, *args, **kwargs)
        self._action_groups.append(group)
        group.add_enzyme = self.add_enzyme
        group.add_fasta = self.add_fasta
        group.add_out = self.add_out
        group.add_vcf = self.add_vcf
        return group

    def add_ms_files(self, help='The corresponding raw mass spectrometry file (ie .raw, .mzml).'):
        self.add_argument('--raw', nargs='?', help=help, type=str)

    def add_processed_ms(self, group=None, help='The corresponding processed mass spectrometry files (ie .msf, .dat).', required=True):
        parser = self if group is None else group
        parser.add_argument('--processed', help=help, type=argparse.FileType('r'), required=required)

    def add_column_function(self, argument,
                            help="The column to apply a function to (if you want to combine fields, sum fields, etc.).",
                            col_argument=None, group=None, col_help="The function to apply.", col_default='concat', parent=True):
        adder = self
        if group:
            adder = group
        from ..utils import ColumnFunctions
        if col_argument is None:
            col_argument = '{}-func'.format(argument)
        if parent:
            adder.add_argument(argument, help=help, type=str)
        adder.add_argument(col_argument, help='{0} Options: {1}'.format(col_help, ', '.join(ColumnFunctions.METHODS)),
                           type=str, choices=ColumnFunctions.METHODS, default=col_default)


