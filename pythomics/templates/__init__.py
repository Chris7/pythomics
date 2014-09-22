import gzip
import argparse
import sys
# from distutils import spawn
from collections import deque
import pythomics.proteomics.config as protein_config


class GenericIterator(object):
    gz = False
    CHUNK_SIZE = 2**16
    UNCONSUMED = ''
    contents = []

    def __init__(self, filename, **kwrds):
        if isinstance(filename, basestring) and filename.endswith('.gz'):
            # if spawn.find_executable('zcat'):
            #     import subprocess
            #     p = subprocess.Popen(['zcat', filename])
            #     from cStringIO import StringIO
            #     self.filename = StringIO(p.communicate()[0])
            # else:
            self.gz = True
            self.filename = gzip.GzipFile(filename)
        elif isinstance(filename, basestring):
            self.filename = open(filename)
        elif isinstance(filename, (file,)):
            if filename.name.endswith('.gz'):
                self.gz = True
                self.filename = gzip.GzipFile(filename.name)
            else:
                self.filename = filename
        else:
            raise TypeError

    def __iter__(self):
        return self
    
    def next(self):
        if self.gz:
            if self.contents:
                return self.contents.popleft()
            new_contents = self.filename.read(self.CHUNK_SIZE)
            if not new_contents:
                if self.UNCONSUMED:
                    uc = self.UNCONSUMED
                    self.UNCONSUMED = ''
                    return self.UNCONSUMED
                raise StopIteration
            if new_contents and new_contents[-1] != '\n':
                new_uc_index = new_contents.rfind('\n')+1
                new_unconsumed = new_contents[new_uc_index:]
                new_contents = new_contents[:new_uc_index]
            else:
                new_unconsumed = ''
            new_contents = self.UNCONSUMED+new_contents
            self.contents = new_contents.split('\n')
            self.UNCONSUMED = new_unconsumed
            self.contents = deque(self.contents)
        else:
            return self.filename.next()


class CustomParser(argparse.ArgumentParser):
    def __init__(self, *args, **kwrds):
        super(CustomParser, self).__init__(*args, **kwrds)
        self.add_argument('-p', help="Threads to run", type=int, default=1)

    def add_enzyme(self):
        self.add_argument('--enzyme', help="The enzyme to cleave with. Also valid is a"
                                           " cleavage pattern such as [KR]|{P}.",
                          choices=protein_config.ENZYMES.keys(), type=str, default='trypsin')
        self.add_argument('--min', help="Minimum cleavage length", type=int, default=7)
        self.add_argument('--max', help="Maximum cleavage length", type=int, default=30)

    def add_delimited_file(self):
        self.add_argument('-t', '--tsv', help="The delimited file.", type=argparse.FileType('r'), required=True)
        self.add_argument('-d', '--delimiter', help="The delimiter for fields.", type=str, default='\t')
        self.add_argument('-c', '--col', help="The column of interest (default: 1).", type=int, default=1)
        self.add_argument('--header', help="The number of headers lines (default: 1).", type=int, default=1)


    def add_fasta(self, help="The fasta file to operate on."):
        self.add_argument('-f', '--fasta', nargs='?', help=help, type=argparse.FileType('r'))

    def add_read_pair(self):
        self.add_argument('--left', help="The left (5') read pairs", nargs='?', type=argparse.FileType('r'))
        self.add_argument('--right', help="The right (3') read pairs", nargs='?', type=argparse.FileType('r'))

    def add_out(self, help='The file to write results to. Leave blank for stdout,'):
        self.add_argument('-o', '--out', nargs='?', help=help, type=argparse.FileType('w'), default=sys.stdout)

    def add_vcf(self, help="The VCF file to use."):
        self.add_argument('--vcf', help=help, type=argparse.FileType('r'))
        self.add_argument('--no-homozygous', help="Don't include homozygous variants (default to include)", action='store_false', default=False)
        self.add_argument('--heterozygous', help="Use heterozygous variants", action='store_true', default=False)
        self.add_argument('--no-snps', help="Don't use SNPs (default to true).", action='store_false', default=False)
        self.add_argument('--dels', help="Use Deletions.", action='store_true', default=False)
        self.add_argument('--ins', help="Use Insertions.", action='store_true', default=False)
        self.add_argument('--individual', help="The Individual to use.", type=int, default=1)
        self.add_argument('--append-chromosome', help="Should 'chr' be appended to the chromosome column?.", action='store_true', default=False)

    def add_argument_group(self, *args, **kwargs):
        group = argparse._ArgumentGroup(self, *args, **kwargs)
        self._action_groups.append(group)
        group.add_enzyme = self.add_enzyme
        group.add_fasta = self.add_fasta
        group.add_out = self.add_out
        group.add_vcf = self.add_vcf
        return group
    
