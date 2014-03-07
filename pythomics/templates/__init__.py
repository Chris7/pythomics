import gzip
import argparse
import sys
import pythomics.proteomics.config as protein_config

class GenericIterator(object):
    def __init__(self, filename, **kwrds):
        if isinstance(filename, (str, unicode)) and filename.endswith('.gz'):
            self.filename = gzip.GzipFile(filename)
        elif isinstance(filename, (str, unicode)):
            self.filename = open(filename)
        elif isinstance(filename, (file)):
            self.filename = filename
        else:
            raise TypeError
        
    def __iter__(self):
        return self
    
    def next(self):
        raise StopIteration


class CustomParser(argparse.ArgumentParser):
    def __init__(self, *args, **kwrds):
        super(CustomParser, self).__init__(*args, **kwrds)

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
        self.add_argument('-f', '--fasta', nargs='?', help=help, type=argparse.FileType('r'), default=sys.stdin)

    def add_out(self, help='The file to write results to. Leave blank for stdout,'):
        self.add_argument('-o', '--out', nargs='?', help=help, type=argparse.FileType('w'), default=sys.stdout)

    def add_vcf(self, help="The VCF file to use."):
        self.add_argument('--vcf', help=help, type=argparse.FileType('r'), default=sys.stdin)
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
    