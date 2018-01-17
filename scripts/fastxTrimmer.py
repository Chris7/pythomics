#!/usr/bin/env python

description = """
This script will trim N's from the ends of a fasta/fastq file so it can be
aligned by tophat (which pukes if there are >5 N's. We remove them from the read
ends only)
"""

import sys, re, os, gzip
from itertools import izip
from multiprocessing import Pool
from pythomics.templates import CustomParser
parser = CustomParser(description = description)
parser.add_fasta()
parser.add_read_pair()
parser.add_out()
parser.add_argument('--min-len', help="The minimum read length reads must be after trimming.", type=int, default=25)
parser.add_argument('--prefix', help="If using paired reads, this is the filename prefix.", type=str)
parser.add_argument('--quality', help='If provided, remove qualities below a given score.', type=int, default=0)
parser.add_argument('--chunk', help='How many reads to submit to each core.', type=int, default=1000)
parser.add_argument('--no-gzip', help='To disable compression with gzip.', action='store_false')
# parser.add_argument('--5partial-match', help='This will trim partial matches at the 3\' end of the sequence if there is a match of at least x nucleotides.', type=int, default=0)
# parser.add_argument('--seed-length', help='The seed length for a match.', type=int, default=0)
# parser.add_argument('--mismatches', help='The number of possible mismatches in a sequence.', type=int, default=3)

start_trim = re.compile(r'^N+')
end_trim = re.compile(r'N+$')
global quality_min
global quality_offset
global paired
global read_min
quality_min = 0
quality_offset = 64
paired = False
read_min = 25

#@profile
def trim_read(entries):
    result = []
    for entry in entries:
        read_name, sequence = entry[:2]
        start_cut = start_trim.match(sequence)
        end_cut = end_trim.search(sequence)
        start_cut = start_cut.end() if start_cut else 0
        end_cut = end_cut.start() if end_cut else len(sequence)
        sequence = sequence[start_cut:end_cut]
        out = [read_name, sequence]
        if len(entry) > 2:
            qual_header, qual_seq = entry[2:]
            qual_seq = qual_seq[start_cut:end_cut]
            if quality_min:
                for i, v in enumerate(qual_seq):
                    if ord(v)-quality_offset > quality_min:
                        break
                start_cut = i
                for i, v in enumerate(reversed(qual_seq)):
                    if ord(v)-quality_offset > quality_min:
                        break
                end_cut = len(qual_seq)-i-1 if i else len(qual_seq)
                qual_seq = qual_seq[start_cut:end_cut]
                out[1] = out[1][start_cut:end_cut]
            out += [qual_header, qual_seq]
        result.append(out)
    return result

#@profile
def main():
    global quality_min
    global quality_offset
    global paired
    global read_min
    args = parser.parse_args()
    read_min = args.min_len
    cores = args.p
    if args.left:
        first_file = args.left
        second_file = args.right
        paired = True
    else:
        first_file = args.fasta
    gzip_out = args.no_gzip
    quality_min = args.quality
    if '.fastq' or '.fq' in first_file:
        import pythomics.parsers.fastq as fastq
        iterator = fastq.FastqIterator
        # figure out the quality
        inum = 0
        for read_head, seq, qual_hed, qual in iterator(first_file):
            inum+=1
            if inum > 50:
                break
            if any([i for i in qual if ord(i) - quality_offset < 0]):
                quality_offset -= 31
    else:
        import pythomics.parsers.fasta as fasta
        iterator = fasta.FastaIterator
    if quality_offset < 0:
        sys.stderr.write('Negative quality scores encountered with both phred64 and phred33 encoding. Please inspect your input.\n')
        return 1
    fasta_one = iterator(first_file)
    if paired:
        fasta_two = iterator(second_file)
    chunk_size = args.chunk
    pool = Pool(cores)
    if paired:
        left_base = '%s_1'%args.prefix
        right_base = '%s_2'%args.prefix
        end_names = ['fasta', 'fq', 'fastq', 'fa']
        try:
            file_type = [i for i in end_names if first_file.name.endswith(i)][0]
        except IndexError:
            base_name = os.path.splitext(first_file.name)[0]
            file_type = [i for i in end_names if base_name.endswith(i)][0]
        if gzip_out:
            o = gzip.open('%s.%s.gz'%(left_base, file_type), 'wb')
            o2 = gzip.open('%s.%s.gz'%(right_base, file_type), 'wb')
        else:
            o = open('%s.%s'%(left_base, file_type), 'wb')
            o2 = open('%s.%s'%(right_base, file_type), 'wb')
    else:
        o = args.out
    # do it in chunks of 1000 reads
    read_list = []
    current_list = []
    if paired:
        it = enumerate(izip(fasta_one, fasta_two))
    else:
        it = enumerate(fasta_one)
    break_point = cores*chunk_size
    for index, entry in it:
        if index >= break_point and not index % break_point and current_list:
            read_list.append(current_list)
            l = []
            current_list = l
            if len(read_list) == cores:
                results = pool.map(trim_read, read_list)
                if paired:
                    for result in results:
                        for read_index, reads in enumerate(result):
                            if read_index % 2:
                                # it's odd, write them both
                                # check length
                                if len(first_read[1]) >= read_min and len(reads[1]) >= read_min:
                                    o.write('%s\n'%'\n'.join(first_read))
                                    o2.write('%s\n'%'\n'.join(reads))
                            else:
                                first_read = reads
                else:
                    for result in results:
                        for read_index, reads in enumerate(result):
                            if len(reads[1]) >= read_min:
                                o.write('%s\n'%'\n'.join(reads))
                read_list = []
        if paired:
            current_list.append(entry[0])
            current_list.append(entry[1])
        else:
            current_list.append(entry)
    read_list.append(current_list)
    results = pool.map(trim_read, read_list)
    if paired:
        for result in results:
            for read_index, reads in enumerate(result):
                if read_index % 2:
                    # it's odd, write them both
                    # check length
                    if len(first_read[1]) >= read_min and len(reads[1]) >= read_min:
                        o.write('%s\n'%'\n'.join(first_read))
                        o2.write('%s\n'%'\n'.join(reads))
                else:
                    first_read = reads
    else:
        for result in results:
            for read_index, reads in enumerate(result):
                if len(reads[1]) >= read_min:
                    o.write('%s\n'%'\n'.join(reads))
    o.flush()
    o.close()
    if paired:
        o2.flush()
        o2.close()

if __name__ == "__main__":
    sys.exit(main())