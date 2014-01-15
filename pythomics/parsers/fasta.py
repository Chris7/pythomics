import os, re
import templates, config

def _reverse_complement(seq):
    return ''.join([config.BASE_PAIR_COMPLEMENTS[i] for i in reversed(seq)])

def _complement(seq):
    return ''.join([config.BASE_PAIR_COMPLEMENTS[i] for i in seq])

def _translate(seq):
    return ''.join([config.CODON_TABLE.get(seq[i:i+3],'') for i in xrange(0,len(seq.upper()),3)])

class FastaIterator(templates.GenericIterator):
    def __init__(self, filename, delimiter='>', **kwrds):
        """
        Optional argument: delimiter -- > default
        """
        super(FastaIterator, self).__init__(filename)
        self.fasta_file = self.filename
        self.fasta_index = False
        self.delimiter = delimiter
        self.delimiter_length = len(delimiter)
        self.parse = kwrds.get('parse', None)
        if self.parse:
            self.parse = re.compile(self.parse)
        #look for an index
        self.fasta_index = kwrds.get('index', None)
        if not self.fasta_index:
            if os.path.exists(os.path.join(self.fasta_file.name, '.fai')):
                self.fasta_index = os.path.join(self.fasta_file.name, '.fai')
            elif os.path.exists(os.path.join(self.fasta_file.name, '.faidx')):
                self.fasta_index = os.path.join(self.fasta_file.name, '.faidx')
            else:
                self.fasta_index = None
            if self.fasta_index:
                self.open_fasta_index(self.fasta_index)
        self.sequence_index = {}
        self.row = None
        
    def __iter__(self):
        return self
    
    def next(self):
        seq = ""
        row = self.row
        # remove new lines and read in our header
        while not row:
            row = self.fasta_file.next()
        if self.parse:
                header = self.parse.match(row)
        else:
            header = row.strip()[self.delimiter_length:]
        #get our sequence
        row = self.fasta_file.next()
        while row and row[0:self.delimiter_length] != self.delimiter:
            seq+=row.strip()
            try:
                row = self.fasta_file.next()
            except StopIteration:
                return header,seq
        self.row = row
        if header and seq:
            return header,seq
        
    def open_fasta_index(self):
        """
        custom type for file made w/ buildFastaIndex, fai for ones made with samtools
        """
        index = self.fasta_index
        try:
            handle = open(index, 'rb')
        except IOError:
            print 'index not found, creating it'
            try:
                self.build_fasta_index()
                return
            except IOError:
                raise IOError("Index File "+self.fasta_index+"can't be found nor created, check file permissions")
        self.sequence_index = {}
        _seq_dict = self.sequence_index
        for row in handle:
            entry = row.strip().split('\t')
            #stored as: {header: length, # of chars to end of this header, length of fasta lines, length of each line including breakchar}
            _seq_dict[entry[0]] = (entry[1], entry[2], entry[3], entry[4])
            
    def get_sequence(self, chrom, start, end, strand='+', indexing=(-1,0)):
        """
        chromosome is entered relative to the file it was built with, so it can be 'chr11' or '11',
        start/end are coordinates, which default to python style [0,1) internally. So positions should be
        entered with (1,1) indexing. This can be changed with the indexing keyword.
        The default is for everything to be relative to the positive strand
        """
        try:
            divisor = int(self.sequence_index[chrom][2])
        except KeyError:
            self.open_fasta_index()
            try:
                divisor = int(self.sequence_index[chrom][2])
            except KeyError:
                return None
        start+=indexing[0]
        end+=indexing[1]
        #go to start of chromosome
        seekpos = int(self.sequence_index[chrom][1])
        #find how many newlines we have
        seekpos+=start+(start)/divisor
        slen = end-start
        endpos=slen+slen/divisor+1 #a hack of sorts but it works and is easy
        self.fasta_file.seek(seekpos, 0)
        output = self.fasta_file.read(endpos)
        output = output.replace('\n', '')
        out = output[:slen]
        if strand == '+' or strand == 1:
            return out
        if strand == '-' or strand == -1:
            return _reverse_complement(out)
    
    def build_fasta_index(self):
        with open('%s.fai'%self.fasta_file.name, 'wb') as o:
            f = self.fasta_file
            if not self.parse:
                chromReg = re.compile(r'%s(.+)'%self.delimiter)
            else:
                chromReg = self.parse
            row = f.readline()
            total_read = 0
            sequence_read = 0
            header_end = {}
            header = None
            with_break = {}
            without_break = {}
            header_order = []
            while row:
                m = chromReg.match(row)
                total_read += len(row)
                if m:
                    header_end[m.group(1)] = total_read
                    if header:
                        #this is always the LAST header
                        self.sequence_index[header] = (sequence_read, header_end[header], without_break[header], with_break[header])
                        sequence_read = 0
                    header = m.group(1)
                    header_order.append(header)
                    sequence_read = 0
                else:
                    if header not in with_break:
                        with_break[header] = len(row)
                    if header not in without_break:
                        without_break[header] = len(row.strip())
                    sequence_read += len(row.strip())
                row = f.readline()
            #for the last one we found
            header_order.append(header)
            self.sequence_index[header] = (sequence_read, header_end[header], without_break[header], with_break[header])
            for header in header_order:
                d = self.sequence_index[header]
                o.write('%s\t%d\t%d\t%d\t%d\n'%(header,d[0],d[1],d[2],d[3]))