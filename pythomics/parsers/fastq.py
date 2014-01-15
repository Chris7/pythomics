import templates

class FastqIterator(templates.GenericIterator):
    def __init__(self, filename):
        """
        Optional argument: delimiter -- > default
        """
        super(FastqIterator, self).__init__(filename)
        self.fastq_file = self.filename
        
    def __iter__(self):
        return self
    
    def next(self):
        seq_header = self.fastq_file.next().strip()
        seq = self.fastq_file.next().strip()
        qual_header = self.fastq_file.next().strip()
        qual_seq = self.fastq_file.next().strip()
        return (seq_header,seq,qual_header,qual_seq)