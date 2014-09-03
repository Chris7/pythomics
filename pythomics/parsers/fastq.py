import pythomics.templates as templates

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
        # return is sequence header, sequence, quality header, quality sequence
        return (self.fastq_file.next().strip(), self.fastq_file.next().strip() ,self.fastq_file.next().strip(), self.fastq_file.next().strip())
