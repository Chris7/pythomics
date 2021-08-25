import os
import re
import subprocess
import sys

import six

from .. import templates
from . import config


def _reverse_complement(seq):
    return "".join([config.BASE_PAIR_COMPLEMENTS.get(i, "N") for i in reversed(seq)])


def _complement(seq):
    return "".join([config.BASE_PAIR_COMPLEMENTS.get(i, "N") for i in seq])


def _translate(seq):
    seq = seq.upper()
    return "".join(
        [
            config.CODON_TABLE.get(seq[i : i + 3], "")
            for i in six.moves.range(0, len(seq.upper()), 3)
        ]
    )


class FastaIterator(templates.GenericIterator):
    def __init__(self, filename, delimiter=">", **kwrds):
        """
        Optional argument: delimiter -- > default
        """
        super(FastaIterator, self).__init__(filename)
        self.fasta_index = False
        self.delimiter = delimiter
        self.delimiter_length = len(delimiter)
        self.parse = kwrds.get("parse", None)
        if self.parse:
            self.parse = re.compile(self.parse)
        # look for an index
        self.fasta_index = kwrds.get("index", None)
        if not self.fasta_index:
            if os.path.exists("%s.%s" % (self.filename, "fai")):
                self.fasta_index = "%s.%s" % (self.filename, "fai")
            elif os.path.exists("%s.%s" % (self.filename, "faidx")):
                self.fasta_index = "%s.%s" % (self.filename, "faidx")
            if self.fasta_index:
                self.open_fasta_index()
            else:
                self.fasta_index = None
        self.sequence_index = {}
        self.row = None

    def _next(self):
        seq = ""
        row = self.row
        # remove new lines and read in our header
        while not row:
            row = six.next(self.handle).decode()

        if self.parse:
            header = self.parse.match(row)
        else:
            header = row.strip()[self.delimiter_length :]
        # get our sequence
        row = six.next(self.handle).decode()
        while row and row[0 : self.delimiter_length] != self.delimiter:
            seq += row.strip()
            try:
                row = six.next(self.handle).decode()
            except StopIteration:
                return header, seq
        self.row = row
        return header, seq

    def open_fasta_index(self):
        index = self.fasta_index
        if index is None or not os.path.exists(index):
            sys.stderr.write("index not found, creating it\n")
            try:
                self.build_fasta_index()
            except IOError:
                raise IOError(
                    "Index File "
                    + self.fasta_index
                    + "can't be found nor created, check file permissions"
                )

        self.sequence_index = {}
        _seq_dict = self.sequence_index
        for row in open(self.fasta_index, "rb"):
            entry = row.decode("utf-8").strip().split("\t")
            # Fasta index spec: http://www.htslib.org/doc/faidx.html#DESCRIPTION
            _seq_dict[entry[0]] = tuple(entry[1:])

    def get_sequence(self, chrom, start, end, strand="+", indexing=(-1, 0)):
        """
        chromosome is entered relative to the file it was built with, so it can be 'chr11' or '11',
        start/end are coordinates, which default to python style [0,1) internally. So positions should be
        entered with (1,1) indexing. This can be changed with the indexing keyword.
        The default is for everything to be relative to the positive strand
        """
        try:
            linewidth = int(self.sequence_index[chrom][2])
        except KeyError:
            self.open_fasta_index()
            try:
                linewidth = int(self.sequence_index[chrom][2])
            except KeyError:
                sys.stderr.write(
                    "%s cannot be found within the fasta index file.\n" % chrom
                )
                return ""
        bytewidth = int(self.sequence_index[chrom][3])

        start += indexing[0]
        end += indexing[1]
        # is it a valid position?
        if start < 0 or end > int(self.sequence_index[chrom][0]):
            raise ValueError(
                "The range %d-%d is invalid. Valid range for this feature is 1-%d."
                % (
                    start - indexing[0],
                    end - indexing[1],
                    int(self.sequence_index[chrom][0]),
                )
            )

        # The offset of the start of the chromosome
        seekpos = int(self.sequence_index[chrom][1])

        # Go to the starting line we are reading from
        seekpos += (start // linewidth) * bytewidth

        # Go to the starting base on this line
        seekpos += start % linewidth
        self.handle.seek(seekpos, 0)

        bases_to_read = end - start
        bytes_to_read = (bases_to_read // linewidth) * bytewidth
        bytes_to_read += bases_to_read % linewidth

        output = self.handle.read(bytes_to_read).decode()
        output = output.replace("\n", "")
        if strand == "+" or strand == 1:
            return output
        if strand == "-" or strand == -1:
            return _reverse_complement(output)

    def _build_index(self):
        out = "{}.fai".format(self.filename)
        header_regex = self.parse or re.compile(r"%s(.+)" % self.delimiter)

        with open(out, "w") as index_file:
            bytes_read = 0

            header_bases = 0
            header_name = None
            first_base = None
            line_bases = None
            byte_line_width = None
            for row in self.handle:
                byte_len = len(row)
                bytes_read += byte_len
                line = row.decode()
                cleaned_line = line.strip()
                if cleaned_line:
                    match = header_regex.match(line)
                    if match:
                        # We're always going to be updating the last seen entry
                        if header_name:
                            index_file.write(
                                "{}\t{}\t{}\t{}\t{}\n".format(
                                    header_name,
                                    header_bases,
                                    first_base,
                                    line_bases,
                                    byte_line_width,
                                )
                            )
                        header_name = match.group(1)
                        header_bases = 0
                        first_base = None
                        line_bases = None
                        byte_line_width = None
                    else:
                        header_bases += len(cleaned_line)
                        if first_base is None:
                            byte_line_width = byte_len
                            first_base = bytes_read - byte_len
                            line_bases = len(cleaned_line)
                else:
                    continue
            # Write the last entry now that we're done
            index_file.write(
                "{}\t{}\t{}\t{}\t{}\n".format(
                    header_name, header_bases, first_base, line_bases, byte_line_width
                )
            )

        return out

    def build_fasta_index(self):
        try:
            subprocess.check_output(["samtools", "faidx", self.filename])
        except OSError:
            self.fasta_index = self._build_index()
        else:
            self.fasta_index = "{}.fai".format(self.filename)
