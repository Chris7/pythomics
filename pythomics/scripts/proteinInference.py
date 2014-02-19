#!/usr/bin/env python

__author__ = 'chris'

description = """
This script will annotate a tab delimited text file with peptides with
corresponding proteins present in an annotation file, and can also
use this annotation to include iBAQ measures.
"""

import argparse, sys, re, csv, copy, decimal
import pythomics.proteomics.config as config
import pythomics.proteomics.digest as digest
import pythomics.parsers.fasta as fasta

parser = argparse.ArgumentParser(description = description)
parser.add_argument('-f', '--fasta', nargs='?', help="The fasta file to cleave.", type=argparse.FileType('r'), default=sys.stdin)
parser.add_argument('-o', '--out', nargs='?', help="The file to write digested products to.", type=argparse.FileType('w'), default=sys.stdout)
parser.add_argument('-t', '--tsv', help="The tab separated file.", type=argparse.FileType('r'), required=True)
parser.add_argument('-d', '--delimiter', help="The delimiter for fields.", type=str, default='\t')
parser.add_argument('-c', '--col', help="The column with peptides (default: 1).", type=int, default=1)
parser.add_argument('--header', help="The number of headers lines (default: 1).", type=int, default=1)
parser.add_argument('-r', '--regex', help="A perl regular expression determining which parts of the header to capture.", type=str)
parser.add_argument('--no-inference', help="Do not append proteins inferred from sequences.", type=str)
group = parser.add_argument_group('iBAQ related options')
group.add_argument('--ibaq', help="Provide to append iBAQ values as well.", action='store_true', default=False)
group.add_argument('--precursors', help="The column with precursor area (defaults to header lines containing 'Precursor').", type=int, default=None)
group.add_argument('--enzyme', help="The enzyme to cleave with.", choices=config.ENZYMES.keys(), type=str, default='trypsin')
group.add_argument('--min', help="Minimum cleavage length", type=int, default=7)
group.add_argument('--max', help="Maximum cleavage length", type=int, default=30)
group.add_argument('--no-normalize', help="Don't normalize iBAQ to total intensity", action='store_false', default=True)
group.add_argument('--case-sensitive', help="Treat peptides as case-sensitive (ie group modifications)", action='store_true', default=False)

def main():
    args = parser.parse_args()
    fasta_file = fasta.FastaIterator(args.fasta)
    peptide_column = args.col-1
    tsv_file = args.tsv
    header_lines = args.header
    delimiter = args.delimiter
    inference = not args.no_inference
    digest_min = args.min
    digest_max = args.max
    sys.stderr.write("Reading in Fasta file.\n")
    normalize = args.no_normalize
    ibaq = args.ibaq
    case_sens = args.case_sensitive
    if ibaq:
        precursor_columns = [int(i) for i in args.precursors.split(',')] if args.precursors else None
        enzyme = digest.Enzyme( enzyme=args.enzyme )
    fasta_headers, protein_sequences = zip(*[(header, sequence) for header, sequence in fasta_file])
    #replace headers with parsed ones
    if args.regex:
        regex = re.compile(args.regex)
        fasta_headers = [regex.search(header) for header in fasta_headers]
        protein_sequences = [protein_sequences[i] for i, v in enumerate(fasta_headers) if v]
        sys.stderr.write("%d header sequences did not match regex and have been discarded.\n" % (len(fasta_headers)-len(protein_sequences)))
        fasta_headers = [' '.join(i.groups()) for i in fasta_headers if i]
    if ibaq:
        ibaq_protein_sequence = copy.deepcopy(protein_sequences)
        cleaved = {}
    protein_sequences = '\n'.join(protein_sequences)
    peptide_history = {}
    with args.out as o:
        writer = csv.writer(o, delimiter=delimiter)
        header = ['Peptide', 'PSMS', 'Total Precursor Area']
        if inference:
            header.append('Proteins')
        if ibaq:
            if normalize:
                header.append('Normalized Precursor Intensity')
            header.append('iBAQ')
        writer.writerow(header)
        with tsv_file as f:
            reader = csv.reader(f, delimiter=delimiter)
            for line_num, entry in enumerate(reader):
                if line_num < header_lines:#we assume the first header line is the one we care about
                    if ibaq:
                        if not precursor_columns:
                            precursor_columns = [i for i, v in enumerate(entry) if 'precursor' in v.lower()]
                        if normalize:
                            normalizations = [0 for i in precursor_columns]
                else:
                    peptide = entry[peptide_column]
                    if not case_sens:
                        peptide = peptide.upper()
                    if ibaq:
                        for n_i, e_i in enumerate(precursor_columns):
                            if entry[e_i]:
                                intensity = decimal.Decimal(entry[e_i])
                                try:
                                    peptide_history[peptide][n_i].add(intensity)
                                except KeyError:
                                    peptide_history[peptide] = {}
                                    for i in xrange(len(precursor_columns)):
                                        peptide_history[peptide][i] = set([])
                                    peptide_history[peptide][n_i].add(intensity)
        if normalize:
            for peptide in peptide_history:
                for i, v in peptide_history[peptide].iteritems():
                    normalizations[i] += sum(v)
        else:
            normalizations = [1 for peptide in peptide_history for intensities in peptide_history[peptide]]
        for index, (peptide, d) in enumerate(peptide_history.iteritems()):
            if not index%1000:
                sys.stderr.write('%d of %d complete.\n' % (index, len(peptide_history)))
            precursor_int = float(sum([sum(d[i]) for i in d]))
            entry = [peptide, sum([len(d[i]) for i in d]), precursor_int]
            if inference or ibaq:
                sites = [match.start() for match in re.finditer(peptide, protein_sequences)]
                indices = [protein_sequences.count('\n', 0, match_position) for match_position in sites]
            if inference:
                matches = ';'.join([fasta_headers[i] for i in indices])
                entry.append(matches)
            if ibaq:
                ibaqs = []
                intensities = [sum(d[i]) for i in d]
                if normalize:
                    precursor_int = sum([intensities[i]/normalizations[i] for i in xrange(len(normalizations))])
                    entry.append(precursor_int)
                for protein_index in indices:
                    peptides = cleaved.get(protein_index,None)
                    if peptides is None:
                        peptides = len(enzyme.cleave(ibaq_protein_sequence[protein_index], min=digest_min, max=digest_max))
                        cleaved[protein_index] = peptides
                    if not peptides:
                        ibaqs.append('NA')
                        continue
                    ibaqs.append(precursor_int/peptides)
                entry.append(';'.join(str(i) for i in ibaqs))
            writer.writerow(entry)

if __name__ == "__main__":
    sys.exit(main())