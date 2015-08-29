#!/usr/bin/env python

description = """
This script will digest a given protein fasta file with the specified enzymes
and summarize how much of the proteome is covered, what residues are missed,
and what isoforms can be uniquely identified.
"""

import sys, copy, re
import pythomics.proteomics.digest as digest
import pythomics.parsers.fasta as fasta
import pythomics.proteomics.config as config
from pythomics.templates import CustomParser

parser = CustomParser(description = description)
parser.add_fasta()
parser.add_out()
parser.add_enzyme(help="Enzyme to use. Pass a list like \"trypsin lysc\" to use multiple enzymes.  "
                    "The order of enzymes will be the order of digestion if digesting in series.")
parser.add_argument('--parallel', help="Should cleavages be done in parallel (default is serial digestion)?", action='store_true')


def main():
    args = parser.parse_args()
    digest_min = args.min
    digest_max = args.max
    enzymes = args.enzyme
    peptides_found = {}
    retained = {}
    total = 0
    proteinMap = {}
    coverageMap = {}
    aas = config.RESIDUE_MASSES.keys()
    aas.sort()
    tlen = 0
    parallel = args.parallel
    for protease_index,protease in enumerate(enzymes):
        if parallel or protease_index == 0:
            fasta_file = fasta.FastaIterator(args.fasta)
        enzyme = digest.Enzyme(enzyme=protease)
        sys.stderr.write('processing %s\n' % protease)
        #if doing in series, this iterator is not reset and will never run
        for header, sequence in fasta_file:
            if protease_index == 0:
                total += 1
                proteinMap[header] = sequence
                tlen += len(sequence)
            for peptide in set(enzyme.cleave(sequence, min=digest_min, max=999999)):
                if len(peptide) > digest_max:
                    #we don't see this one
                    if not parallel:
                        try:
                            retained[header].add(peptide)
                        except KeyError:
                            retained[header] = set([peptide])
                else:
                    #we see this one
                    try:
                        peptides_found[peptide].add(header)
                    except KeyError:
                        peptides_found[peptide] = set([header])
                    try:
                        coverageMap[header].add(peptide)
                    except KeyError:
                        coverageMap[header] = set([peptide])
        if not parallel and protease_index > 0:
            for header in retained:
                sequences = copy.deepcopy(retained[header])
                for sequence in sequences:
                    for peptide in set(enzyme.cleave(sequence, min=digest_min, max=999999)):
                        if len(peptide) > digest_max:
                            if not parallel:
                                retained[header].add(peptide)
                        else:
                            try:
                                peptides_found[peptide].add(header)
                            except KeyError:
                                peptides_found[peptide] = set([header])
                            try:
                                coverageMap[header].add(peptide)
                            except KeyError:
                                coverageMap[header] = set([peptide])
        sys.stderr.write('%d total peptides after digesting with %s\n' % (len(peptides_found),protease))
        if parallel:
            args.fasta.seek(0)
    unique_proteins = set([])
    for peptide in peptides_found:
        if len(peptides_found[peptide]) == 1:
            unique_proteins |= peptides_found[peptide]
    with args.out as o:
        o.write('Protein\tDetectable Length\tTotal Length\tCoverage%%\tUnique ID\t%s\n' % '\t'.join(aas))
        sys.stderr.write('%d proteins found out of %d total proteins in database.\n' % (len(coverageMap),total))
        sys.stderr.write('%d of these detectable proteins may be uniquely identified.\n' % (len(unique_proteins)))
        #figure out coverage
        covered = {}
        found_proteins = set([])
        inum=0
        for peptide in peptides_found:
            inum+=1
            if inum % 50000 == 0:
                sys.stderr.write('%d peptides processed\n' % inum)
            for header in peptides_found[peptide]:
                found_proteins.add(header)
                sequence = proteinMap[header]
                found = covered.get(header,set(xrange(len(sequence))))
                sites = [match.start() for match in re.finditer(peptide, sequence)]
                for match_position in sites:
                    found -= set(xrange(match_position,match_position+len(peptide)))
                covered[header] = found
        avg_cov = 0
        missed_len = 0
        detected = 0
        for header in coverageMap:
            total_len = len(proteinMap[header])
            found_len = total_len-len(covered.get(header,[]))
            perc_cov = float(found_len)/float(total_len)
            o.write('%s\t%d\t%d\t%d\t%s\t' % (header,found_len,total_len,perc_cov*100.0, str(header in unique_proteins)))
            #what aa's do we miss
            aas_missed = ''.join(proteinMap[header][i] for i in covered[header])
            missed = [aas_missed.count(j) for j in aas]
            missed_len += sum(missed)
            o.write('%s\n' % '\t'.join([str(i) for i in missed]))
            avg_cov += perc_cov
            if header in found_proteins:
                detected += 1
        sys.stderr.write('average coverage is %0.4f over entire proteome\n' % (float(tlen-missed_len)/float(tlen)))
        sys.stderr.write('average coverage is %0.4f over detected proteins\n' % (avg_cov/detected))
            
                    
if __name__ == "__main__":
    sys.exit(main())