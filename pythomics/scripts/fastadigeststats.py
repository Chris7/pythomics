#!/usr/bin/env python

"""
This script will digest a given protein fasta file with the specified enzymes
and summarize how much of the proteome is covered, what residues are missed,
and what isoforms can be uniquely identified.

Chris Mitchell - February 3, 2014
"""

import argparse, sys, copy, re
import pythomics.proteomics.digest as digest
import pythomics.parsers.fasta as fasta
import pythomics.proteomics.config as config

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--file', nargs='?', help="The fasta file to process.", type=argparse.FileType('r'), default=sys.stdin)
parser.add_argument('-o', '--out', nargs='?', help="The file to write summary to.", type=argparse.FileType('w'), default=sys.stdout)
parser.add_argument('--min', help="Minimum cleavage length", type=int, default=7)
parser.add_argument('--max', help="Maximum cleavage length", type=int, default=30)
parser.add_argument('--parallel', help="Should cleavages be done in parallel (default is serial digestion)?", action='store_true')
parser.add_argument('--series', help="Should cleavages be done in series? (default)", action='store_true')
parser.add_argument('--enzyme', help="Enzyme to use. Pass a command separated list (no spaces); "
                    "the order of enzymes will be the order of digestion if digesting in series.", type=str, default='trypsin')


def main():
    args = parser.parse_args()
    digest_min = args.min
    digest_max = args.max
    enzymes = args.enzyme.split(',')
    peptides_found = {}
    retained = {}
    total = 0
    proteinMap = {}
    coverageMap = {}
    aas = config.RESIDUE_MASSES.keys()
    aas.sort()
    tlen = 0
    if args.parallel and args.series:
        print 'Unable to do both parallel and sequential digestion in a single run'
        return 1
    if not args.parallel:
        args.series = True
    for protease_index,protease in enumerate(enzymes):
        if args.parallel or protease_index == 0:
            fasta_file = fasta.FastaIterator(args.file)
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
                    if args.series:
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
        if not args.parallel and protease_index > 0:
            for header in retained:
                sequences = copy.deepcopy(retained[header])
                for sequence in sequences:
                    for peptide in set(enzyme.cleave(sequence, min=digest_min, max=999999)):
                        if len(peptide) > digest_max:
                            if args.series:
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
        if args.parallel:
            args.file.seek(0)
    unique_proteins = set([])
    for peptide in peptides_found:
        if len(peptides_found[peptide]) == 1:
            unique_proteins |= peptides_found[peptide]
    with args.out as o:
        o.write('Protein\tDetectable Length\tTotal Length\tCoverage%%\tUnique ID\t%s\n' % '\t'.join(aas))
        sys.stderr.write('%d proteins found out of %d total proteins in database\n' % (len(unique_proteins),total))
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