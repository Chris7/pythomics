#!/usr/bin/env python

__author__ = 'chris'

description = """
This script will annotate a tab delimited text file with peptides with
corresponding proteins present in an annotation file, and can also
use this annotation to include iBAQ measures.
"""

import argparse, sys, re, csv, copy, decimal
from collections import Counter
from pythomics.templates import CustomParser
import pythomics.proteomics.config as config
import pythomics.proteomics.digest as digest
import pythomics.parsers.fasta as fasta

parser = CustomParser(description = description)
parser.add_fasta(help="The fasta file to match peptides against.")
parser.add_out(help="The name of the file you wish to create with results appended.")
parser.add_argument('--peptide_out', nargs='?', help="The file to write digested products to.", type=argparse.FileType('w'), default=sys.stdout)
parser.add_argument('--protein_out', nargs='?', help="The file to write grouped products to.", type=argparse.FileType('w'), default=sys.stdout)
parser.add_delimited_file()
parser.add_argument('-r', '--regex', help="A perl regular expression determining which parts of the header to capture.", type=str)
parser.add_argument('--inferred-name', help="The name you want to assign for protein inference (in case you are regexing for gene names or something).", type=str, default='Proteins')
parser.add_argument('--no-inference', help="Do not append proteins inferred from sequences.", action='store_false', default=False)
group = parser.add_argument_group('iBAQ related options')
group.add_argument('--ibaq', help="Provide to append iBAQ values as well (requires protein inference).", action='store_true', default=False)
group.add_argument('--precursors', help="The column with precursor area (defaults to header lines containing 'Precursor').", default=None)
parser.add_enzyme()
group.add_argument('--no-normalize', help="Don't normalize iBAQ to total intensity", action='store_false', default=True)
protein_group = parser.add_argument_group('Protein Grouping Options')
protein_group.add_argument('--unique-only', help="Only group proteins with unique peptides", action='store_true', default=False)
protein_group.add_argument('--position', help="Write the position of the peptide matches.", action='store_true', default=False)
protein_group.add_argument('--case-sensitive', help="Treat peptides as case-sensitive (ie separate modified peptides)", action='store_true', default=False)
protein_group.add_argument('--modification-site', help="Write the position in the parent protein of the modification (requires case-sensitive and modifications being lower-cased).", action='store_true', default=False)




def main():
    args = parser.parse_args()
    fasta_file = fasta.FastaIterator(args.fasta)
    peptide_column = args.col-1
    tsv_file = args.tsv
    out_file = args.out
    header_lines = args.header
    delimiter = args.delimiter
    inference = not args.no_inference
    inferred_name = args.inferred_name
    digest_min = args.min
    digest_max = args.max
    sys.stderr.write("Reading in Fasta file.\n")
    normalize = args.no_normalize
    ibaq = args.ibaq
    case_sens = args.case_sensitive
    mod_site = args.modification_site
    unique = args.unique_only
    out_position = args.position
    precursor_columns = [int(i) for i in args.precursors.split(',')] if args.precursors else None
    if ibaq:
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
    pep_count = 0
    pep_set = set([])
    with tsv_file as f:
        reader = csv.reader(f, delimiter=delimiter)
        for line_num, entry in enumerate(reader):
            if line_num < header_lines:#we assume the first header line is the one we care about
                if not precursor_columns:
                        precursor_columns = [i for i, v in enumerate(entry) if 'precursor' in v.lower()]
                if ibaq:
                    if normalize:
                        normalizations = [0 for i in precursor_columns]
            else:
                peptide = entry[peptide_column]
                pep_count+=1
                if not case_sens:
                    peptide = peptide.upper()
                pep_set.add(peptide)
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
        if ibaq and normalize:
            for peptide in peptide_history:
                for i, v in peptide_history[peptide].iteritems():
                    normalizations[i] += sum(v)
        else:
            normalizations = [1 for peptide in peptide_history for intensities in peptide_history[peptide]]
    protein_grouping = {}
    peptide_grouping = {}
    stats = {'peptides': pep_count}
    stats['peptides_found'] = len(pep_set)
    proteins_mapped = set([])
    with args.peptide_out as o:
        writer = csv.writer(o, delimiter=delimiter)
        header = ['Peptide', 'PSMS', 'Total Precursor Area']
        if inference:
            header.append(inferred_name)
        if ibaq:
            if normalize:
                header.append('Normalized Precursor Intensity')
            header.append('iBAQ')
        if out_position:
            header.append('Peptide %s Position'%inferred_name)
        if mod_site:
            header.append('Modification Positions')
        writer.writerow(header)
        for index, (peptide, d) in enumerate(peptide_history.iteritems()):
            try:
                peptide_dict = peptide_grouping[peptide]
            except KeyError:
                peptide_dict = {}
                peptide_grouping[peptide] = peptide_dict
            if not index%1000:
                sys.stderr.write('%d of %d complete.\n' % (index, len(peptide_history)))
            precursor_int = float(sum([sum(d[i]) for i in d]))
            entry = [peptide, sum([len(d[i]) for i in d]), precursor_int]
            if 'inference' not in peptide_dict:
                peptide_dict['inference'] = {'proteins': False}
                if inference or ibaq:
                    sites = [match.start() for match in re.finditer(peptide.upper(), protein_sequences)]
                    if out_position or mod_site:
                         indices = [(protein_sequences.count('\n', 0, match_position), match_position-protein_sequences[:match_position].rfind('\n')) for match_position in sites]
                    else:
                         indices = [(protein_sequences.count('\n', 0, match_position), 0) for match_position in sites]
                if inference:
                    proteins = [fasta_headers[i[0]] for i in indices]
                    if mod_site:
                        start_positions = [i[1] for i in indices]
                    proteins_mapped|=set(proteins)
                    matches = ';'.join(proteins)
                    if unique:
                        proteins = list(set(proteins))
                    if not unique or len(proteins) == 1:
                        for protein_index, protein in enumerate(proteins):
                            try:
                                protein_grouping[protein][peptide] = d
                            except KeyError:
                                protein_grouping[protein] = {peptide: d}
                    entry.append(matches)
                    peptide_dict['inference']['proteins'] = matches
                    if mod_site:
                        mod_site_additions = []
                        for start_position, protein in zip(start_positions, proteins):
                            mod_site_addition = []
                            for j,k in enumerate(peptide):
                                if k.islower():
                                    mod_site_addition.append('%s:%d'%(k,start_position+j))
                            mod_site_additions.append('%s(%s)'%(protein,','.join(mod_site_addition)))
                        peptide_dict['inference']['mod_sites'] = ';'.join(mod_site_additions)
                    if out_position:
                        peptide_dict['inference']['matched_positions'] = ','.join(str(i[1]) for i in indices)
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
            if out_position:
                # entry.append(','.join(str(i[1]) for i in indices))
                entry.append(peptide_dict['inference']['matched_position'])
            if mod_site:
                # entry.append(';'.join(mod_site_additions))
                entry.append(peptide_dict['inference']['mod_sites'])
            writer.writerow(entry)
    stats['proteins_mapped'] = len(proteins_mapped)
    if inference:
        with args.protein_out as o:
            writer = csv.writer(o, delimiter=delimiter)
            header = [inferred_name, 'Peptides', 'Total Precursor Area']
            if mod_site:
                header.append('Modification Positions')
            if ibaq:
                if normalize:
                    header.append('Normalized Precursor Intensity')
                header.append('iBAQ')
            writer.writerow(header)
            for protein in protein_grouping:
                entry = [protein]
                intensities = []
                precursor_int = 0
                peptide_psm_count = []
                mods = set([])
                for peptide in protein_grouping[protein]:
                    if mod_site:
                        peptide_dict = peptide_grouping.get(peptide, False)
                        if peptide_dict:
                            mod_proteins = peptide_dict['inference']['mod_sites']
                            for mod_protein in mod_proteins.split(';'):
                                #mod protein looks like:
                                #WBGene00004829(y:467,k:471);WBGene00019361(m:68);WBGene00019361(m:118);WBGene00019361(m:68);WBGene00020808(m:261);WBGene00020808(m:156)
                                mod_prots = mod_protein.split(';')
                                for mod_prot_ in mod_prots:
                                    mod_prot, mod_prot_sites = mod_prot_.split('(')
                                    if mod_prot == protein:
                                        for mod_prot_site in mod_prot_sites[:-1].split(','):
                                            if mod_prot_site:
                                                mod_aa, mod_prot_site = mod_prot_site[:-1].split(':')
                                                mods.add((mod_aa, mod_prot_site))
                    d = protein_grouping[protein][peptide]
                    peptide_psm_count.append((peptide,sum([len(d[i]) for i in d])))
                    intensities += [sum(d[i]) for i in d]
                    if ibaq and normalize:
                        precursor_int += sum([intensities[i]/normalizations[i] for i in xrange(len(normalizations))])
                entry.append(';'.join(['%s(%s)' % (i,j) for i,j in peptide_psm_count]))
                entry.append(sum(intensities))
                if mod_site:
                    mods = list(mods)
                    mods.sort(key=lambda x:x[1])
                    entry.append(';'.join(['%s%s'%(i,j) for i,j in mods]))
                if ibaq:
                    if normalize:
                        entry.append(precursor_int)
                    peptides = cleaved.get(protein_index,None)
                    if peptides:
                        entry.append(precursor_int/peptides)
                    else:
                        entry.append('NA')
                writer.writerow(entry)
    tsv_file = open(tsv_file.name)
    with tsv_file as f:
        reader = csv.reader(f, delimiter=delimiter)
        mod_stats = {}
        with out_file as o:
            out_writer = csv.writer(o, delimiter=delimiter)
            total_mods = Counter()
            for line_num, entry in enumerate(reader):
                if line_num < header_lines:#we assume the first header line is the one we care about
                    if inference:
                        entry.append(inferred_name)
                    if out_position:
                        entry.append('Peptide %s Position'%inferred_name)
                    if mod_site:
                        entry.append('Modification Position')
                else:
                    peptide = entry[peptide_column]
                    if not case_sens:
                        peptide = peptide.upper()
                    d = peptide_grouping.get(peptide,False)
                    total_mods.update([k for k in peptide if k.islower()])
                    if d:
                        if inference:
                            entry.append(d['inference']['proteins'])
                        if out_position:
                            entry.append(d['inference']['matched_position'])
                        if mod_site:
                            mod_proteins = d['inference']['mod_sites']
                            peptide_mods = {}
                            mod_entry = []
                            for mod_protein in mod_proteins.split(';'):
                                #mod protein looks like:
                                mod_prots = mod_protein.split(';')
                                for mod_prot_ in mod_prots:
                                    if not mod_prot_:
                                        continue
                                    mod_prot, mod_prot_sites = mod_prot_.split('(')
                                    for mod_prot_site in mod_prot_sites[:-1].split(','):
                                        if mod_prot_site:
                                            mod_aa, mod_prot_site = mod_prot_site.split(':')
                                            try:
                                                peptide_mods[mod_prot].add((mod_aa, mod_prot_site))
                                            except KeyError:
                                                peptide_mods[mod_prot] = set([(mod_aa, mod_prot_site)])
                                            try:
                                                mod_stats[mod_aa].add((mod_prot, mod_prot_site))
                                            except KeyError:
                                                mod_stats[mod_aa] = set([(mod_prot, mod_prot_site)])
                            for mod_prot, mods in peptide_mods.iteritems():
                                modl = list(mods)
                                modl.sort(key=lambda x: x[1])
                                mod_entry.append('%s(%s)'%(mod_prot, ' '.join(['%s:%s'%(i,j) for i,j in modl])))
                            entry.append(';'.join(mod_entry))
                out_writer.writerow(entry)
        stats['modifications'] = mod_stats
    # write stats
    sys.stderr.write('Peptides Searched: %s\n'%stats['peptides'])
    sys.stderr.write('Unique Peptides Found: %s\n'%stats['peptides_found'])
    sys.stderr.write('%s Mapped to: %s\n'%(inferred_name, stats['proteins_mapped']))
    sys.stderr.write('Modifications:\n')
    for site, sites in stats['modifications'].iteritems():
        sys.stderr.write('  %s: %s found with %d potential sites (%d mappings)\n'%(site, total_mods[site], len(sites), len(set([i[0] for i in sites]))))


if __name__ == "__main__":
    sys.exit(main())
