#!/usr/bin/env python

__author__ = 'chris'

description = """
This script will annotate a tab delimited text file with peptides with
corresponding proteins present in an annotation file, and can also
use this annotation to include iBAQ measures.
"""

import argparse, sys, re, csv, copy, decimal, itertools, os
from multiprocessing import Pool, Value
from collections import Counter
from pythomics.templates import CustomParser
import pythomics.proteomics.config as config
import pythomics.proteomics.digest as digest
import pythomics.parsers.fasta as fasta
from pythomics.utils import ColumnFunctions

parser = CustomParser(description = description)
parser.add_fasta(help="The fasta file to match peptides against.")
parser.add_out(help="The name of the file you wish to create with results appended.")
parser.add_argument('--peptide-out', nargs='?', help="The file to write digested products to.", type=argparse.FileType('w'), default=os.devnull)
parser.add_argument('--protein-out', nargs='?', help="The file to write grouped products to.", type=argparse.FileType('w'), default=os.devnull)
parser.add_argument('--mod-out', nargs='?', help="The file to write a modification-centric summary to.", type=argparse.FileType('w'), default=os.devnull)
parser.add_argument('--mod-col', help="The column to append for modifications (if you want to report ratios, etc.).", type=str)
parser.add_argument('--mod-col-func', help="The function to apply to grouped entries in modification columns.", type=str, default='concat')
parser.add_argument('--strict', help='For numeric operations, fail if types are incorrect (converting NA to a float for instance).', action='store_true')
parser.add_delimited_file()
parser.add_argument('-r', '--regex', help="A perl regular expression determining which parts of the header to capture.", type=str)
parser.add_argument('--inferred-name', help="The name you want to assign for protein inference (in case you are regexing for gene names or something).", type=str, default='Proteins')
parser.add_argument('--no-inference', help="Do not append proteins inferred from sequences.", action='store_false')
group = parser.add_argument_group('iBAQ related options')
group.add_argument('--ibaq', help="Provide to append iBAQ values as well (requires protein inference).", action='store_true')
group.add_argument('--precursors', help="The column with precursor area (defaults to header lines containing 'Precursor').", default=None)
parser.add_enzyme()
group.add_argument('--no-normalize', help="Don't normalize iBAQ to total intensity", action='store_false')
protein_group = parser.add_argument_group('Protein Grouping Options')
protein_group.add_argument('--unique-only', help="Only group proteins with unique peptides", action='store_true')
protein_group.add_argument('--position', help="Write the position of the peptide matches.", action='store_true')
protein_group.add_argument('--case-sensitive', help="Treat peptides as case-sensitive (ie separate modified peptides)", action='store_true')
mod_group = parser.add_argument_group('Peptide Modification Options')
mod_group.add_argument('--modification-site', help="Write the position in the parent protein of the modification (requires case-sensitive and modifications being lower-cased).", action='store_true')
motif_group = mod_group.add_argument_group('Motif Options')
motif_group.add_argument('--motifs', help="Enable creation of motifs for each modification.", action='store_true')
motif_group.add_argument('--motif-window', help="The width of the motif window (how many residues to go out from each modification).", type=int, default=10)
motif_group.add_argument('--motif-unique', help="Only output motifs where the peptide mapping is unambiguous.", action='store_true')
motif_group.add_argument('--motif-out', help="Where to save the file with motifs. Default: --out file with _motif suffix.", type=str)


global protein_sequences
global fasta_headers
cache = {}
proteins = []
fasta_headers = []
protein_sequences = ''

def mapper(peptides):
    peptidestring = '(%s)'%'|'.join([i.upper() for i in peptides])
    matches = [(i.group(1), i.start()) for i in re.finditer(peptidestring, protein_sequences)]
    matches.sort(key=lambda x: x[0])
    groups = itertools.groupby(matches, key=lambda x: x[0])
    matched = {}
    for peptide, peptide_grouped in groups:
        pg = list(peptide_grouped)
        indices = [(protein_sequences.count('\n', 0, match), match-protein_sequences[:match].rfind('\n')) for peptide, match in pg]
        found_proteins = list(set([fasta_headers[i[0]] for i in indices]))
        matched[peptide] = {
            'proteins': found_proteins,
            'positions': [i[1] for i in indices],
            'indices': [i[0] for i in indices],
            'matches': [match for peptide, match in pg]
        }
    return matched



def main():
    global protein_sequences
    global fasta_headers
    peptides_mapped = Value('i', 0)
    args = parser.parse_args()
    cores = args.p
    fasta_file = fasta.FastaIterator(args.fasta)
    try:
        peptide_column = int(args.col)-1
    except:
        peptide_column = None
    tsv_file = args.tsv
    out_file = args.out
    header_lines = args.header
    delimiter = args.delimiter
    inference = args.no_inference
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
    if mod_site:
        case_sens = True
        inference = True
    precursor_columns = [int(i) for i in args.precursors.split(',')] if args.precursors else None
    if ibaq:
        enzyme = digest.Enzyme( enzyme=args.enzyme )
    fasta_headers, protein_sequences = zip(*[(header.replace(';', ''), sequence) for header, sequence in fasta_file])
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
    mod_grouping = {} # ordered by protein, site, type
    pep_count = 0
    pep_set = set([])
    mod_col = args.mod_col
    correlator = ColumnFunctions(args)
    motif_search = args.motifs
    if motif_search:
        motif_window = args.motif_window
        motif_unique = args.motif_unique
        if args.motif_out:
            motif_out = open(args.motif_out, 'wb')
        elif args.out:
            motif_out = open('{0}_motif'.format(args.out.name), 'wb')
        else:
            sys.stderr.write("You must provide an output name for motif-out if you are piping to stdout\n")
            return -1
    mod_col_func = getattr(correlator, args.mod_col_func, correlator.concat)
    with tsv_file as f:
        reader = csv.reader(f, delimiter=delimiter)
        for line_num, entry in enumerate(reader):
            if line_num < header_lines:# we assume the first header line is the one we care about
                if peptide_column is None:
                    for i,v in enumerate(entry):
                        if v.lower() == args.col.lower():
                            peptide_column = i
                            break
                if mod_col is not None and mod_col.isdigit():
                    mod_col = int(mod_col)-1
                elif mod_col is not None:
                    for i,v in enumerate(entry):
                        if v.lower() == args.mod_col.lower():
                            mod_col = i
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
                            peptide_history[peptide]['intensities'][n_i].add(intensity)
                        except KeyError:
                            peptide_history[peptide] = {'intensities': {}}
                            for i in xrange(len(precursor_columns)):
                                peptide_history[peptide]['intensities'][i] = set([])
                            peptide_history[peptide]['intensities'][n_i].add(intensity)
                if mod_col is not None:
                    try:
                        peptide_history[peptide]['mod_col'] = entry[mod_col]
                    except KeyError:
                        peptide_history[peptide] = {'mod_col': entry[mod_col], 'intensities': {}}
        if ibaq and normalize:
            for peptide in peptide_history:
                for i, v in peptide_history[peptide]['intensities'].iteritems():
                    normalizations[i] += sum(v)
        else:
            normalizations = [1 for peptide in peptide_history for intensities in peptide_history[peptide]['intensities']]

    # map our peptides is a multi-cored manner
    pool = Pool(cores)
    # get our matches
    plen = 100#len(peptide_history)/cores+1
    peptides = peptide_history.keys()
    # break into groups of 100 (empirically gives fastest mapping)
    subpeptides = [peptides[n:n+plen] for n in xrange(0, len(peptides), plen)]
    results = pool.map(mapper, subpeptides)
    mapped_peptides = dict((k, v) for d in results for (k, v) in d.items())

    protein_grouping = {}
    peptide_grouping = {}
    stats = {'peptides': pep_count}
    stats['peptides_found'] = len(pep_set)
    proteins_mapped = set([])
    peptide_out = []
    empty_dict = {'proteins': '', 'positions': [], 'indices': [], 'matches': []}
    for index, (peptide, d) in enumerate(peptide_history.iteritems()):
        try:
            peptide_dict = peptide_grouping[peptide]
        except KeyError:
            peptide_dict = {'intensities': {}}
            peptide_grouping[peptide] = peptide_dict
        if not index%1000:
            sys.stderr.write('%d of %d complete.\n' % (index, len(peptide_history)))
        mapped_info = mapped_peptides.get(peptide.upper(), empty_dict)
        precursor_int = float(sum([sum(d['intensities'][i]) for i in d['intensities']]))
        entry = [peptide, sum([len(d['intensities'][i]) for i in d['intensities']]), precursor_int]
        if 'inference' not in peptide_dict:
            peptide_dict['inference'] = {'proteins': False}
            # if inference or ibaq:
                # sites = [match.start() for match in re.finditer(peptide.upper(), protein_sequences)]
                # if out_position or mod_site:
                #      indices = [(protein_sequences.count('\n', 0, match_position), match_position-protein_sequences[:match_position].rfind('\n')) for match_position in sites]
                # else:
                #      indices = [(protein_sequences.count('\n', 0, match_position), 0) for match_position in sites]
            if inference:
                proteins = mapped_info['proteins']#[fasta_headers[i[0]] for i in indices]
                if mod_site:
                    start_positions = mapped_info['positions']#[i[1] for i in indices]
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
                    motifs_found = {}
                    find_motif = False
                    if motif_search and len(proteins) == 1 or not motif_unique:
                        find_motif = True
                    for start_position, protein in zip(start_positions, proteins):
                        mod_site_addition = []
                        for j, k in enumerate(peptide):
                            if k.islower():
                                mod_pos = start_position+j
                                mod_key = '%s:%d'%(k, mod_pos)
                                if find_motif:
                                    motif_sequence = protein_sequences[mapped_info['matches'][0]+j-motif_window:mapped_info['matches'][0]+j+motif_window+1]
                                    motif_pos = motif_window
                                    # remove any newlines to the left of us
                                    cut = motif_sequence[:motif_pos].rfind('\n')
                                    if cut != -1:
                                        motif_sequence = motif_sequence[cut+1:]
                                        motif_pos -= (cut+1)
                                    cut = motif_sequence[motif_pos+1:].rfind('\n')
                                    if cut != -1:
                                        motif_sequence = motif_sequence[:motif_pos+cut]
                                    motifs_found[mod_key] = motif_sequence
                                mod_site_addition.append(mod_key)
                                try:
                                    mod_grouping[protein][mod_key]['values'].add(d['mod_col'])
                                    mod_grouping[protein][mod_key]['peptides'].add(peptide)
                                except KeyError:
                                    try:
                                        mod_grouping[protein][mod_key] = {'values': set([d['mod_col']]), 'peptides': set([peptide])}
                                    except KeyError:
                                        mod_grouping[protein] = {mod_key: {'values': set([d['mod_col']]), 'peptides': set([peptide])}}
                        mod_site_additions.append('%s(%s)'%(protein,','.join(mod_site_addition)))
                    peptide_dict['inference']['mod_sites'] = ';'.join(mod_site_additions)
                    peptide_dict['inference']['motifs'] = motifs_found
                if out_position:
                    peptide_dict['inference']['matched_positions'] = ','.join(str(i) for i in start_positions)
        if ibaq:
            ibaqs = []
            intensities = [sum(d['intensities'][i]) for i in d['intensities']]
            if normalize:
                precursor_int = sum([intensities[i]/normalizations[i] for i in xrange(len(normalizations))])
                entry.append(precursor_int)
            for protein_index in mapped_info['indices']:
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
        for i in peptide_out:
            writer.writerow(i)
    with motif_out as o:
        writer = csv.writer(o, delimiter=delimiter)
        header = ['Residue', 'Motif']
        if inference:
            header.append(inferred_name)
        writer.writerow(header)
        for peptide, peptide_dict in peptide_grouping.iteritems():
            for motif_key, motif in peptide_dict['inference'].get('motifs', {}).iteritems():
                writer.writerow([motif_key, motif, peptide_dict['inference']['proteins']])
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
                                    mod_prot, mod_prot_sites = mod_prot_.rsplit('(', 1)
                                    if mod_prot == protein:
                                        # print mod_prot_sites
                                        for mod_prot_site in mod_prot_sites[:-1].split(','):
                                            if mod_prot_site:
                                                mod_aa, mod_prot_site = mod_prot_site[:-1].split(':')
                                                mods.add((mod_aa, mod_prot_site))
                    d = protein_grouping[protein][peptide]
                    peptide_psm_count.append((peptide,sum([len(d['intensities'][i]) for i in d['intensities']])))
                    intensities += [sum(d['intensities'][i]) for i in d['intensities']]
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
                                    mod_prot, mod_prot_sites = mod_prot_.rsplit('(', 1)
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
    if args.mod_out:
        with args.mod_out as o:
            writer = csv.writer(o, delimiter=delimiter)
            header = ['Site', inferred_name, 'Peptide']
            if mod_col:
                header.append(args.mod_col)
            writer.writerow(header)
            #mod_grouping[protein] = {'%s%d'%(k, mod_pos): {'values': set([d['mod_col']]), 'peptides': set([peptide])}}
            for protein, sites_dict in mod_grouping.iteritems():
                for site, site_dict in sites_dict.iteritems():
                    entry = [site, protein, ';'.join(site_dict.get('peptides'))]
                    if mod_col:
                        entry.append(mod_col_func(site_dict.get('values',[])))
                    writer.writerow(entry)
    # write stats
    sys.stderr.write('Peptides Searched: %s\n'%stats['peptides'])
    sys.stderr.write('Unique Peptides Found: %s\n'%stats['peptides_found'])
    sys.stderr.write('%s Mapped to: %s\n'%(inferred_name, stats['proteins_mapped']))
    sys.stderr.write('Modifications:\n')
    for site, sites in stats['modifications'].iteritems():
        sys.stderr.write('  %s: %s found with %d potential sites (%d mappings)\n'%(site, total_mods[site], len(sites), len(set([i[0] for i in sites]))))


if __name__ == "__main__":
    sys.exit(main())
