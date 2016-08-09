#!/usr/bin/env python

__author__ = 'chris'

description = """
This script will annotate a tab delimited text file with peptides with
corresponding proteins present in an annotation file, and can also
use this annotation to include iBAQ measures.
"""

import argparse, sys, csv, copy, decimal, itertools, os, operator
try:
    import re2 as re
except ImportError:
    import re
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
parser.add_argument('--strict', help='For numeric operations, fail if types are incorrect (converting NA to a float for instance).', action='store_true')
parser.add_delimited_file(cols=['--peptide-col'], col_default='Peptide')
parser.add_argument('-r', '--regex', help="A perl regular expression determining which parts of the header to capture.", type=str)
parser.add_argument('--inferred-name', help="The name you want to assign for protein inference (in case you are regexing for gene names or something).", type=str, default='Proteins')
parser.add_argument('--no-inference', help="Do not append proteins inferred from sequences.", action='store_true')
parser.add_argument('--no-equality', help="Do not consider Leucine and Isoleucine equal for peptide mapping.", action='store_true')
ibaq_group = parser.add_argument_group('iBAQ related options')
ibaq_group.add_argument('--ibaq', help="Provide to append iBAQ values as well (requires protein inference).", action='store_true')
ibaq_group.add_argument('--precursors', help="The column with precursor area (defaults to header lines containing 'Precursor').", type=str)
parser.add_column_function('', col_argument='--ibaq-function', group=ibaq_group, col_help="The function to apply to groups of iBAQ values (for multiple peptide matches).", parent=False)
ibaq_group.add_argument('--non-redundant', help="Use only non-redundant theoretical tryptic peptides for the iBAQ denominator.", action='store_true')
parser.add_enzyme(help="The enzyme used to digest the sample.")
ibaq_group.add_argument('--normalize', help="Normalize iBAQ to total intensity of column (useful for comparing multiple samples).", action='store_true')
protein_group = parser.add_argument_group('Protein Grouping Options')
protein_group.add_argument('--unique-only', help="Only group proteins with unique peptides", action='store_true')
protein_group.add_argument('--position', help="Write the position of the peptide matches.", action='store_true')
protein_group.add_argument('--case-sensitive', help="Treat peptides as case-sensitive (ie separate modified peptides)", action='store_true')
mod_group = parser.add_argument_group('Peptide Modification Options')
mod_group.add_argument('--mod-out', nargs='?', help="The file to write a modification-centric summary to.", type=argparse.FileType('w'), default=None)
mod_group.add_argument('--modification-site', help="Write the position in the parent protein of the modification (requires case-sensitive and modifications being lower-cased).", action='store_true')
parser.add_column_function('--mod-col', help="The column containing modification information.", group=mod_group)
motif_group = mod_group.add_argument_group('Motif Options')
motif_group.add_argument('--motifs', help="Enable creation of motifs for each modification.", action='store_true')
motif_group.add_argument('--motif-window', help="The width of the motif window (how many residues to go out from each modification).", type=int, default=10)
motif_group.add_argument('--motif-unique', help="Only output motifs where the peptide mapping is unambiguous.", action='store_true')
motif_group.add_argument('--motif-out', help="Where to save the file with motifs. Default: --out file with _motif suffix.", type=str)


global protein_sequences
global protein_sequences_converted
global fasta_headers
global il_convert
global maps_done
maps_done=0
cache = {}
proteins = []
fasta_headers = []
protein_sequences = ''
IBAQ_NORMALIZATION = decimal.Decimal(1e6)
# the number of peptides we map per core
peptides_per_core = 100

def progress_update(a, b):
    sys.stderr.write('\r{0:2.2f}% Completed'.format(a/b*100))
    sys.stderr.flush()

def progress_finish():
    sys.stderr.write('\r100% Completed')
    sys.stderr.flush()
    sys.stderr.write('\n')

def make_unique(l):
    seen = set()
    return [x for x in l if x not in seen and not seen.add(x)]
    #return list(OrderedDict.fromkeys(l))

# from profilestats import profile

# @profile(print_stats=10, dump_stats=True)
def mapper(peptides):
    #pep_format = r'(.+?)\t[^\n]+?({})'
    # pep_format = r'^([^\n]+)\t[^\n]+?({})'
    # matches = [{'peptide': i.group(2), 'peptide_start': i.start(2), 'accession': i.group(1)}
    #            for j in peptides
    #            for i in re.finditer(pep_format.format(j.upper().replace('L', '!').replace('!', '[IL]')), protein_sequences, re.M)]
    sequences = protein_sequences.replace('L', 'I')
    matches = []

    for j in peptides:
        pep_match = j.upper().replace('L', 'I')
        pos_start = 0
        pos = sequences.find(pep_match)
        while pos != -1:
            peptide_offset = pos_start+pos
            accession_end = sequences[:peptide_offset].rfind('\t')
            accession_start = sequences[:accession_end].rfind('\n')+1
            accession = protein_sequences[accession_start:accession_end]
            d = {'peptide': j, 'peptide_start': peptide_offset, 'accession': accession}
            matches.append(d)
            pos_start += pos+1
            pos = sequences[pos_start:].find(pep_match)
    matches.sort(key=operator.itemgetter('peptide', 'peptide_start'))
    groups = itertools.groupby(matches, key=lambda x: x['peptide'])
    matched = {}
    for peptide, peptide_grouped in groups:
        pg = list(peptide_grouped)
        pg.sort(key=lambda x: x['peptide_start'])
        accessions = [group['accession'] for group in pg]
        found_proteins = make_unique(accessions)
        matched[peptide] = {
            'proteins': found_proteins,
            'positions': [group['peptide_start']-protein_sequences[:group['peptide_start']].rfind('\t') for group in pg],
            'matches': [group['peptide_start'] for group in pg],
            'accessions': accessions,
            'unique': True,
        }
    # handle the unmatched ones
    for peptide in set(peptides)-set(matched.keys()):
        matched[peptide] = {
            'proteins': [],
            'positions': [],
            'matches': [],
            'accessions': [],
            'unique': True
        }
    return matched

def main():
    global protein_sequences
    global fasta_headers
    global il_convert
    peptides_mapped = Value('i', 0)
    args = parser.parse_args()
    cores = args.p
    fasta_file = fasta.FastaIterator(args.fasta)
    peptide_column = args.peptide_col
    try:
        peptide_index = int(peptide_column)-1
        peptide_column = peptide_index
    except ValueError:
        peptide_index = None
    tsv_file = args.tsv
    il_convert = not args.no_equality
    out_file = args.out
    header_lines = args.header
    delimiter = args.delimiter
    inference = not args.no_inference
    inferred_name = args.inferred_name
    digest_min = args.min
    digest_max = args.max
    normalize = args.normalize
    ibaq = args.ibaq
    ibaq_redunant = not args.non_redundant
    case_sens = args.case_sensitive
    mod_site = args.modification_site
    unique = args.unique_only
    out_position = args.position
    if mod_site:
        case_sens = True
        inference = True
    precursor_columns = [i for i in args.precursors.split(',')] if args.precursors else None
    if ibaq:
        enzyme = digest.Enzyme(enzyme=args.enzyme[0] if isinstance(args.enzyme, list) else args.enzyme)
    sys.stderr.write("Reading in Fasta file.\n")
    fasta_headers, protein_sequences = zip(*[(header.replace(';', ''), sequence) for header, sequence in fasta_file])

    #replace headers with parsed ones
    if args.regex:
        regex = re.compile(args.regex)
        fasta_headers = [regex.search(header) for header in fasta_headers]
        protein_sequences = [protein_sequences[i] for i, v in enumerate(fasta_headers) if v]
        fasta_headers = [' '.join(i.groups()) for i in fasta_headers if i]
        sys.stderr.write('{0} header sequences did not match regex {1} and have been discarded.\n'.format(len(fasta_headers)-len(protein_sequences), args.regex))
    if ibaq:
        ibaq_protein_sequence = {header: sequence for header, sequence in zip(fasta_headers, protein_sequences)}
        cleaved = {}
    protein_sequences = '\n'.join(['{}\t{}'.format(header, sequence) for header, sequence in zip(fasta_headers, protein_sequences)])
    # protein_sequences = '\n'.join(protein_sequences)
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
            sys.stderr.write("You must provide an output name for motif-out if you are piping to stdout.\n")
            return -1
    mod_col_func = getattr(correlator, args.mod_col_func, correlator.concat)
    ibaq_col_func = getattr(correlator, args.ibaq_function, correlator.concat)
    with tsv_file as f:
        reader = csv.reader(f, delimiter=delimiter)
        for line_num, entry in enumerate(reader):
            if line_num < header_lines:# we assume the first header line is the one we care about
                if peptide_index is None:
                    for i,v in enumerate(entry):
                        if v.lower() == args.peptide_col.lower():
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
                try:
                    precursor_columns = [int(i) for i in precursor_columns]
                except ValueError:
                    precursor_columns = [entry.index(i) for i in precursor_columns]
                normalizations = [0 for i in precursor_columns]
            else:
                peptide = entry[peptide_column]
                pep_count+=1
                if not case_sens:
                    peptide = peptide.upper()
                pep_set.add(peptide)
                if peptide not in peptide_history:
                    peptide_history[peptide] = {
                        'intensities': dict([(i, set([])) for i in xrange(len(precursor_columns))]) if precursor_columns is not None else {},
                    }
                if precursor_columns:
                    for n_i, e_i in enumerate(precursor_columns):
                        if entry[e_i]:
                            try:
                                intensity = decimal.Decimal(entry[e_i])
                            except decimal.InvalidOperation:
                                intensity = decimal.Decimal(0)
                            peptide_history[peptide]['intensities'][n_i].add(intensity)
                if mod_col is not None:
                    peptide_history[peptide]['mod_col'] = entry[mod_col]
        if ibaq and normalize and precursor_columns:
            for peptide in peptide_history:
                for i, v in peptide_history[peptide]['intensities'].iteritems():
                    normalizations[i] += sum(v)
        else:
            normalizations = [decimal.Decimal(1) for i in normalizations]

    # map our peptides is a multi-cored manner
    pool = Pool(cores)
    # get our matches

    peptides = list(set([i.upper() for i in peptide_history.keys()]))
    # break into groups of 100 (empirically gives fastest mapping)
    subpeptides = [peptides[n:n+peptides_per_core] for n in xrange(0, len(peptides), peptides_per_core)]
    if n < len(peptides):
        subpeptides.extend(peptides[n+peptides_per_core:])
    num_peps = len(peptides)

    progress_finish()

    sys.stderr.write('Mapping Peptides.\n')
    results = pool.map_async(mapper, subpeptides)
    results.wait()
    mapped_peptides = dict((k, v) for d in results.get() for (k, v) in d.items())
    sys.stderr.write('\nPeptides mapped.\n')

    protein_grouping = {}
    peptide_grouping = {}
    stats = {'peptides': pep_count}
    stats['peptides_found'] = len(pep_set)
    proteins_mapped = set([])
    peptide_out = []
    empty_dict = {'proteins': '', 'positions': [], 'accessions': [], 'matches': [], 'unique': True}
    for index, (peptide, d) in enumerate(peptide_history.iteritems()):
        try:
            peptide_dict = peptide_grouping[peptide]
        except KeyError:
            peptide_dict = {'intensities': {}}
            peptide_grouping[peptide] = peptide_dict
        if not index%100:
            progress_update(index, len(peptide_history))
        mapped_info = mapped_peptides.get(peptide.upper(), empty_dict)
        precursor_int = float(sum([sum(d['intensities'][i]) for i in d['intensities']]))
        entry = [peptide, sum([len(d['intensities'][i]) for i in d['intensities']]), precursor_int]
        if 'inference' not in peptide_dict:
            peptide_dict['inference'] = {'proteins': ''}
            if inference:
                proteins = mapped_info['proteins']
                accessions = mapped_info['accessions']
                start_positions = mapped_info['positions'] if mod_site else []
                proteins_mapped|=set(proteins)
                if unique:
                    proteins = make_unique(proteins)
                    if len(proteins) > 1:
                        mapped_info['unique'] = False
                matches = ';'.join(proteins)
                peptide_dict['inference']['proteins'] = matches
                if not unique or mapped_info['unique']:
                    entry.append(matches)
                else:
                    entry.append('')
                for protein_index, protein in enumerate(proteins):
                    try:
                        protein_grouping[protein][peptide] = d
                    except KeyError:
                        protein_grouping[protein] = {peptide: d}

                if mod_site:
                    mod_site_additions = []
                    motifs_found = {}
                    find_motif = False
                    if motif_search and (len(proteins) == 1 or not motif_unique):
                        find_motif = True
                    for start_position, protein in zip(start_positions, accessions):
                        mod_site_addition = []
                        for j, k in enumerate(peptide):
                            if k.islower():
                                mod_pos = start_position+j
                                mod_key = '%s:%d'%(k, mod_pos)
                                if find_motif:
                                    motif_sequences = [protein_sequences[i+j-motif_window:i+j+motif_window+1] for i in mapped_info['matches']]
                                    motif_pos = motif_window
                                    # remove any newlines to the left of us
                                    for motif_sequence in motif_sequences:
                                        cut = motif_sequence[:motif_pos].rfind('\t')
                                        if cut != -1:
                                            motif_sequence = motif_sequence[cut+1:]
                                            motif_pos -= (cut+1)
                                        cut = motif_sequence[motif_pos+1:].rfind('\t')
                                        if cut != -1:
                                            motif_sequence = motif_sequence[:motif_pos+cut]
                                        found = motifs_found.get(mod_key, [])
                                        motifs_found[mod_key] = make_unique(found+[motif_sequence])
                                mod_site_addition.append(mod_key)
                                if mod_col or mod_site:
                                    try:
                                        mod_values = mod_grouping[protein][mod_key]['values']
                                        mod_peptides = mod_grouping[protein][mod_key]['peptides']
                                        if mod_col:
                                            mod_values.append(d['mod_col'])
                                        mod_grouping[protein][mod_key]['values'] = make_unique(mod_values)
                                        mod_grouping[protein][mod_key]['peptides'] = make_unique(mod_peptides+[peptide])
                                    except KeyError:
                                        try:
                                            mod_grouping[protein][mod_key] = {'values': make_unique([d['mod_col']]) if mod_col else '', 'peptides': make_unique([peptide])}
                                        except KeyError:
                                            mod_grouping[protein] = {mod_key: {'values': make_unique([d['mod_col']]) if mod_col else '', 'peptides': make_unique([peptide])}}
                        mod_site_additions.append('%s(%s)'%(protein,','.join(mod_site_addition)))
                    peptide_dict['inference']['mod_sites'] = ';'.join(mod_site_additions)
                    peptide_dict['inference']['motifs'] = motifs_found
                peptide_dict['inference']['matched_positions'] = ','.join(str(i) for i in start_positions)
        if ibaq:
            ibaqs = []
            intensities = [sum(d['intensities'][i]) for i in d['intensities']]
            try:
                precursor_int = sum([intensities[i]/normalizations[i] for i in xrange(len(normalizations))])
            except decimal.InvalidOperation:
                precursor_int = 0
            entry.append(precursor_int)
            for protein_index in mapped_info['accessions']:
                peptides = cleaved.get(protein_index, None)
                if peptides is None:
                    if ibaq_redunant:
                        peptides = sum([len(enzyme.cleave(ibaq_protein_sequence[protein_accession], min=digest_min, max=digest_max)) for protein_accession in mapped_info['accessions']])
                    else:
                        peptides = len(set([peptide for tryptic_peptides in [enzyme.cleave(ibaq_protein_sequence[possible_protein_index], min=digest_min, max=digest_max) for possible_protein_index in mapped_info['indices']]
                                            for peptide in tryptic_peptides]))
                    cleaved[protein_index] = peptides
                if not peptides:
                    ibaqs.append(0)
                    continue
                # this divides the precursor intensity of the given peptide by the number of theoretically
                #  possible cleaved peptides per protein.
                # If the user is grouping things at a higher level, say the gene level this will output the ibaq
                # per each mapped isoform if that gene has isoforms.
                # if peptide.upper() == 'HMSFHAHVR':
                #     import pdb; pdb.set_trace();
                ibaqs.append(precursor_int/peptides if peptides and precursor_int else 0)
            peptide_dict['inference']['iBAQ'] = ibaq_col_func([int(IBAQ_NORMALIZATION*i) for i in ibaqs]) if ibaqs else 0
            entry.append(peptide_dict['inference']['iBAQ'] if not unique or mapped_info['unique'] else '')
        if out_position:
            entry.append(peptide_dict['inference'].get('matched_positions', '') if not unique or mapped_info['unique'] else '')
        if mod_site:
            entry.append(peptide_dict['inference'].get('mod_sites', '') if not unique or mapped_info['unique'] else '')
        if motif_search:
            entry.append(';'.join(['{}({})'.format(motif_site,';'.join(motifs)) for motif_site,motifs in peptide_dict['inference'].get('motifs',{}).iteritems()]))
        peptide_out.append(entry)
    progress_finish()
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
        if motif_search:
            header.append('Motif')
        writer.writerow(header)
        for i in peptide_out:
            writer.writerow(i)
    if motif_search:
        with motif_out as o:
            writer = csv.writer(o, delimiter=delimiter)
            header = ['Residue', 'Motif']
            if inference:
                header.append(inferred_name)
            writer.writerow(header)
            for peptide, peptide_dict in peptide_grouping.iteritems():
                for motif_key, motifs in peptide_dict['inference'].get('motifs', {}).iteritems():
                    writer.writerow([motif_key, ';'.join(motifs), peptide_dict['inference']['proteins']])
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
                                        for mod_prot_site in mod_prot_sites[:-1].split(','):
                                            if mod_prot_site:
                                                mod_aa, mod_prot_site = mod_prot_site[:-1].split(':')
                                                mods.add((mod_aa, mod_prot_site))
                    d = protein_grouping[protein][peptide]
                    if not unique or mapped_peptides.get(peptide, {}).get('unique'):
                        peptide_psm_count.append((peptide,sum([len(d['intensities'][i]) for i in d['intensities']])))
                        intensities += [sum(d['intensities'][i]) for i in d['intensities']]
                        if ibaq and normalize:
                            try:
                                precursor_int += sum([intensities[i]/normalizations[i] for i in xrange(len(normalizations))])
                            except decimal.InvalidOperation:
                                pass
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
                    ibaq_value = [int(IBAQ_NORMALIZATION*precursor_int/peptides) if peptides and precursor_int else 0]
                    entry.append(ibaq_col_func(ibaq_value))
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
                    if ibaq:
                        entry.append('iBAQ')
                else:
                    peptide = entry[peptide_column]
                    if not case_sens:
                        peptide = peptide.upper()
                    d = peptide_grouping.get(peptide,False)
                    total_mods.update([k for k in peptide if k.islower()])
                    if d:
                        if inference:
                            entry.append(d['inference']['proteins'] if not unique or mapped_peptides.get(peptide, {}).get('unique') else '')
                        if out_position:
                            entry.append(d['inference']['matched_positions'] if not unique or mapped_peptides.get(peptide, {}).get('unique') else '')
                        if mod_site:
                            mod_proteins = d['inference']['mod_sites']
                            peptide_mods = {}
                            mod_entry = []
                            if not unique or mapped_peptides.get(peptide, {}).get('unique'):
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
                        if ibaq:
                            entry.append(d['inference'].get('iBAQ', 0) if not unique or mapped_peptides.get(peptide, {}).get('unique') else 0)
                # if peptide.upper() == 'HYNEAVKR':
                #     import pdb; pdb.set_trace();
                out_writer.writerow(entry)
        stats['modifications'] = mod_stats
    mod_out = args.mod_out if args.mod_out else open(os.path.join('{}_mods'.format(out_file.name)), 'wb')
    with mod_out as o:
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
    if stats['modifications']:
        sys.stderr.write('Modifications:\n')
        for site, sites in stats['modifications'].iteritems():
            sys.stderr.write('  %s: %s found with %d potential sites (%d mappings)\n'%(site, total_mods[site], len(sites), len(set([i[0] for i in sites]))))


if __name__ == "__main__":
    sys.exit(main())
