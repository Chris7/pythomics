__author__ = 'Chris Mitchell'

import sys
import re
from collections import OrderedDict
from pythomics.genomics import config

GENOTYPE_SPLIT = re.compile(r'[\|/]')

class VCFFile(object):
    def __init__(self, filename):
        """VCF File container
        
        This class is designed to be somewhat consistent with the C++
        implementation of vcftools. It contains our higher level
        information such as metadata, samples, the format field,
        and can contain entries as well. Thus it can serve as a
        simple iterator or an in-memory VCF instance for rapid queries.
        
        """
        self.filename = filename
        self.meta = VCFMeta()
        self.n_individuals = 0
        self.individuals = OrderedDict() # maps individual name to index
        self.entries = {}

    def add_extra(self, entry):
        """Parses the VCF extra info that is usually at the top of the file"""
        return self.meta.add_extra(entry)

    def add_info(self, entry):
        """Parses the VCF Info field and returns a VCFMeta object"""
        return self.meta.add_info(entry)

    def add_filter(self, entry):
        """Parses the VCF Filter field and returns a VCFMeta object"""
        return self.meta.add_filter(entry)

    def add_format(self, entry):
        """Parses the VCF Format field and returns a VCFMeta object"""
        return self.meta.add_format(entry)

    def add_header(self, entry):
        """Parses the VCF Header field and returns the number of samples in the VCF file"""
        info = entry.split('\t')
        self.n_individuals = len(info)-9
        for i,v in enumerate(info[9:]):
            self.individuals[v] = i
        return self.n_individuals > 0

    def add_contig(self, entry):
        """Not Implemented"""
        return True

    def add_alt(self, entry):
        """This adds an alternative allele as specified by the VCF meta information and returns a VCFMeta object"""
        return self.meta.add_alt(entry)

    def parse_entry(self, row):
        """Parse an individual VCF entry and return a VCFEntry which contains information about
        the call (such as alternative allele, zygosity, etc.)

        """
        var_call = VCFEntry(self.individuals)
        var_call.parse_entry(row)
        return var_call

    def add_entry(self, row):
        """This will parse the VCF entry and also store it within the VCFFile. It will also
        return the VCFEntry as well.

        """
        var_call = VCFEntry(self.individuals)
        var_call.parse_entry( row )
        self.entries[(var_call.chrom, var_call.pos)] = var_call
        return var_call

    def get_header(self, individual=-1):
        """Returns the vcf header

        """
        type_map = dict([(val,key) for key,val in self.meta.type_map.iteritems()])
        extra = '\n'.join(['##{0}'.format(i) for i in self.meta.extra])
        info = '\n'.join(['##INFO=<ID={0},Number={1},Type={2},Description={3}>'.format(key, val.get('num_entries','.'), type_map.get(val.get('type', '')), val.get('description')) for key,val in self.meta.info.iteritems()])
        filter = '\n'.join(['##FILTER=<ID={0},Description={1}>'.format(key, val.get('description','.')) for key,val in self.meta.filter.iteritems()])
        format = '\n'.join(['##FORMAT=<ID={0},Number={1},Type={2},Description={3}>'.format(key, val.get('num_entries','.'), type_map.get(val.get('type', '')), val.get('description')) for key,val in self.meta.format.iteritems()])
        alt = '\n'.join(['##ALT=<ID={0},Description={1}>'.format(key, val.get('description','.')) for key,val in self.meta.alt.iteritems()])
        header = '\t'.join(['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT'])
        if individual is not None:
            if individual == -1:
                individual = '\t'.join(self.individuals.keys())
            else:
                if isinstance(individual, int):
                    for i, v in self.individuals.iteritems():
                        if v == individual:
                            individual = i
                            break
            header += '\t'+individual
        return '\n'.join([extra, info, filter, format, alt, header])


    def get_individual(self, individual=None):
        """Returns the index of the individual
        """

        return self.individuals.get(individual, None)


class VCFMeta(object):

    def __init__(self):
        """VCF Metainfo Container
        
        This contains information in the metainfo of a VCF File
        such as the format fields, ids, and filters.
    
        """
        self.info = OrderedDict()
        self.filter = OrderedDict()
        self.format = OrderedDict()
        self.alt = OrderedDict()
        self.extra = []
        self.type_map = {'Integer': int, 'Float': float, 'Numeric': float, 'Character': str, 'String': str, 'Flag': int}

    def add_extra(self, entry):
        self.extra.append(entry[2:])
        return True

    def add_info(self, entry):
        """Parse and store the info field"""
        entry = entry[8:-1]
        info = entry.split(',')
        if len(info) < 4:
            return False
        for v in info:
            key, value = v.split('=', 1)
            if key == 'ID':
                self.info[value] = {}
                id_ = value
            elif key == 'Number':
                if value == 'A' or value == 'G':
                    value = -1
                self.info[id_]['num_entries'] = value
            elif key == 'Type':
                self.info[id_]['type'] = self.type_map[value]
            elif key == 'Description':
                self.info[id_]['description'] = value
                if len(info) > 4:
                    self.info[id_]['description'] += '; '.join(info[4:])
                break
        return True

    def add_filter(self, entry):
        """Parse and store the filter field"""
        entry = entry[10:-1]
        info = entry.split(',')
        if len(info) < 2:
            return False
        for v in info:
            key, value = v.split('=', 1)
            if key == 'ID':
                self.filter[value] = {}
                id_ = value
            elif key == 'Description':
                self.filter[id_]['description'] = value
                if len(info) > 2:
                    self.info[id_]['description'] += '; '.join(info[2:])
        return True

    def add_format(self, entry):
        """Parse and store the format field"""
        entry = entry[10:-1]
        info = entry.split(',')
        if len(info) < 4:
            return False
        for v in info:
            key, value = v.split('=', 1)
            if key == 'ID':
                self.format[value] = {}
                id_ = value
            elif key == 'Number':
                if value == 'A' or value == 'G':
                    value = -1
                self.format[id_]['num_entries'] = value
            elif key == 'Type':
                self.format[id_]['type'] = self.type_map[value]
            elif key == 'Description':
                self.format[id_]['description'] = value
                if len(info) > 4:
                    self.format[id_]['description'] += '; '.join(info[4:])
                break
        return True

    def add_alt(self, entry):
        """Parse and store the alternative allele field"""
        entry = entry[7:-1]
        info = entry.split(',')
        if len(info) < 2:
            return False
        for v in info:
            key, value = v.split('=', 1)
            if key == 'ID':
                self.alt[value] = {}
                id_ = value
            elif key == 'Description':
                self.alt[id_]['description'] = value
                if len(info) > 4:
                    self.alt[id_]['description'] += '; '.join(info[4:])
                break
        return True

class VCFEntry(object):
    def __init__(self, individuals=None):
        """The VCFEntry object holds our sample information such as the
        reference genotype, the alternative genotypes, phase information, etc.

        """
        # These are all indexed by the samples index so we can do out of order
        # indexing without having to build empty lists
        self.phase = {}
        self.genotype = OrderedDict()
        self.genome_quality = {}
        self.depth = {}
        self.individuals = individuals
        self.gfilter = ""
        self.alt = None
        self.GT = -1
        self.GQ = -1
        self.DP = -1
        self.FT = -1
        self.ploidy = 2

    def __str__(self):
        """Returns the VCF entry as it appears in the vcf file minus sample info"""
        return '\t'.join([self.chrom, self.pos, self.id, '.' if not self.ref else self.ref,
                          '.' if not self.alt else ','.join(self.alt), self.qual,
                          'PASS' if self.passed else ';'.join(self.passed),
                          self.info,
                          ':'.join(self.format)])

    def sample_string(self, individual=-1):
        """Returns the VCF entry as it appears in the vcf file"""
        base = str(self)
        extra = self.get_sample_info(individual=individual)
        extra = [':'.join([str(j) for j in i]) for i in zip(*extra.values())]
        return '\t'.join([base, '\t'.join(extra)])

    def get_sample_info(self, individual=-1):
        """Returns the sample info of a given sample or all by default

        """
        if isinstance(individual, str):
            individual = self.individuals[individual]
        extra = OrderedDict()
        for format_ in self.format:
            index = getattr(self, format_)
            if index != -1:
                if format_ == 'GT':
                    d = self.genotype
                elif format_ == 'GQ':
                    d = self.genome_quality
                elif format_ == 'DP':
                    d = self.depth
                if individual == -1:
                    if len(d) != len(self.samples):
                        [self.parse_sample(i) for i in xrange(len(self.samples))]
                    extra[format_] = [d[i] for i in xrange(len(d))]
                else:
                    if individual not in d:
                        self.parse_sample(individual)
                    extra[format_] = [d[individual]]
        return extra

    def is_homozygous(self, individual=None):
        """This will give a boolean list corresponding to whether each individual
        is homozygous for the alternative allele.

        """
        if individual is not None:
            if isinstance(individual, str):
                individual = self.individuals[individual]
            alts = self.genotype[individual]
            return [sum(alts) == len(alts)] if sum(alts) > 0 else [False]
        else:
            return [sum(alts) == len(alts) if sum(alts) > 0 else False for i, alts in self.genotype.iteritems()]

    def is_heterozygous(self, individual=None):
        """This will give a boolean list corresponding to whether each individual
        is heterozygous for the alternative allele.

        """
        if individual is not None:
            if isinstance(individual, str):
                individual = self.individuals[individual]
            alts = self.genotype[individual]
            return [sum(alts) != len(alts)] if sum(alts) > 0 else [False]
        else:
            return [sum(alts) != len(alts) if sum(alts) > 0 else False for i, alts in self.genotype.iteritems()]


    def get_alt(self, individual=0, nucleotides_only=True):
        """Returns the alternative alleles of the individual as a list"""
        #not i.startswith(',') is put in to handle cases like <DEL:ME:ALU> where we have no alternate allele
        #but some reference
        if isinstance(individual, str):
            individual = self.individuals[individual]
        if nucleotides_only:
            return [self.alt[i-1].replace('.', '') for i in self.genotype[individual] if i > 0 and not self.alt[i-1].startswith('<')]
        else:
            return [self.alt[i-1].replace('.', '') for i in self.genotype[individual] if i > 0]

    def get_alt_length(self, individual=0):
        """Returns the number of basepairs of each alternative allele"""
        if isinstance(individual, str):
            individual = self.individuals[individual]
        return [len(self.alt[i-1].replace('.','')) for i in self.genotype[individual] if i > 0 and not self.alt[i-1].startswith('<')]

    def get_alt_lengths(self):
        """Returns the longest length of the variant. For deletions, return is negative,
        SNPs return 0, and insertions are +. None return corresponds to no variant in interval
        for specified individual

        """
        #this is a hack to store the # of individuals without having to actually store it
        out = []
        for i in xrange(len(self.genotype)):
            valid_alt = self.get_alt_length(individual=i)
            if not valid_alt:
                out.append(None)
            else:
                out.append(max(valid_alt)-len(self.ref))
        return out

    def has_snp(self, individual=0):
        """Returns a boolean list of SNP status, ordered by samples"""
        if isinstance(individual, str):
            individual = self.individuals[individual]
        alts = self.get_alt(individual=individual)
        if alts:
            return [i != self.ref and len(i) == len(self.ref) for i in alts]
        return [False]

    def has_insertion(self, individual=0):
        """Returns a boolean list of insertion status, ordered by samples"""
        if isinstance(individual, str):
            individual = self.individuals[individual]
        alts = self.get_alt(individual=individual)
        if alts:
            return [i != self.ref and len(i) > len(self.ref) for i in alts]
        return [False]

    def has_deletion(self, individual=0):
        """Returns a boolean list of deletion status, ordered by samples"""
        if isinstance(individual, str):
            individual = self.individuals[individual]
        alts = self.get_alt(individual=individual)
        if alts:
            return [i != self.ref and len(i) < len(self.ref) for i in alts]
        return [False]

    def has_variant(self, individual=0):
        if isinstance(individual, str):
            individual = self.individuals[individual]
        if individual not in self.genotype:
            self.parse_sample(individual)
        gt = self.genotype[individual]
        return not gt.count('0') == len(gt)-1

    def parse_entry(self, entry):
        """This parses a VCF row and stores the relevant information"""
        entry = entry.split('\t')
        self.chrom, self.pos, self.id, self.ref, alt_, self.qual, filter_, info, self.format = entry[:9]
        self.samples = entry[9:]
        self.alt = alt_.split(',')
        if filter_ == 'PASS' or filter_ == '.':
            self.passed = True
        else:
            self.passed = filter_.split(';')
        self.info = info
        # currently unused
        #if info != '.':
            #info_l = info.split(';')
            #self.info = [v.split('=') if '=' in v else (v,1) for v in info_l]
        self.format = self.format.split(':')
        if 'GT' in self.format:
            self.GT = self.format.index('GT')
        if 'GQ' in self.format:
            self.GQ = self.format.index('GQ')
        if 'DP' in self.format:
            self.DP = self.format.index('DP')
        if 'FT' in self.format:
            self.FT = self.format.index('FT')

    def parse_sample(self, individual=-1):
        #on to the samples
        self.extra_sample_info = {}
        if isinstance(individual, str):
            individual = self.individuals[individual]
        sample = self.samples[individual]
        for info_n, sample_info in enumerate(sample.split(':')):
            if info_n == self.GT:
                if len(sample_info) == 3 and ('|' in sample_info or '/' in sample_info):
                    self.genotype[individual] = [int(i) for i in GENOTYPE_SPLIT.split(sample_info)]
            elif info_n == self.GQ:
                if not sample_info:
                    self.genome_quality[individual] = None
                elif sample_info == '.':
                    sample_info = -1
                else:
                    sample_info = int(sample_info)
                if sample_info > 99:
                    sample_info = 99
                self.genome_quality[individual] = sample_info
            elif info_n == self.DP:
                if not sample_info or sample_info == '.':
                    self.depth[individual] = -1
                else:
                    self.depth[individual] = int(sample_info)
            elif info_n == self.FT:
                #not supported, I haven't encountered this yet
                pass
            else:
                try:
                    self.extra_sample_info[individual][info_n] = sample_info
                except KeyError:
                    self.extra_sample_info[individual] = {info_n: sample_info}

class SequenceObject(object):
    chrom = None
    start = None
    end = None

    def __len__(self):
        return self.end-self.start


class BedObject(SequenceObject):
    name = None
    strand = None

    def parse(self, entry):
        if not isinstance(entry, list):
            entry = entry.split('\t')
        for order, attribute in zip(config.BED_ORDER, entry):
            setattr(self, order, attribute)


class GFFObject(SequenceObject):
    def __init__(self, info_delimiter=';', key_delimiter='=', quotechar='', attribute_delimiter=','):
        self.info_delimiter = info_delimiter
        self.key_delimiter = key_delimiter
        self.quotechar = quotechar
        #the attribute_delimiter is set but not used generally since commas are abused
        #and so few GFF/GFF3/GTF files actually follow the specification.
        #For parent/child relationships this is used however.
        self.attribute_delimiter = attribute_delimiter

    def parse_entry(self, row):
        self.seqid, self.source, self.feature_type, self.start, self.end, \
        self.score, self.strand, self.phase, info = row.split('\t')
        self.attributes = {}
        self.start = int(self.start)
        self.end = int(self.end)
        for entry in info.split(self.info_delimiter):
            entry = entry.strip()
            if not entry:
                continue
            try:
                key,value = entry.split(self.key_delimiter)
            except ValueError:
                sys.stderr.write('Error in info field on entry %s\n' % entry)
                continue
            if self.quotechar:
                if value.startswith(self.quotechar) and value.endswith(self.quotechar):
                    value = value[len(self.quotechar):-1*len(self.quotechar)]
            self.attributes[key] = value

    def __str__(self):
        attributes = self.info_delimiter.join(['%s%s%s%s%s' % (key,self.key_delimiter,self.quotechar,value,self.quotechar)
                               for key,value in self.attributes.iteritems()])
        general = '\t'.join([self.seqid, self.source, self.feature_type, str(self.start), str(self.end), self.score,
                             self.strand, self.phase, attributes])
        return general


class GFFFeature(object):
    def __init__(self, id):
        self.features = set([])
        self.id = id

    def add_child(self, child):
        """Children are GFFFeatures and are defined when added. This is done to avoid memory overheads
        that may be incurred by GFF files that have millions of rows.

        """
        child_id = getattr(child, 'id', None)
        if child_id:
            if not hasattr(self, 'children'):
                self.children = {}
            if child_id not in self.children:
                self.children[child_id] = child

    def get_children(self):
        """A dictionary of GFF features which are children of this feature

        :return: A dictionary of children GFF Features
        """
        return {} if not hasattr(self, 'children') else self.children

    def parts(self):
        """A generator consisting GFFObjects within this feature

        """
        for i in self.features:
            yield i

    def __str__(self):
        rep = [str(i) for i in self.parts()]
        return '\n'.join(rep)