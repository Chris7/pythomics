__author__ = 'Chris Mitchell'

import sys

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
        self.individuals = {}
        self.entries = {}

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
        for i in xrange(self.n_individuals):
            self.individuals[i] = []
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
        var_call = VCFEntry(self.n_individuals)
        var_call.parse_entry(row)
        return var_call

    def add_entry(self, row):
        """This will parse the VCF entry and also store it within the VCFFile. It will also
        return the VCFEntry as well.

        """
        var_call = VCFEntry(self.n_individuals)
        var_call.parse_entry( row )
        self.entries[(var_call.chrom, var_call.pos)] = var_call
        return var_call

class VCFMeta(object):

    def __init__(self):
        """VCF Metainfo Container
        
        This contains information in the metainfo of a VCF File
        such as the format fields, ids, and filters.
    
        """
        self.info = {}
        self.filter = {}
        self.format = {}
        self.type_map = {'Integer': int, 'Float': float, 'Numeric': float, 'Character': str, 'String': str, 'Flag': int}

    def add_info(self, entry):
        """Parse and store the info field"""
        entry = entry[8:-1]
        info = entry.split(',')
        if len(info) < 4:
            return False
        for v in info:
            key, value = v.split('=')
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
            key, value = v.split('=')
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
            key, value = v.split('=')
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
            key, value = v.split('=')
            if key == 'ID':
                self.format[value] = {}
                id_ = value
            elif key == 'Description':
                self.format[id_]['description'] = value
                if len(info) > 4:
                    self.format[id_]['description'] += '; '.join(info[4:])
                break
        return True

class VCFEntry(object):
    def __init__(self, individuals):
        """The VCFEntry object holds our sample information such as the
        reference genotype, the alternative genotypes, phase information, etc.

        """
        self.passed = [False for __ in xrange(individuals)]
        self.phase = ['/' for __ in xrange(individuals)]
        self.genotype = [[-1,-1]  for __ in xrange(individuals)]
        self.genome_quality = [-1 for __ in xrange(individuals)]
        self.depth = [-1 for __ in xrange(individuals)]
        self.gfilter = ""
        self.alt = None
        self.info = {}
        self.GT = -1
        self.GQ = -1
        self.DP = -1
        self.FT = -1
        self.ploidy = 2

    def __str__(self):
        """Returns the VCF entry as it appears in the vcf file, minus sample info currently"""
        info = self.info.keys()
        values = [self.info[i] for i in info]
#         if self.GT != -1:
#         if self.GQ != -1:
#         if self.DP != -1:
#         if self.FT != -1:
#         sample_info = self.genotype[n]
        return '\t'.join([self.chrom, self.pos, self.id, '.' if not self.ref else self.ref,
                          '.' if not self.alt else ','.join(self.alt), self.qual,
                          'PASS' if self.passed else ';'.join(self.passed),
                          ';'.join(['%s=%s' % (i,j) for i,j in zip(info,values)]),
                          ';'.join(self.format)])

    def is_homozygous(self, sample = None):
        """This will give a boolean list corresponding to whether each individual
        is homozygous for the alternative allele.

        """
        if sample:
            pass
        else:
            return [self.alt[i-1] == self.alt[j-1] if i > 0 and j > 0 else False for i,j in self.genotype]

    def is_heterozygous(self, sample = None):
        """This will give a boolean list corresponding to whether each individual
        is heterozygous for the alternative allele.

        """
        if sample:
            pass
        else:
            return [True if ((i == 0 and j > 0) or (i > 0 and j == 0)) else False for i,j in self.genotype]

    def get_alt(self, individual=0, nucleotides_only=True):
        """Returns the alternative alleles of the individual as a list"""
        #not i.startswith(',') is put in to handle cases like <DEL:ME:ALU> where we have no alternate allele
        #but some reference
        if nucleotides_only:
            return [self.alt[i-1].replace('.','') for i in self.genotype[individual] if i > 0 and not self.alt[i-1].startswith('<')]
        else:
            return [self.alt[i-1].replace('.','') for i in self.genotype[individual] if i > 0]

    def get_alt_length(self, individual=0):
        """Returns the number of basepairs of each alternative allele"""
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
        alts = self.get_alt(individual=individual)
        if alts:
            return [i != self.ref and len(i) == len(self.ref) for i in alts]
        return [False]

    def has_insertion(self, individual=0):
        """Returns a boolean list of insertion status, ordered by samples"""
        alts = self.get_alt(individual=individual)
        if alts:
            return [i != self.ref and len(i) > len(self.ref) for i in alts]
        return [False]

    def has_deletion(self, individual=0):
        """Returns a boolean list of deletion status, ordered by samples"""
        alts = self.get_alt(individual=individual)
        if alts:
            return [i != self.ref and len(i) < len(self.ref) for i in alts]
        return [False]

    def parse_entry(self, entry):
        """This parses a VCF row and stores the relevant information"""
        entry = entry.split('\t')
        self.chrom, self.pos, self.id, self.ref, alt_, self.qual, filter_, info, self.format = entry[:9]
        samples = entry[9:]
        self.alt = alt_.split(',')
        if filter_ == 'PASS' or filter_ == '.':
            self.passed = True
        else:
            self.passed = filter_.split(';')
        if info != '.':
            info_l = info.split(';')
            for v in info_l:
                if '=' in v:
                    key, value = v.split('=')
                else:
                    key, value = v, "1"
                self.info[key] = value
        self.format = self.format.split(':')
        if 'GT' in self.format:
            self.GT = self.format.index('GT')
        if 'GQ' in self.format:
            self.GQ = self.format.index('GQ')
        if 'DP' in self.format:
            self.DP = self.format.index('DP')
        if 'FT' in self.format:
            self.FT = self.format.index('FT')
        #on to the samples
        self.extra_sample_info = {}
        for n, sample in enumerate(samples):
            for info_n, sample_info in enumerate(sample.split(':')):
                if info_n == self.GT:
                    if len(sample_info) == 3 and (sample_info[1] == '|' or sample_info[1] == '/'):
                        self.phase[n] = sample_info[1]
                        self.genotype[n] = [int(i) if i != '.' else -1 for i in sample_info.split(sample_info[1])]
                elif info_n == self.GQ:
                    if not sample_info:
                        self.genome_quality[n] = None
                    elif sample_info == '.':
                        sample_info = -1
                    else:
                        sample_info = int(sample_info)
                    if sample_info > 99:
                        sample_info = 99
                    self.genome_quality[n] = sample_info
                elif info_n == self.DP:
                    if not sample_info or sample_info == '.':
                        self.depth[n] = -1
                    else:
                        self.depth[n] = int(sample_info)
                elif info_n == self.FT:
                    #not supported, I haven't encountered this yet
                    pass
                else:
                    self.extra_sample_info[info_n] = sample_info


class GFFObject(object):
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

    def __len__(self):
        return self.end-self.start

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