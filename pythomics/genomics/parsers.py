__author__ = 'Chris Mitchell'

import sys, copy
import pythomics.templates as templates
import pythomics.genomics.structures as structure
import pythomics.genomics.config as config

class VCFMixin(object):
    vcf_file = None

    def get_header(self, individual=-1):
        return self.vcf_file.get_header(individual=individual)

    def get_individual(self, individual=None):
        if individual is None:
            return self.vcf_file.individuals.keys()
        return self.vcf_file.get_individual(individual=individual)


class VCFIterator(templates.GenericIterator, VCFMixin):
    
    def __init__(self, filename, sample=None):
        """An iterator over the VCF file format
        
        This reads VCF files and has been tested on VCF 4.0.
        The returned items are VCFEntry
        
        """
        super(VCFIterator, self).__init__(filename)
        #process our meta information
        row = self.filename.next().strip()
        self.vcf_file = structure.VCFFile(filename)
        self.inum=0
        while row[:2] == '##':
            if row.startswith('##INFO='):
                assert(self.vcf_file.add_info(row)==True)
            elif row.startswith('##FILTER='):
                assert(self.vcf_file.add_filter(row))
            elif row.startswith('##FORMAT='):
                assert(self.vcf_file.add_format(row))
            elif row.startswith('##CONTIG='):
                assert(self.vcf_file.add_contig(row))
            elif row.startswith('##ALT='):
                assert(self.vcf_file.add_alt(row))
            elif row.startswith('##'):
                assert(self.vcf_file.add_extra(row))
            row = self.filename.next().strip()
        #got to the end of meta information, we now have the header
        assert(self.vcf_file.add_header(row))
        self.sample = sample
            
    def __iter__(self):
        return self

    def get_vcf_entry(self):
        row = self.filename.next().strip()
        while not row or row.startswith('#'):
            row = self.filename.next().strip()
        self.inum+=1
        if self.inum%100000 == 0:
            sys.stderr.write('Processed %d VCF entries\n' % self.inum)
        vcf_entry = self.vcf_file.parse_entry(row)
        if self.sample is not None:
            vcf_entry.parse_sample(individual=self.sample)
        return vcf_entry
    
    def next(self):
        vcf_entry = self.get_vcf_entry()
        if self.sample is not None:
            while vcf_entry and not vcf_entry.has_variant(self.sample):
                vcf_entry = self.get_vcf_entry()
        return vcf_entry


class GFFIterator(templates.GenericIterator):
    def __init__(self, filename, info_delimiter=';', key_delimiter='=', quotechar=''):
        super(GFFIterator, self).__init__(filename)
        self.info_delimiter=info_delimiter
        self.key_delimiter=key_delimiter
        self.quotechar=quotechar

    def __iter__(self):
        return self

    def next(self):
        row = self.filename.next()
        while not row or row[0] == '#':#skip blanks and comments
            row = self.filename.next()
        ob = structure.GFFObject(info_delimiter=self.info_delimiter,
                                 key_delimiter=self.key_delimiter,
                                 quotechar=self.quotechar)
        ob.parse_entry(row)
        return ob


class BedIterator(templates.GenericIterator):
    def __init__(self, filename):
        super(BedIterator, self).__init__(filename)

    def __iter__(self):
        return self

    def next(self):
        row = self.filename.next()
        while not row:#skip blanks
            row = self.filename.next()
        ob = structure.BedObject()
        ob.parse(row)
        return ob


class BamIterator(templates.GenericIterator):
    pass


class SamIterator(templates.GenericIterator):
    pass


class PileupIterator(templates.GenericIterator):
    pass


class GFFReader(templates.GenericIterator):
    def __init__(self, filename, fast_filter=None, info_delimiter=';', key_delimiter='=', quotechar='',
                 tag_map={'ID': 'ID', 'Parent': 'Parent'}, preset=None):
        """This will read an entire GFF/GTF/GFF3 file and establish parent-child relations and
         allow for querying based on attributes such as transcripts, gene ids, and allow writing
         of useful file types such as fasta files of coding sequences, etc. It works at the
         "feature" level, which allows for grouping of exons of a gene and so on.

         The implementation is liberal because GFF/GFF3/GTF standards are so poorly followed. So
         errors such as a missing parent feature will be reported, but not break usage.

        :param filename: GFF file to read
        :param fast_filter: A list of attributes to keep in memory a reference to for fast querying
        :param info_delimiter: The character to split the info field on (default: ";")
        :param key_delimiter: The character to split the key-value pairings of the info field on (default: "=")
        :param quotechar: The character surrounding values in the info field (default: None)
        :param tag_map: A dictionary which can substitute the ID and Parent tags for alternatives
        :param preset: A string identifying various presets (choices: cufflinks)
        """
        super(GFFReader, self).__init__(filename)
        self.positions = {}
        self.fast_attributes = {}
        self.filters = set(fast_filter) if fast_filter else []
        GFF_TAGS = copy.deepcopy(config.GFF_TAGS)
        self.feature_map = {}
        child_list = []
        if preset and preset.lower() == 'cufflinks':
            key_delimiter=' '
            quotechar='"'
            self.id_tag = 'transcript_id'
            self.parent_tag = 'gene_id'#'transcript_id'#it's weird, but we're our own parent
        else:
            self.id_tag = tag_map['ID']
            self.parent_tag = tag_map['Parent']
        GFF_TAGS.add(self.id_tag)
        GFF_TAGS.add(self.parent_tag)
        for gff_object in GFFIterator(filename, info_delimiter=info_delimiter,
                                      key_delimiter=key_delimiter, quotechar=quotechar):
            chrom, start, end = gff_object.seqid, int(gff_object.start), int(gff_object.end)
            try:
                self.positions[chrom][(start, end)].append(gff_object)
            except KeyError:
                try:
                    self.positions[chrom][(start, end)] = [gff_object]
                except KeyError:
                    self.positions[chrom] = {(start, end): [gff_object]}
            fast_lookup = [fast_access for fast_access in self.filters for attribute in gff_object.attributes if fast_access in attribute]
            for fast_access in fast_lookup:
                try:
                    self.fast_attributes[fast_access].append(gff_object)
                except KeyError:
                    self.fast_attributes[fast_access] = [gff_object]
            known_tags = set([attribute for attribute in gff_object.attributes if attribute in GFF_TAGS])
            if self.id_tag in known_tags:
                gff_id = gff_object.attributes.get(self.id_tag, None)
                if gff_id:
                    try:
                        gff_feature = self.feature_map[gff_id]
                    except KeyError:
                        gff_feature = structure.GFFFeature(gff_id)
                        self.feature_map[gff_id] = gff_feature
                    gff_feature.features.add(gff_object)
                    #This is kept in the known_tags valid ID block, if we don't use IDs we shouldn't
                    #be using parent-child relations either
                    if self.parent_tag in known_tags:
                        child_list.append(gff_object)
            else:
                #this is a poorly formed GFF file and feature relations cannot be clearly established
                gff_id = str(len(self.feature_map))
                gff_feature = structure.GFFFeature(gff_id)
                gff_feature.features.add(gff_object)
                self.feature_map[gff_id] = gff_feature
        #Add parent-child relationships here since all GFF objects will be made at this point
        for child_gff_object in child_list:
            #get parents from the gff entries
            parent_ids = child_gff_object.attributes[self.parent_tag].split(',')
            for parent_id in parent_ids:
                parent_feature = self.feature_map.get(parent_id, None)
                if not parent_feature:
                    #missing the parent feature, make it but complain about it
                    parent_feature = structure.GFFFeature(parent_id)
                    self.feature_map[parent_id] = parent_feature
                    sys.stderr.write('Missing Parent Feature: %s\n' % parent_id)
                parent_feature.add_child(self.feature_map.get(child_gff_object.attributes.get(self.id_tag, None)))


    def get_attribute(self, attribute, value=None, features=False):
        """This returns a list of GFF objects (or GFF Features) with the given attribute and if supplied, those
        attributes with the specified value

        :param attribute: The 'info' field attribute we are querying
        :param value: Optional keyword, only return attributes equal to this value
        :param features: Optional keyword, return GFF Features instead of GFF Objects
        :return: A list of GFF objects (or GFF features if requested)
        """
        if attribute in self.filters:
            valid_gff_objects = self.fast_attributes[attribute] if not value else\
                [i for i in self.fast_attributes[attribute] if i.attributes.get(attribute, False) == value]
            if features:
                valid_ids = [gff_object.attributes.get(self.id_tag, None) for gff_object in valid_gff_objects]
                return [self.feature_map[gff_id] for gff_id in valid_ids if gff_id]
            else:
                return valid_gff_objects
        else:
            valid_gff_objects = [gff_object for gff_feature in self.feature_map.values()
                              for gff_object in gff_feature.features
                              if gff_object.attributes.get(attribute, False)]
            valid_gff_objects = valid_gff_objects if not value else [gff_object for gff_object in valid_gff_objects
                                                                     if gff_object.attributes[attribute] == value]
            if features:
                valid_ids = [gff_object.attributes.get(self.id_tag, None) for gff_object in valid_gff_objects]
                return [self.feature_map[gff_id] for gff_id in valid_ids if gff_id]
            else:
                return valid_gff_objects

    def contains(self, seqid, start, end, overlap=True):
        """This returns a list of GFF objects which cover a specified location.

        :param seqid: The landmark identifier (usually a chromosome)
        :param start: The 1-based position of the start of the range we are querying
        :param end: The 1-based position of the end of the range we are querying
        :param overlap: A boolean value, if true we allow features to overlap the query range.
        For instance, overlap=True with the range (5,10), will return a GFF object
        spanning from (8,15). overlap=False will only return objects fully containing the range.
        :return: A list of GFF objects
        """
        d = self.positions.get(seqid,[])
        if overlap:
            return [gff_object for gff_start, gff_end in d
                    for gff_object in d[(gff_start, gff_end)]
                    if not (end <= gff_start or start >= gff_end)]
        else:
            return [gff_object for gff_start, gff_end in d
                    for gff_object in d[(gff_start, gff_end)]
                    if (gff_start <= start and gff_end >= end)]

    def get_feature(self, feature_id):
        """This returns a GFF Feature object corresponding to the
        provided id.

        :param feature_id: The ID of the feature we wish to fetch
        :return: A GFF Feature object or None if it does not exist.
        """
        return self.feature_map.get(feature_id)

    def get_features(self):
        return self.feature_map.iteritems()

class VCFReader(templates.GenericIterator, VCFMixin):
    def __init__(self, filename, append_chromosome=False, store_positions=True, sample=None):
        """This will read an entire VCF file and allow for querying based on attributes such
         location, individuals.

        :param filename: GFF file to read

        """
        super(VCFReader, self).__init__(filename)
        self.positions = {}
        self.append_chromosome = append_chromosome
        self.vcf_iterator = VCFIterator(filename, sample=sample)
        self.vcf_file = self.vcf_iterator.vcf_file
        self.vcf_entries = []
        for vcf_entry in self.vcf_iterator:
            chrom, start = vcf_entry.chrom, vcf_entry.pos
            if append_chromosome:
                chrom = 'chr%s' % chrom
            if store_positions:
                start = int(start)
                furthest_end = max([abs(i) if i is not None else None for i in vcf_entry.get_alt_lengths()])
                if furthest_end is not None:
                    end = start+abs(furthest_end)
                    try:
                        self.positions[chrom][(start, end)].append(vcf_entry)
                    except KeyError:
                        try:
                            self.positions[chrom][(start, end)] = [vcf_entry]
                        except KeyError:
                            self.positions[chrom] = {(start, end): [vcf_entry]}
            self.vcf_entries.append(vcf_entry)

    def get_sample(self, individual=0):
        if isinstance(individual, str):
            individual = self.get_individual(individual=individual)
        for entry in self.vcf_entries:
            if entry.has_variant(individual=individual):
                yield entry


    def contains(self, chrom, start, end, overlap=True):
        """This returns a list of VCFEntry objects which cover a specified location.

        :param chrom: The landmark identifier (usually a chromosome)
        :param start: The 1-based position of the start of the range we are querying
        :param end: The 1-based position of the end of the range we are querying
        :param overlap: A boolean value, if true we allow features to overlap the query range.
        For instance, overlap=True with the range (5,10), will return a VCFEntry object
        spanning from (8,15). overlap=False will only return objects fully containing the range.
        :return: A list of VCFEntry objects
        """
        d = self.positions.get(chrom,[])
        if overlap:
            return [vcf_entry for vcf_start, vcf_end in d
                    for vcf_entry in d[(vcf_start, vcf_end)]
                    if not (end < vcf_start or start > vcf_end)]
        else:
            return [vcf_entry for vcf_start, vcf_end in d
                    for vcf_entry in d[(vcf_start, vcf_end)]
                    if (vcf_start <= start and vcf_end >= end)]

    def remove_variants(self, variants):
        """Remove a list of variants from the positions we are scanning"""
        chroms = set([i.chrom for i in variants])
        for chrom in chroms:
            if self.append_chromosome:
                chrom = 'chr%s' % chrom
            to_delete = [pos for pos in self.positions[chrom] if pos in variants]
            for pos in to_delete:
                del self.positions[chrom][pos]