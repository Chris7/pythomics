__author__ = 'Chris Mitchell'

import sys
import pythomics.templates as templates
import pythomics.genomics.structures as structure
import pythomics.genomics.config as config


class VCFIterator(templates.GenericIterator):
    
    def __init__(self, filename):
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
            row = self.filename.next().strip()
        #got to the end of meta information, we now have the header
        assert(self.vcf_file.add_header(row))
            
    def __iter__(self):
        return self
    
    def next(self):
        row = self.filename.next()
        while not row:
            row = self.filename.next()
        self.inum+=1
        if self.inum%100000 == 0:
            sys.stderr.write('Processed %d VCF entries\n' % self.inum)
        return self.vcf_file.parse_entry( row.strip() )


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
    pass


class BamIterator(templates.GenericIterator):
    pass


class SamIterator(templates.GenericIterator):
    pass


class PileupIterator(templates.GenericIterator):
    pass


class GFFReader(templates.GenericIterator):
    def __init__(self, filename, fast_filter=None, info_delimiter=';', key_delimiter='=', quotechar='',
                 tag_map={'ID': 'ID', 'Parent': 'Parent'}):
        """This will read an entire GFF/GTF/GFF3 file and establish parent-child relations and
         allow for querying based on attributes such as transcripts, gene ids, and allow writing
         of useful file types such as fasta files of coding sequences, etc. It works at the
         "feature" level, which allows for grouping of exons of a gene and so on.

         The implementation is liberal because GFF/GFF3/GTF standards are so poorly followed. So
         errors such as a missing parent feature will be reported, but not break usage.

        """
        super(GFFReader, self).__init__(filename)
        self.positions = {}
        self.fast_attributes = {}
        self.filters = set(fast_filter) if fast_filter else []
        GFF_TAGS = config.GFF_TAGS
        id_tag = tag_map['ID']
        parent_tag = tag_map['Parent']
        self.feature_map = {}
        child_list = []
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
            if id_tag in known_tags:
                gff_id = gff_object.attributes.get('ID',None)
                if gff_id:
                    try:
                        gff_feature = self.feature_map[gff_id]
                    except KeyError:
                        gff_feature = structure.GFFFeature(gff_id)
                        gff_feature.features.add(gff_object)
                        self.feature_map[gff_id] = gff_feature
                    #This is kept in the known_tags valid ID block, if we don't use IDs we shouldn't
                    #be using parent-child relations either
                    if parent_tag in known_tags:
                        child_list.append(gff_feature)
            else:
                #this is a poorly formed GFF file and feature relations cannot be clearly established
                gff_id = str(len(self.feature_map))
                gff_feature = structure.GFFFeature(gff_id)
                gff_feature.features.add(gff_object)
                self.feature_map[gff_id] = gff_feature
        #Add parent-child relationships here since all GFF objects will be made at this point
        for child_gff_feature in child_list:
            #get parents from the gff entries
            parent_ids = [feature.attributes['Parent'] for feature in child_gff_feature.features
                          if feature.attributes.has_key('Parent')]
            parent_features = set([j for i in parent_ids for j in i.split(',')])
            for parent_feature_name in parent_features:
                parent_feature = self.feature_map.get(parent_feature_name, None)
                if parent_feature:
                    parent_feature.add_child(child_gff_feature)
                else:
                    sys.stderr.write('Missing Parent Feature: %s\n' % parent_feature_name)

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
                valid_ids = [gff_object.attributes.get('ID', None) for gff_object in valid_gff_objects]
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
                valid_ids = [gff_object.attributes.get('ID', None) for gff_object in valid_gff_objects]
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
                    if (gff_start >= start and gff_end >= end)]

    def get_feature(self, feature_id):
        """This returns a GFF Feature object corresponding to the
        provided id.

        :param feature_id: The ID of the feature we wish to fetch
        :return: A GFF Feature object or None if it does not exist.
        """
        return self.feature_map.get(feature_id, None)
