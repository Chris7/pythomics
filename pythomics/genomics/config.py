__author__ = 'Chris Mitchell'

GFF_TAGS = set(['ID', 'Name', 'Alias', 'Parent', 'Target', 'Gap',
                'Derives_from', 'Note', 'Dbxref', 'Ontology_term',
                'Is_circular'])


BED_ORDER = ['chrom', 'start', 'end', 'name', 'score', 'strand',]