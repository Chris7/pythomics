class VCFFile(object):
    """
    This is designed to be somewhat consistent with the C++ vcftools
    """
    def __init__(self, filename):
        self.filename = filename
        self.meta = VCFMeta()
        self.n_individuals = 0
        self.individuals = {}
        self.entries = {}
        
    def add_info(self, entry):
        return self.meta.add_info(entry)
        
    def add_filter(self, entry):
        return self.meta.add_filter(entry)
        
    def add_format(self, entry):
        return self.meta.add_format(entry)
    
    def add_header(self, entry):
        info = entry.split('\t')
        self.n_individuals = len(info)-9
        for i in xrange(self.n_individuals):
            self.individuals[i] = [] 
        return self.n_individuals > 0
    
    def add_contig(self, entry):
        return True
    
    def add_entry(self, entry):
        var_call = VCFEntry(self.n_individuals)
        var_call.parse_entry( entry )
        self.entries[(var_call.chrom, var_call.pos)] = var_call
        return self.entries[(var_call.chrom, var_call.pos)]
    
class VCFMeta(object):
    def __init__(self):
        self.info = {}
        self.filter = {}
        self.format = {}
        self.type_map = {'Integer': int, 'Float': float, 'Numeric': float, 'Character': str, 'String': str, 'Flag': int}
        
    def add_info(self, entry):
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
    
class VCFEntry(object):
    def __init__(self, individuals):
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
    
    def parse_entry(self, entry):
        entry = entry.split('\t')
        self.chrom, self.pos, self.id, self.ref, self.alt, self.qual, filter_, info, self.format = entry[:9]
        samples = entry[9:]
        if self.alt == '.':
            self.alt = None
        if filter_ == 'PASS':
            self.passed = True
        else:
            if filter_ != '.':
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