class ColumnFunctions(object):
    strict = False
    METHODS = ['concat', 'mean', 'median', 'var', 'std', 'sum', 'count']

    def __init__(self, parser_args):
        if parser_args.strict:
            self.strict = True

    def process_list(self, l, ordered=False):
        if ordered and isinstance(l, set):
            l = list(l)
        if self.strict:
            try:
                l = [float(i) for i in l if i != '']
            except ValueError:
                import sys
                sys.stderr.write('Invalid entry found in list %s.\n'%l)
                raise ValueError
        else:
            nl = []
            for i in l:
                try:
                    nl.append(float(i))
                except ValueError:
                    pass
            l = nl
        return l

    def concat(self, l):
        if not len(l):
            return 'NA'
        return ';'.join([str(i) for i in l])

    def mean(self, l):
        if not len(l):
            return 'NA'
        l = self.process_list(l)
        return sum(l)/float(len(l))

    def median(self, l):
        if not len(l):
            return 'NA'
        elif len(l) == 1:
            if isinstance(l, set):
                return iter(l).next()
            else:
                return l[0]
        l = self.process_list(l, ordered=True)
        l_sorted = sorted(l)
        mid = len(l_sorted)/2
        if len(l_sorted)%2:
            # odd
            return l_sorted[mid]
        else:
            return self.mean(l_sorted[mid-1:mid+1])

    def var(self, l):
        if not len(l):
            return 'NA'
        elif len(l) == 1:
            return 0
        l = self.process_list(l)
        sample_mean = self.mean(l)
        dev = 0
        for x in l:
            dev += (x - sample_mean)*(x - sample_mean)
        return dev/(len(l)-1)

    def std(self, l):
        if not len(l):
            return 'NA'
        elif len(l) == 1:
            return 0
        l = self.process_list(l)
        import math
        return math.sqrt(self.var(l))

    def sum(self, l):
        if not len(l):
            return 'NA'
        l = self.process_list(l)
        return sum(l)

    def count(self, l):
        l = self.process_list(l)
        return len(l)
