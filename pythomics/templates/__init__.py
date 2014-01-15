import gzip

class GenericIterator(object):
    def __init__(self, filename, **kwrds):
        if not filename:
            print 'error'
            raise KeyError
        if isinstance(filename, (str, unicode)) and filename.endswith('.gz'):
            self.filename = gzip.GzipFile(filename)
        elif isinstance(filename, (str, unicode)):
            self.filename = open(filename)
        elif isinstance(filename, (file)):
            self.filename = filename
        else:
            print 'no idea what we got'
            raise KeyError
        
    def __iter__(self):
        return self
    
    def next(self):
        raise StopIteration