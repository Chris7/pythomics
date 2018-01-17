import six

from .. import templates

class FastqIterator(templates.GenericIterator):
    def _next(self):
        # return is sequence header, sequence, quality header, quality sequence
        _next = lambda: six.next(super(FastqIterator, self))
        return (_next(), _next(), _next(), _next())
