import collections
import re

import six

from . import config


class Enzyme(object):
    def __init__(self, enzyme="trypsin", pattern=None):
        if pattern:
            rules = pattern.split(",")
        else:
            rules = config.ENZYMES[enzyme]
            rules = rules.split(",") if rules else [None]
        self.patterns = []
        exclude = re.compile(r"\{([^\}])")
        for rule in rules:
            if rule:
                cterm, nterm = rule.split("|")
                cterm = cterm.replace("X", "A-Z")
                nterm = nterm.replace("X", "A-Z")
                no_cleave = exclude.search(nterm)
                if no_cleave:
                    # incorporate negative lookhead
                    regex = re.compile("(%s(?!%s))" % (cterm, no_cleave.group(1)))
                else:
                    regex = re.compile("(%s)" % rule)  # it's already in regex form
            else:
                regex = None
            self.patterns.append(regex)

    def cleave(self, sequence, min=7, max=30, unique=False):
        peptides = []
        for regex in self.patterns:
            if regex:
                try:
                    digests = [i for i in regex.split(sequence)]
                except ValueError:
                    digests = [sequence]
            else:
                digests = [sequence]
            digests = [
                "".join(i)
                for i in six.moves.zip_longest(
                    digests[0::2], digests[1::2], fillvalue=""
                )
            ]
            if unique:
                peptides += list(
                    collections.OrderedDict.fromkeys(
                        [i for i in digests if min <= len(i) <= max]
                    )
                )
            else:
                peptides += [i for i in digests if min <= len(i) <= max]
        return peptides
