"""
Author: Chris Mitchell (chris.mit7@gmail.com)
Copyright (C) 2012 Chris Mitchell

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""
__author__ = 'chris'
from pythomics.proteomics import config

class Chromatogram(object):
    pass


class ScanObject(object):
    """
    A scan object to store peaklist information in.
    Attributes:
    title, charge, mass, scans(list), rt
    """
    def __init__(self):
        self.scans = []
        self.rt = 0.0
        self.ms_level = None

    def writeScan(self, o):
        o.write('BEGIN IONS\n')
        o.write('TITLE=%s\n'%self.title)
        try:
            o.write('RTINSECONDS=%f\n'%self.rt)
        except AttributeError:
            pass
        o.write('PEPMASS=%s\n'%self.mass)
        o.write('CHARGE=%s\n'%self.charge)
        for i in self.scans:
            o.write('%f\t%f\n'%i)
        o.write('END IONS\n\n')

class PeptideObject(ScanObject):
    """
    An enhanced scan object that can store peptide information as well
    attributes:
    mods (set item), peptide, expect, id, acc(accession)
    matched -> dict with keys and lists as values:
    keys: labels(like y1+-h20), m/z, intensity, error, series(y,a,b,...), start, end, losses(neutral losses), charge
    for msf files: spectrumId, confidence, rank
    """
    def __init__(self):
        super(PeptideObject, self).__init__()
        self.mods = set([])
        self.peptide = ""
        self.hit = 0

    def addModification(self, aa,position, modMass, modType):
        """
        !!!!MODIFICATION POSITION IS 0 BASED!!!!!!
        Modifications are stored internally as a tuple with this format:
        (amino acid modified, index in peptide of amino acid, modification type, modification mass)
        ie (M, 7, Oxidation, 15.9...)
        such as: M35(o) for an oxidized methionine at residue 35
        """
        #clean up xtandem
        if not modType:
            #try to figure out what it is
            tmass = abs(modMass)
            smass = str(tmass)
            prec = len(str(tmass-int(tmass)))-2
            precFormat = '%'+'0.%df'%prec
            modType = ""
            masses = config.MODIFICATION_MASSES
            for i in masses:
                if tmass in masses[i] or smass == precFormat%masses[i][0]:
                    #found it
                    modType = i
            if not modType:
                sys.stderr.write('mod not found %s\n'%modMass)
        self.mods.add((aa,str(position),str(modMass),str(modType)))

    def getModifications(self):
        return '|'.join([','.join(i) for i in self.mods])
