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
import sys

class Chromatogram(object):
    pass


class ScanObject(object):
    """
    A scan object to store peaklist information in.
    Attributes:
    title, charge, mass, scans(list), rt
    rawId: The source of the scan in the raw file (ie the source -- mzML, raw, etc.)
    """
    def __init__(self):
        self.scans = []
        self.rt = 0.0
        self.ms_level = None
        self._rawId = None
        self.parent = None

        self.centroid = False

        #SRM/MRM parameters
        self.product_ion = 0

    @property
    def rawId(self):
        if self._rawId is None:
            return getattr(self, 'id', None)
        return self._rawId

    @rawId.setter
    def rawId(self, value):
        self._rawId = value

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
        self.acc = ''# accession information

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
            # modType = ""
            # masses = config.MODIFICATION_MASSES
            # for i in masses:
            #     if tmass in masses[i] or smass == precFormat%masses[i][0]:
            #         #found it
            #         modType = i
            # if not modType:
            #     sys.stderr.write('mod not found %s\n'%modMass)
        self.mods.add((aa,str(position),str(modMass),str(modType)))

    def getModifications(self):
        return '|'.join([','.join(i) for i in self.mods])

    @property
    def modifiedPeptide(self):
        peptide = list(self.peptide)
        for i in self.mods:
            aa, pos = i[0], i[1]
            if aa not in config.RESIDUE_MASSES:
                continue
            pos = int(pos)
            assert peptide[pos].upper() == aa.upper(), 'Amino acid {} in is not equal to modification amino acid {}'.format(peptide[pos], aa)
            peptide[pos] = aa.lower()
        return ''.join(peptide)

    @property
    def modifiedMass(self):
        mass = self.mass
        for i in self.mods:
            mass += float(i[2])
        return mass

    def getTheorMass(self):
        aa_info = [config.RESIDUE_MASSES[i.upper()] for i in self.peptide if i in config.RESIDUE_MASSES]
        peptide_mass = sum([i[0] for i in aa_info])+config.MODIFICATION_MASSES['h'][0]+config.MODIFICATION_MASSES['oh'][0]
        modifications = sum([float(i[2]) for i in self.mods])
        mass = peptide_mass+modifications+float(self.charge)*config.PROTON
        return mass/float(self.charge)