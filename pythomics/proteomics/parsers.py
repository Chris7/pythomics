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
__author__ = 'Chris Mitchell'

import sys, copy, struct, base64, gzip
from os import path
import pythomics.templates as templates
from pythomics.proteomics.structures import PeptideObject, ScanObject
from pythomics.proteomics import config
import re, os, sqlite3, StringIO, zipfile, time
try:
    from lxml import etree
except ImportError:
    import xml.etree.cElementTree as etree

try:
    import msparser
except ImportError:
    print 'msparser not found, Mascot DAT files unable to be parsed'

#regex for common use
scanSplitter = re.compile(r'[\t\s]')
distillerParse = re.compile(r'_DISTILLER_RAWFILE\[(\d+)\]=\(1\)(.+)')
lastSplit = re.compile(r'.+[/\\](.+)')


class MZMLIterator(templates.GenericIterator):
    def __init__(self, filename):
        """An iterator over the mzML file format.
         I just need a parser right now so the specs might not be exactly correct

        The returned items are ScanObjects
        
        """
        super(MZMLIterator, self).__init__(filename)
        try:
            if isinstance(self.filename, gzip.GzipFile):
                print 'gzip not supported'
                dom1 = etree.parse(self.filename).findall('{http://psi.hupo.org/ms/mzml}spectrum')#, tag=('{http://psi.hupo.org/ms/mzml}spectrum', '{http://psi.hupo.org/ms/mzml}indexedmzML'))
            else:
                dom1 = etree.iterparse(self.filename, tag=('{http://psi.hupo.org/ms/mzml}spectrum', '{http://psi.hupo.org/ms/mzml}indexedmzML'))
            self.lxml = True
        except NameError:
            self.lxml = False
            print 'MZML parsing unavailable: lxml is required to parse mzML files'
            return
        if self.lxml:
            self.spectra = dom1
            self.scans = {}
            self.ra = {}
        else:
            #this isn't implemented
            self.nest = 0
        self.startIter = True
        self.db = None
            
    def __iter__(self):
        return self

    def _get_scan_from_string(self, value):
        return value[value.find('scan=')+5:]

    def parselxml(self, spectra, full=False, namespace='{http://psi.hupo.org/ms/mzml}'):
        if spectra.tag == '{0}indexedmzML'.format(namespace):
            # read our index in
            self.ra = dict([(self._get_scan_from_string(i.values()[0]), i.text) for i in spectra.findall('{0}indexList/{0}index/'.format(namespace))])
        elif spectra.tag == '{0}spectrum'.format(namespace):
            scanObj = ScanObject()
            if full:
                mzmls, intensities = spectra.findall('{0}binaryDataArrayList/'.format(namespace))
                mzmls = base64.b64decode(mzmls.find('{0}binary'.format(namespace)).text)
                intensities = base64.b64decode(intensities.find('{0}binary'.format(namespace)).text)
                mzmls = [struct.unpack('d', mzmls[i:i+8])[0] for i in xrange(0,len(mzmls),8)]
                intensities = [struct.unpack('f', intensities[i:i+4])[0] for i in xrange(0,len(intensities),4)]
                scanObj.scans = [(i,j) for i,j in zip(mzmls, intensities)]
            spectra_info = dict(zip(spectra.keys(),spectra.values()))
            spectra_params = dict([(i.get('name'), i.get('value')) for i in spectra.findall('{0}cvParam'.format(namespace))])
            scan_info = dict([(i.get('name'), i.get('value')) for i in spectra.findall('{0}scanList/{0}scan/{0}cvParam'.format(namespace))])
            precursor_info = dict([(i.get('name'), i.get('value')) for i in spectra.findall('{0}precursorList/{0}precursor/{0}selectedIonList/{0}selectedIon/{0}cvParam'.format(namespace))])
            # print 'spectra info', spectra_info
            # print 'spectra params', spectra_params
            # print 'our scan info', scan_info
            # print 'our precursor info', precursor_info
            ms_level = spectra_params.get('ms level', 0)
            charge = precursor_info.get('charge state', 0)
            precursor_ion = precursor_info.get('selected ion m/z', 0)
            precursor_intensity = precursor_info.get('peak intensity', 0)
            rt = scan_info.get('scan start time', 0)
            rt_length = scan_info.get('ion injection time', 0)
            scanObj.ms_level = int(ms_level)
            scanObj.mass = float(precursor_ion)
            scanObj.charge = charge
            scanObj.title = dict([i.split('=') for i in spectra.get('id').split(' ')]).get('scan', 'No Title')
            scanObj.rt = float(rt)
            return scanObj
        else:
            raise StopIteration
    
    def next(self):
        if self.spectra:
            spectra = self.spectra.next()
        else:
            raise StopIteration
        return self.parselxml(spectra[1])

    def getScan(self, id):
        try:
            if isinstance(self.filename, gzip.GzipFile):
                return self.parselxml(spectra, full=True, namespace='')
            else:
                if not self.ra:
                    [i for i in self]
                self.filename.seek(int(self.ra[str(id)]))
                entry = self.filename.readline()
                row = [entry]
                while '</spectrum>' not in entry:
                    entry = self.filename.readline()
                    row.append(entry)
                spectra = etree.fromstring(''.join(row))
            return self.parselxml(spectra, full=True, namespace='')
        except IndexError:
            print 'wtf'

class PepXMLIterator(templates.GenericIterator):
    def __init__(self, filename):
        """An iterator over the pepXML file format.

        The returned items are ScanObjects

        """
        super(PepXMLIterator, self).__init__(filename)
        try:
            #find our mzml first
            mzml = etree.iterparse(self.filename, tag=('{http://regis-web.systemsbiology.net/pepXML}msms_run_summary',)).next()[1]
            info = dict(mzml.items())
            try:
                self.mzml = MZMLIterator('%s%s'%(info['base_name'], info['raw_data_type']))
            except IOError:
                # maybe we moved it, try this
                xml_path, junk = os.path.split(filename)
                junk, file_name = os.path.split(info['base_name'])
                new_path = os.path.join(xml_path, file_name)
                try:
                    self.mzml = MZMLIterator('%s%s'%(new_path, info['raw_data_type']))
                except IOError:
                    # one more try for a gz
                    self.mzml = MZMLIterator('%s%s.gz'%(new_path, info['raw_data_type']))
            self.filename.seek(0)
            dom1 = etree.iterparse(self.filename, tag=('{http://regis-web.systemsbiology.net/pepXML}spectrum_query',))
            self.lxml = True
        except NameError:
            self.lxml = False
            print 'pepXML parsing unavailable: lxml is required to parse mzML files'
            return
        if self.lxml:
            self.spectra = dom1
            self.scans = {}
            self.ra = {}
        else:
            #this isn't implemented
            self.nest = 0
        self.ptm_regex = re.compile(r'PTMProphet\_[A-Za-z]+([0-9\.\-]+)')
        self.db = None

    def __iter__(self):
        return self

    def _get_scan_from_string(self, value):
        return value[value.find('scan=')+5:]

    def parselxml(self, spectra, full=False, namespace='{http://regis-web.systemsbiology.net/pepXML}'):
        if spectra.tag == '{0}spectrum_query'.format(namespace):
            pepObj = PeptideObject()
            pep_info = dict(spectra.items())
            title = pep_info.get('spectrum', False)
            start_scan = str(int(pep_info.get('start_scan')))
            if not title:
                title = pep_info.get('start_scan', 'Unknown')
            pepObj.title = title
            pepObj.id = title
            pepObj.charge = pep_info.get('assumed_charge', 0)
            pepObj.mass = float(pep_info.get('precursor_neutral_mass', 0))
            rt = float(pep_info.get('retention_time_sec', 0))
            rt = '{0}.{1:02f}'.format(int(rt/60), int(rt%60))
            pepObj.rt = rt
            # get our search result
            search_result = spectra.find('{0}search_result/'.format(namespace))
            search_info = dict(search_result.items())
            rank = search_info.get('hit_rank', 0)
            peptide = search_info.get('peptide', 'X')
            analysis = search_result.findall('{0}analysis_result'.format(namespace))
            pprophet_expect = False
            iprophet_expect = False
            ptmprophet_info = {}
            for i in analysis:
                info = dict(i.items())
                analysis_type = info.get('analysis', False)
                # print analysis_type
                if analysis_type == 'peptideprophet':
                    results = i.find('{0}peptideprophet_result'.format(namespace))
                    prophet_info = dict(results.items())
                    pprophet_expect = prophet_info.get('probability', 0)
                elif analysis_type == 'interprophet':
                    results = i.find('{0}interprophet_result'.format(namespace))
                    inter_info = dict(results.items())
                    iprophet_expect = inter_info.get('probability', 0)
                elif analysis_type == 'ptmprophet':
                    results = i.find('{0}ptmprophet_result'.format(namespace))
                    ptm_info = dict(results.items())
                    probs = results.findall('{0}mod_aminoacid_probability'.format(namespace))
                    # print ptm_info
                    mod_mass = self.ptm_regex.match(ptm_info['ptm']).group(1)
                    for prob in probs:
                        prob_info = dict(prob.items())
                        mod_pos = int(prob_info['position'])-1
                        try:
                            ptmprophet_info[mod_pos][mod_mass] = prob_info['probability']
                        except KeyError:
                            ptmprophet_info[mod_pos] = {mod_mass: prob_info['probability']}
            if iprophet_expect:
                pepObj.expect = iprophet_expect
            elif pprophet_expect:
                pepObj.expect = pprophet_expect
            else:
                pepObj.expect = 1
            accession = [search_info.get('protein', '')]
            alt_proteins = search_result.findall('{0}alternative_protein'.format(namespace))
            accession += [i.get('protein') for i in alt_proteins]
            pepObj.peptide = peptide
            pepObj.acc = ';'.join(accession)
            if ptmprophet_info:
                for position, modification in ptmprophet_info.iteritems():
                    mod_aa = peptide[position]
                    for mass, prob in modification.iteritems():
                        pepObj.addModification(mod_aa, position, mass, '%s(%s)'%(mass, prob))
            else:
                mods = search_result.findall('{0}mod_aminoacid_mass'.format(namespace))
                for i in mods:
                    mod_info = dict(i.items())
                    # print i.items()
                    # print mod_info
                    mod_pos = int(mod_info['position'])-1
                    mod_aa = peptide[mod_pos]
                    pepObj.addModification(mod_aa, mod_pos, mod_info['mass'], mod_info['mass'])
            self.scans[title] = pepObj
            # get our mz/intensities
            scanInfo = self.mzml.getScan(start_scan)
            pepObj.scans = scanInfo.scans
            return pepObj
        else:
            raise StopIteration

    def next(self):
        if self.spectra:
            spectra = self.spectra.next()
        else:
            raise StopIteration
        return self.parselxml(spectra[1])

    def getScan(self, id):
        try:
            return self.scans[id]
        except IndexError:
            print 'wtf'

class XTandemXML(templates.GenericIterator):
    """
    Parser for X!Tandem XML Files.
    """
    def __init__(self, filename, **kwrds):
        super(XTandemXML, self).__init__(filename)
        #parse in our X!Tandem xml file
        if kwrds and kwrds['exclude']:
            exclude = set(kwrds['exclude'])
        else:
            exclude = set()
        try:
            dom1 = etree.parse(self.filename)
            self.lxml = True
        except NameError:
            self.lxml = False
            print 'XTandem parsing unavailable: lxml is required to parse X!tandem xml files due to the namespaces employed'
            return
        if self.lxml:
            self.group = dom1.findall("group")
            self.scans = {}
            self.num = len(self.group)
        else:
            #this isn't implemented
            self.nest = 0
        self.startIter = True
        self.db = None

    def __iter__(self):
        return self

    def parselxml(self, group):
        try:
            expect = group.attrib["expect"]
        except KeyError:
            self.next()
        subnote = list(group.iter("note"))
        for i in subnote:
            if (i.attrib["label"] == "Description"):
                experiment = i.text.strip()
        charge = group.attrib["z"]
        premass = group.attrib["mh"]
        rt = group.attrib["rt"]
        proteins = list(group.iter("protein"))
        fullProtein = ""
        for protein in proteins:
            scanObj = PeptideObject()
            scanObj.charge = charge
            scanObj.mass =premass*int(charge)
            if rt:
                scanObj.rt = float(rt)
            sgroup = group.iter("group")
            for i in sgroup:
                #This is horridly inefficient...
                if "fragment ion mass spectrum" in i.attrib["label"]:
                    ab = i.iter('{http://www.bioml.com/gaml/}Xdata')
                    for j in ab:
                        mzIter = j.iter('{http://www.bioml.com/gaml/}values')
                        for k in mzIter:
                            mz = [mval for mval in k.text.strip().replace('\n',' ').split(' ')]
                    ab = i.iter('{http://www.bioml.com/gaml/}Ydata')
                    for j in ab:
                        mzIter = j.iter('{http://www.bioml.com/gaml/}values')
                        for k in mzIter:
                            inten = [mval for mval in k.text.strip().replace('\n',' ').split(' ')]
                    for j,k in zip(mz,inten):
                        scanObj.scans.append((float(j),float(k)))
            domain = list(protein.iter("domain"))[0]#we only have one domain per protein instance
            note = list(protein.iter("note"))[0]#same here
            mods = list(protein.iter("aa"))#we can have multiple modifications
            if not self.db:
                files = list(protein.iter("file"))
                self.db = files[0].attrib["URL"]
            id = domain.attrib["id"]
            start = domain.attrib["start"]
            end = domain.attrib["end"]
            peptide = domain.attrib["seq"].replace(' ','')#for some reason it can have spaces
            pExpect = domain.attrib["expect"]
            for mod in mods:
                scanObj.addModification(mod.attrib["type"],int(mod.attrib["at"])-1,float(mod.attrib["modified"]), False)
            scanObj.peptide = peptide
            scanObj.expect = float(pExpect)
            scanObj.id = id
            scanObj.title = id
            scanObj.acc = note.text
            self.scans[id] = scanObj
        return scanObj

    def next(self):
        if self.group:
            group = self.group.pop(0)
        else:
            raise StopIteration
        return self.parselxml(group)

    def getScan(self, id):
        try:
            return self.scans[id]
        except IndexError:
            print 'wtf'

    def getProgress(self):
        return len(self.scans)*100/self.num

class MGFIterator(templates.GenericIterator):
    def __init__(self, filename, **kwrds):
        super(MGFIterator, self).__init__(filename, **kwrds)
        #load our index
        indexFile=''.join(path.splitext(self.filename.name)[0])+'.mgfi'
        self.scanSplit = re.compile(r'[\s\t]')
        self.rand = True
        self.titleMap = {}
        self.ra = {}
        self.tparse = re.compile(r'TITLE=(\d+),(\d+): Scan (\d+) \(rt=(.+)\)')
        self.openMGFIndex(indexFile)

    def openMGFIndex(self, path):
        try:
            f = open(path, 'rb')
            for row in f:
                entry = row.strip().split('\t')
                self.ra[entry[0]] = (int(entry[1]),int(entry[2]))
            try:
                self.epos = entry[2]
            except UnboundLocalError:
                self.epos=1
        except IOError:
            print 'building index for:',path
            if os.path.exists(path):
                raise Exception("Index file path: %s found, but appears incomplete"%path)
            f = open(path, 'wb')
            while True:
                try:
                    self.next()
                except StopIteration:
                    break
            for i in self.ra:
                f.write('%s\t%d\t%d\n'%(i,self.ra[i][0],self.ra[i][1]))
            self.f.seek(0)

    def getScan(self, title):
        """
        allows random lookup
        """
        if self.ra.has_key(title):
            self.f.seek(self.ra[title][0],0)
            toRead = self.ra[title][1]-self.ra[title][0]
            info = self.f.read(toRead)
            scan = self.parseScan(info)
        else:
            return None
        return scan

    def parseScan(self, scan):
        """
        All input follows the BEGIN IONS row and ends before END IONS
        """
        setupScan = True
        foundCharge = False
        foundMass = False
        foundTitle = False
        scanObj = ScanObject()
#        print 'stuff in scan',scan
        for row in scan.split('\n'):
            if not row:
                continue
            entry = row.strip().split('=')
            if len(entry) >= 2:
                if entry[0] == 'PEPMASS':
                    scanObj.mass = float(entry[1])
                    foundMass = True
                elif entry[0] == 'CHARGE':
                    scanObj.charge = entry[1]
                    foundCharge = True
                elif entry[0] == 'TITLE':
#                    if self.titleMap:
#                        pos = entry[1].find(',')
#                        title = self.titleMap[int(entry[1][:entry[1].find(',')])]
#                    else:
                    title = '='.join(entry[1:])
                    foundTitle = True
                    scanObj.title = title
                elif entry[0] == 'RTINSECONDS':
                    scanObj.rt = float(entry[1])
            else:
                mz,intensity = self.scanSplit.split(row.strip())
                scanObj.scans.append((float(mz),float(intensity)))
        if foundCharge and foundMass and foundTitle:
            return scanObj
        return None

    def __iter__(self):
        return self

    def next(self):
        row = self.filename.readline()
        if row == '':
            raise StopIteration
        setupScan = False
        scanInfo = ""
        while row:
            if '_DISTILLER' in row:
                if row:
                    m = distillerParse.match(row)
                    if m:
                        self.titleMap[int(m.group(1))+1] = m.group(2)
                pass
            elif 'BEGIN IONS' in row:
                if self.rand:
                    pStart=self.f.tell()
                setupScan=True
#                newScan=True
            elif 'END IONS' in row:
                scan = self.parseScan(scanInfo)
                if scan:
                    if self.rand:
                        self.ra[scan.title] = (pStart,pos)
                    return scan
                return None
            elif setupScan:
                scanInfo+=row
            pos = self.f.tell()
            row = self.f.readline()

    def getProgress(self):
        return self.f.tell()*100/self.epos

class MascotDATIterator(templates.GenericIterator):
    def __init__(self, filename):
        super(MascotDATIterator, self).__init__(filename)
        self.specParse = re.compile(r'Spectrum(\d+)')
        self.afterSlash = re.compile(r'.+[\\/](.+)$')
        self.f = self.filename.name
        resfile = msparser.ms_mascotresfile(self.f)
        self.hit=0
        results = msparser.ms_peptidesummary(resfile, msparser.ms_mascotresults.MSRES_MAXHITS_OVERRIDES_MINPROB,0,0,"",0,5,"",msparser.ms_peptidesummary.MSPEPSUM_USE_HOMOLOGY_THRESH)
        hasFdr, closeFdr, minP = results.getThresholdForFDRAboveHomology(0.01)
        self.resfile = msparser.ms_mascotresfile(self.f)
        self.results = msparser.ms_peptidesummary(self.resfile, msparser.ms_mascotresults.MSRES_SHOW_SUBSETS|msparser.ms_mascotresults.MSRES_CLUSTER_PROTEINS,minP,9999999,"",0,5,"",msparser.ms_peptidesummary.MSPEPSUM_USE_HOMOLOGY_THRESH)
        self.sparams = msparser.ms_searchparams(self.resfile)
        self.quantMethod = False
        if self.sparams.getQUANTITATION():
            self.quantFile = msparser.ms_quant_configfile()
            self.quantFile.setSchemaFileName(''.join(["http://www.matrixscience.com/xmlns/schema/quantitation_2","quantitation_2.xsd"," http://www.matrixscience.com/xmlns/schema/quantitation_1","quantitation_1.xsd"]))
            if self.resfile.getQuantitation(self.quantFile):
                if not self.quantFile.isValid():#true if errors
                    print 'quant file invalid'
                else:
                    self.quantMethod = self.quantFile.getMethodByName(self.sparams.getQUANTITATION())
                    if self.quantMethod:
                        errors = self.quantFile.validateDocument()
                        if errors:
                            print 'quant errors',errors
                        else:
                            self.qstring = self.quantMethod.getProtocol().getType()
                            self.quantMethod = msparser.ms_quant_method(self.qstring)
        self.vMods = msparser.ms_modvector()
        self.uModFile = msparser.ms_umod_configfile()
        if not self.resfile.getUnimod(self.uModFile):
            self.uModFile = msparser.ms_umod_configfile("unimod.xml", "unimod_2.xsd")
        if self.quantMethod:
            self.quantSubs = {}
            for i in xrange(self.quantMethod.getNumberOfModificationGroups()):
                modGroup = self.quantMethod.getModificationGroupByNumber(i)
                for j in modGroup.getNumberOfLocalDefinitions():
                    print 'qunt sub thing',modGroup.getLocalDefinitions(j)
                    print 'get from the unimod xml later'
            for i in xrange(self.quantMethod.getNumberOfComponents()):
                objComponent = self.quantMethod.getComponentByNumber(i)
                for j in xrange(objComponent.getNumberOfModificationGroups()):
                    modGroup = objComponent.getModificationGroupByNumber(j)
                    for k in modGroup.getNumberOfLocalDefinitions():
                        print 'qunt sub thing',modGroup.getLocalDefinitions(k)
                        print 'get from the unimod xml later'
        self.massFile = msparser.ms_masses(self.uModFile)
        self.modFile = msparser.ms_modfile(self.uModFile, msparser.ms_umod_configfile.MODFILE_FLAGS_ALL)
        if not self.modFile.isValid():
            print 'parser error in mod file, add actual errors later'
            return
        self.massType = msparser.MASS_TYPE_MONO
        if self.sparams.getMASS() == "average":
            self.massType = msparser.MASS_TYPE_AVE
        vInd = 1
        modText = self.resfile.getSectionValueStr(msparser.ms_mascotresfile.SEC_MASSES, "delta%d"%vInd)
        while modText:
            modMass, modName = modText.split(',')
            objMod = self.modFile.getModificationByName(modName)
            if objMod:
                self.vMods.appendModification(objMod)
            else:
                print 'Mod',modName,'not found'
            vInd+=1
            modText = self.resfile.getSectionValueStr(msparser.ms_mascotresfile.SEC_MASSES, "delta%d"%vInd)
        self.vecFixed = msparser.ms_modvector();
        vInd = 1
        modText = self.resfile.getSectionValueStr(msparser.ms_mascotresfile.SEC_MASSES, "FixedMod%d"%vInd)
        while modText:
            modMass,modName = modText.split(',')
            objMod = self.modFile.getModificationByName(modName);
            if objMod:
                self.vecFixed.appendModification(objMod);
            else:
              print 'Mod',modName,'not found'
            vInd+=1
            modText = self.resfile.getSectionValueStr(msparser.ms_mascotresfile.SEC_MASSES, "FixedMod%d"%vInd)
        seci = 1
        seckey = self.resfile.enumerateSectionKeys(msparser.ms_mascotresfile.SEC_MASSES,seci)
        self.groupMasses = {}
        self.elementalMasses = {}
        self.vmMass = {}
        secParse = re.compile(r'delta(\d+)')
        while seckey:
            secMatch = secParse.search(seckey)
            if secMatch:
                secMatchInd = int(secMatch.group(1))
                if secMatchInd>9:
                    secMatchInd = chr(secMatchInd+55)
                self.vmMass[secMatchInd] = self.resfile.getSectionValueStr(msparser.ms_mascotresfile.SEC_MASSES, seckey).split(',')
            elif seckey == 'Ignore%d'%seci:
                pass
            elif seckey == 'FixedMod%d'%seci:
                pass
            elif seckey == 'FixedModResidual%d'%seci:
                pass
            elif seckey == 'FixedModNeutralLoss%d'%seci:
                pass
            elif seckey == 'NeutralLoss%d'%seci:
                pass
            elif seckey == 'NeutralLoss%d_master'%seci:
                pass
            elif seckey == 'NeutralLoss%d_slave'%seci:
                pass
            elif seckey == 'NeutralLoss%d_'%seci:
                pass
            elif seckey == 'ReqPepNeutralLoss%d'%seci:
                pass
            elif seckey == 'PepNeutralLoss%d'%seci:
                pass
            else:
                secval = self.resfile.getSectionValueStr(msparser.ms_mascotresfile.SEC_MASSES, seckey)
                if len(seckey) == 1 or '_term' in seckey:
                    self.groupMasses[seckey.lower()] = secval
                else:
                    self.elementalMasses[seckey.lower()] = secval
            seci+=1
            seckey = self.resfile.enumerateSectionKeys(msparser.ms_mascotresfile.SEC_MASSES,seci)
        if not self.elementalMasses['electron']:
            self.elementalMasses['electron'] = 0.000549
        modind=1
        mod = self.resfile.getSectionValueStr(msparser.ms_mascotresfile.SEC_MASSES, "delta%d"%modind)
        while mod:
            print mod
            modind+=1
            mod = self.resfile.getSectionValueStr(msparser.ms_mascotresfile.SEC_MASSES, "delta%d"%modind)
        self.aahelper = msparser.ms_aahelper(self.resfile, "enzymes.txt")
        self.aahelper.setMasses(self.massFile)
        self.aahelper.setAvailableModifications(self.vecFixed,self.vMods)

        self.numHits = self.results.getNumberOfHits()
        self.cRules = set([msparser.ms_fragmentationrules.FRAG_IMMONIUM,msparser.ms_fragmentationrules.FRAG_INTERNAL_YA,msparser.ms_fragmentationrules.FRAG_INTERNAL_YB])
        self.fragRuleMapping = {msparser.ms_fragmentationrules.FRAG_IMMONIUM: 'immonium ion',
                                msparser.ms_fragmentationrules.FRAG_A_SERIES: 'a ion',
                                msparser.ms_fragmentationrules.FRAG_A_MINUS_NH3: 'a ion - NH3',
                                msparser.ms_fragmentationrules.FRAG_A_MINUS_H2O: 'a ion - H2O',
                                msparser.ms_fragmentationrules.FRAG_B_SERIES: 'b ion',
                                msparser.ms_fragmentationrules.FRAG_B_MINUS_NH3: 'b ion - NH3',
                                msparser.ms_fragmentationrules.FRAG_B_MINUS_H2O: 'b ion - H2O',
                                msparser.ms_fragmentationrules.FRAG_C_SERIES: 'c ion',
                                msparser.ms_fragmentationrules.FRAG_X_SERIES: 'x ion',
                                msparser.ms_fragmentationrules.FRAG_Y_SERIES: 'y ion',
                                msparser.ms_fragmentationrules.FRAG_Y_MINUS_NH3: 'y ion - NH3',
                                msparser.ms_fragmentationrules.FRAG_Y_MINUS_H2O: 'y ion - H2O',
                                msparser.ms_fragmentationrules.FRAG_Z_SERIES: 'z ion',
                                msparser.ms_fragmentationrules.FRAG_INTERNAL_YA: 'internal ya ion',
                                msparser.ms_fragmentationrules.FRAG_INTERNAL_YB: 'internal yb ion',
                                msparser.ms_fragmentationrules.FRAG_Z_PLUS_1: 'z+1 ion',
                                msparser.ms_fragmentationrules.FRAG_D_SERIES: 'd ion',
                                msparser.ms_fragmentationrules.FRAG_V_SERIES: 'v ion',
                                msparser.ms_fragmentationrules.FRAG_W_SERIES: 'w ion',
                                msparser.ms_fragmentationrules.FRAG_Z_PLUS_2: 'z+2 ion'
                                }
        self.scanMap = {}

    def __iter__(self):
        return self

    def next(self):
        self.hit+=1
        return self.parseScan(self.hit)

    def getScan(self, title,hit,rank):
        return self.scanMap[int(hit),int(rank)]
#        scan = self.parseScan(int(hit),peprank=int(rank),full=True)
#        return scan

    def parseScan(self, hit,full=True,peprank=False):
        prot = self.results.getHit(hit)
        if prot:
            num_peps = prot.getNumPeptides()
            rms = prot.getRMSDeltas(self.results)
            for i in range(1, 1+ num_peps) :
                query = prot.getPeptideQuery(i)
                p     = prot.getPeptideP(i)
                if p != -1 and query != -1:
                    pep = self.results.getPeptide(query, p)
                    rank = pep.getRank()
                    if peprank and rank != peprank:
                        continue
                    if not pep:
                        continue
                    pstart = prot.getPeptideStart(i)
                    scanObj = PeptideObject()
                    scanObj.hit = hit
                    scanObj.rank = rank
                    peptide = pep.getPeptideStr()
                    scanObj.peptide = peptide
                    charge = pep.getCharge()
                    scanObj.charge = charge
                    mass = self.resfile.getObservedMrValue(query)
                    scanObj.mass = float(mass)
                    vmods = pep.getVarModsStr()
                    for aanum,residue in enumerate(vmods):
                        try:
                            residue = int(residue)
                        except ValueError:
                            #we're a character
                            residue = 9+ord(residue)-64
                        if residue:
                            modName = self.sparams.getVarModsName(residue)
                            modMass = self.sparams.getVarModsDelta(residue)
                            scanObj.addModification(peptide[aanum], pstart-1+aanum, modMass, modName)
                        #self.sparams.getVarModsNeutralLosses(residue),self.sparams.getVarModsPepNeutralLoss(residue)
                    svec = msparser.VectorString()
                    ivec    = msparser.vectori()
                    prot.getSimilarProteins(svec,ivec)
                    simProteins = set([prot.getAccession()])
                    for acc,db in zip(svec,ivec):
                        simProteins.add(acc)
                    proteinGroups = ';'.join(simProteins)
                    scanObj.acc = proteinGroups
                    modString = self.results.getReadableVarMods(query,p)
                    inputQuery = msparser.ms_inputquery(self.resfile,query)
                    stitle = inputQuery.getStringTitle(True)
                    if full:
                        for ionIndex in xrange(1,4):
                            for peakIndex in xrange(inputQuery.getNumberOfPeaks(ionIndex)):
                                scanObj.scans.append((inputQuery.getPeakMass(ionIndex,peakIndex),float(inputQuery.getPeakIntensity(ionIndex,peakIndex))))
                        ionsUsed = pep.getSeriesUsedStr()
                        aaHelperStr = pep.getComponentStr()
                        if not aaHelperStr:
                            aaHelperStr = "1"
                        rules = [int(ruleNum) for ruleNum in self.sparams.getRULES().split(',')];
                        rules.sort()
                        chargeRule = 0
                        matched = {'labels':[],'m/z': [], 'intensity': [], 'error': [], 'series': [], 'start': [], 'end': [], 'losses': [], 'charge': []}
                        for rule in (int(arule) for arule in rules):
#                            print 'looking at rule',rule
#                            print 'crules',self.cRules
                            if rule == 1:
                                chargeRule = 1
                            elif rule == 2 or rule == 3:
                                if (rule == 2 and charge >1) or (rule == 3 and charge > 2):
                                    chargeRule = 2
                            else:
                                for chargeState in xrange(charge+1):
#                                    print 'checking chargestate',chargeState,'with chargerule',chargeRule,'on rule',rule
                                    if (chargeState == 1 and chargeRule == 1) or (chargeState == 2 and chargeRule == 2) and rule not in self.cRules:
                                        frags = msparser.ms_fragmentvector()
                                        errs = msparser.ms_errs()
                                        #print pep, rule, chargeState, 0, pep.getMrCalc(), self.sparams.getMassType(), frags, errs
                                        self.aahelper.setMasses(self.massFile)
                                        self.aahelper.setAvailableModifications(self.vecFixed,self.vMods)
                                        fres = self.aahelper.calcFragmentsEx(pep, rule, chargeState, 0, pep.getMrCalc(), self.sparams.getMassType(), frags, errs)
                                        if not fres:
                                            print 'no fragments made'
                                            print frags
                                            print errs.getLastErrorString()
                                        else:
                                            frags.addExperimentalData(self.resfile, query)
                                            for frag in (frags.getFragmentByNumber(fragIndex) for fragIndex in xrange(frags.getNumberOfFragments())):
                                                if frag.getMatchedIonMass() > 0:
                                                    if frag.isInternal():
                                                        startSite = frag.getStart()
                                                        endSite = frag.getEnd()
                                                    else:
                                                        startSite = frag.getColumn()
                                                        endSite = startSite+len(peptide)
                                                    fmass = frag.getMatchedIonMass()
                                                    fint = frag.getMatchedIonIntensity()
                                                    error = frag.getMatchedIonMass()-frag.getMass()
                                                    matched['labels'].append(frag.getLabel())
                                                    matched['m/z'].append(fmass)
                                                    matched['intensity'].append(fint)
                                                    matched['error'].append(error)
                                                    matched['series'].append(frag.getSeriesName())
                                                    matched['start'].append(startSite)
                                                    matched['end'].append(endSite)
                                                    matched['losses'].append(frag.getNeutralLoss())
                                                    matched['charge'].append(chargeState)
                                            try:
                                                fragname = self.fragRuleMapping[rule]
                                            except KeyError:
                                                fragname = ""
                                            print 'fragment name',fragname
                        scanObj.matched = matched
#                        print matched
                    specId = self.specParse.search(stitle).group(1)
                    scanObj.id = specId
                    self.scanMap[hit,rank] = scanObj
                    return scanObj
        else:
            del self.resfile
            del self.results
            del self.sparams
            raise StopIteration
        return None

class ThermoMSFIterator(templates.GenericIterator):
    def __init__(self, filename, clvl=3,srank=1):
        if isinstance(filename,(str,unicode)):
            self.f = open(filename, 'rb')
        else:
            raise Exception(TypeError,"Unknown Type of filename -- must be a file path")
        self.conn = sqlite3.connect(filename, check_same_thread=False)
        #self.conn.row_factory = sqlite3.Row
        self.cur = self.conn.cursor()
        sql = 'select * from fileinfos'
        self.cur.execute(sql)
        self.fileMap = {}
        self.sFileMap = {}
        self.scans = []
        for i in self.cur.fetchall():
            self.fileMap[str(i[0])]=str(i[1])
            self.sFileMap[i[0]]=lastSplit.search(i[1]).group(1)
        #modification table
        sql = 'select a.AminoAcidModificationID,a.ModificationName, a.DeltaMass from aminoacidmodifications a'
        self.cur.execute(sql)
        self.modTable = {}
        for i in self.cur.fetchall():
            self.modTable[i[0]] = (i[1],i[2])
        #We fetch all modifications here for temporary storage because it is VERY expensive to query peptide by peptide (3 seconds per 100 on my 500 MB test file, with 300,000 scans that's horrid)
        sql = 'select pam.PeptideID, GROUP_CONCAT(pam.AminoAcidModificationID), GROUP_CONCAT(pam.Position) from peptidesaminoacidmodifications pam GROUP BY pam.PeptideID'
        self.cur.execute(sql)
        self.mods = {}
        for i in self.cur.fetchall():
            self.mods[i[0]] = i[1:]
        try:
            sql = 'select COUNT(distinct p.SpectrumID) from peptides p where p.PeptideID IS NOT NULL and p.ConfidenceLevel >= %d and p.SearchEngineRank <= %d'%(clvl,srank)
            self.nrows = self.conn.execute(sql).fetchone()[0]
            #sql = 'select sp.spectrum,p.ConfidenceLevel,p.SearchEngineRank,p.Sequence,p.PeptideID,pp.ProteinID,p.SpectrumID from spectra sp left join peptides p on (p.SpectrumID=sp.UniqueSpectrumID) left join peptidesproteins pp on (p.PeptideID=pp.PeptideID) where p.PeptideID IS NOT NULL'
            sql = 'select GROUP_CONCAT(p.ConfidenceLevel),GROUP_CONCAT(p.SearchEngineRank),GROUP_CONCAT(p.Sequence),GROUP_CONCAT(p.PeptideID), GROUP_CONCAT(pp.ProteinID), p.SpectrumID, sh.Charge, sh.RetentionTime, sh.FirstScan, sh.LastScan, mp.FileID from peptides p join peptidesproteins pp on (p.PeptideID=pp.PeptideID) left join spectrumheaders sh on (sh.SpectrumID=p.SpectrumID) left join masspeaks mp on (sh.MassPeakID=mp.MassPeakID) where p.PeptideID IS NOT NULL and p.ConfidenceLevel >= %d and p.SearchEngineRank <= %d GROUP BY p.SpectrumID'%(clvl,srank)
            self.cur.execute(sql)
        except sqlite3.OperationalError:
            print 'error'
            sql = 'select COUNT(distinct p.SpectrumID) from peptides p where p.PeptideID IS NOT NULL and p.ConfidenceLevel >= %d'%clvl
            self.nrows = self.conn.execute(sql).fetchone()[0]
            #sql = 'select sp.spectrum,p.ConfidenceLevel,p.ConfidenceLevel,p.Sequence,p.PeptideID,pp.ProteinID,p.SpectrumID from spectra sp left join peptides p on (p.SpectrumID=sp.UniqueSpectrumID) left join peptidesproteins pp on (p.PeptideID=pp.PeptideID) where p.PeptideID IS NOT NULL'
            sql = 'select GROUP_CONCAT(p.ConfidenceLevel),GROUP_CONCAT(p.ConfidenceLevel),GROUP_CONCAT(p.Sequence),GROUP_CONCAT(p.PeptideID), GROUP_CONCAT(pp.ProteinID), p.SpectrumID, sh.Charge, sh.RetentionTime, sh.FirstScan, sh.LastScan, mp.FileID from peptides p join peptidesproteins pp on (p.PeptideID=pp.PeptideID) left join spectrumheaders sh on (sh.SpectrumID=p.SpectrumID) left join masspeaks mp on (sh.MassPeakID=mp.MassPeakID) where p.PeptideID IS NOT NULL and p.ConfidenceLevel >= %d GROUP BY p.SpectrumID'%clvl
            self.cur.execute(sql)
        self.index = 0

    def getScan(self, title, specId, peptide):
        """
        get a random scan
        """
        sql = "select sp.Spectrum, p.Sequence, p.PeptideID from spectrumheaders sh left join spectra sp on (sp.UniqueSpectrumID=sh.UniqueSpectrumID) left join peptides p on (sh.SpectrumID=p.SpectrumID) where sh.SpectrumID = %d and p.Sequence = '%s'"%(int(specId),peptide)
        self.cur.execute(sql)
        i = self.cur.fetchone()
        if not i:
            return None
        return self.parseFullScan(i)

    def __iter__(self):
        return self

    def parseScan(self, i):
#sql = 'select GROUP_CONCAT(p.ConfidenceLevel),GROUP_CONCAT(p.SearchEngineRank),GROUP_CONCAT(p.Sequence),GROUP_CONCAT(p.PeptideID), GROUP_CONCAT(pp.ProteinID), p.SpectrumID, sh.Charge, sh.RetentionTime, sh.FirstScan, sh.LastScan, mp.FileID from peptides p join peptidesproteins pp on (p.PeptideID=pp.PeptideID) left join spectrumheaders sh on (sh.SpectrumID=p.SpectrumID) left join masspeaks mp on (sh.MassPeakID=mp.MassPeakID) where p.PeptideID IS NOT NULL and p.ConfidenceLevel = 1 and p.SearchEngineRank = 1 GROUP BY p.SpectrumID'
        objs = []
        self.index+=1
        added = set([])#for some reason redundant scans appear
        for confidence, searchRank, sequence, pepId, proId in zip(i[0].split(','),i[1].split(','),i[2].split(','),i[3].split(','),i[4].split(',')):
            if (sequence,pepId,proId) in added:
                continue
            else:
                added.add((sequence,pepId,proId))
            scanObj = PeptideObject()
            try:
                mods = self.mods[int(pepId)]
                print mods
                for modId, modPosition in zip(mods[0].split(','),mods[1].split(',')):
                    modEntry = self.modTable[int(modId)]
                    scanObj.addModification(sequence[int(modPosition)], modPosition, modEntry[1], modEntry[0])
            except KeyError:
                pass
            scanObj.peptide = sequence
            scanObj.rank = searchRank
            scanObj.confidence = confidence
            scanObj.acc = proId
            scanObj.charge = i[6]
            fName = self.sFileMap[i[10]]
            fScan = i[8]
            lScan = i[9]
            sid = '%s.%s.%s'%(fName, fScan,lScan)
            scanObj.title = sid
            scanObj.id = sid
            scanObj.spectrumId=i[5]
            objs.append(scanObj)
        return objs

    def parseFullScan(self, i):
        """
        parses scan info for giving a Spectrum Obj for plotting. takes significantly longer since it has to unzip/parse xml
        """
        scanObj = PeptideObject()
        sInfo = i[0]
        fp = StringIO.StringIO(sInfo)
        zf = zipfile.ZipFile(fp, 'r')
        peptide = str(i[1])
        pid=i[2]
        sql = 'select aam.ModificationName,pam.Position,aam.DeltaMass from peptidesaminoacidmodifications pam left join aminoacidmodifications aam on (aam.AminoAcidModificationID=pam.AminoAcidModificationID) where pam.PeptideID=%s'%pid
        for row in self.conn.execute(sql):
            scanObj.addModification(peptide[row[1]], str(row[1]), str(row[2]), row[0])
        scanObj.peptide = peptide
        for j in zf.namelist():
            msInfo = zf.read(j)
            msStr = msInfo.split('\n')
            #sIO = StringIO.StringIO('\n'.join(msStr[1:]))
            stage = 0
            #this is dirty, but unfortunately the fastest method at the moment
            print msInfo
            for row in msStr:
                if stage == 0:
                    if 'FileID' in row:
                        finfo = row.split('"')
                        fileName = int(finfo[1])
                        #msScanSum = finfo[5]
                        fName = self.sFileMap[fileName]
                        scanIndex = finfo.index(' ScanNumber=')+1
                        sid = '%s.%s.%s'%(fName, finfo[scanIndex],finfo[scanIndex])
                        #sid = '%s.%s.%s'%(fName, finfo[2],finfo[2])
                        scanObj.title = sid
                        scanObj.id = sid
                        stage=1
                elif stage == 1:
                    if 'PrecursorInfo' in row:
                        finfo = row.split('"')
                        charge = finfo[3]
                        smass = finfo[5]
                        scanObj.charge = charge
                        scanObj.mass = float(smass)
                        stage=2
                elif stage == 2:
                    if '<PeakCentroids>' in row:
                        stage = 3
                elif stage == 3:
                    #we just grab the ms/ms peaks at the moment
                    if 'Peak X' in row:
                        finfo = row.split('"')
                        scanObj.scans.append((float(finfo[1]),float(finfo[3])))
                    elif '</PeakCentroids>' in row:
                        break
        if msInfo:
	    print scanObj.scans
            return scanObj
        else:
            return None

    def next(self):
        if not self.scans:
            i = self.cur.fetchone()
            #we go by groups
            if not i:
                raise StopIteration
            self.scans = self.parseScan(i)
            if not self.scans:
                raise StopIteration
        return self.scans.pop(0)

    def getProgress(self):
        return self.index*100/self.nrows

class indexFolder():
    def __init__(self, folder):
        for root,dirs,files in os.walk(folder):
            for fileName in files:
                if '.mgf' in fileName and '.mgfi' not in fileName:
                    mgfIterator(os.path.join(root,fileName))

class mgfParser(object):
    def __init__(self, filename):
        if type(filename) == type('string'):
            self.f = open(filename)
        else:
            self.f = filename
        self.scans = {}
        newScan=False
        distiller=False
        setupScan = False
        dmap = {}
        dparse = re.compile(r'_DISTILLER_RAWFILE\[(\d+)\]=.+\\(.+)')
        tparse = re.compile(r'TITLE=(\d+),(\d+): Scan (\d+) \(rt=(.+)\)')
        for row in self.f:
            if '_DISTILLER' in row:
                distiller=True
                m = dparse.match(row)
                if m:
                    fname = m.group(2)
                    pos = fname.lower().find('.raw')
                    if pos != -1:
                        fname = fname[:pos]
                    dmap[m.group(1)] = fname
                else:
                    continue
            elif 'BEGIN IONS' in row:
                scanObj = scanObject()
                setupScan=True
                newScan=True
            elif 'END IONS' in row:
                self.scans[title] = scanObj
                newScan = False
            elif setupScan:
                entry = row.strip().split('=')
                if len(entry) >= 2:
                    if entry[0] == 'PEPMASS':
                        scanObj.mass = float(entry[1])
                    elif entry[0] == 'CHARGE':
                        scanObj.charge = entry[1]
                    elif entry[0] == 'TITLE':
                        if distiller:
                            m = tparse.match(row)
#                            print '%s.%s.%s'%(dmap[int(m.group(1))-1],m.group(3),m.group(3))
#                            print m, m.groups()
                        else:
                            title = entry[1]
                            scanObj.title = entry[1]
                else:
                    row.strip()
                    mz,inten = row.strip().split('\t')
                    scanObj.scans.append((float(mz),float(inten)))
                    setupScan=False
            elif newScan and not setupScan:
                mz,inten = row.strip().split('\t')
                scanObj.scans.append((float(mz),float(inten)))
        print 'parsed',len(self.scans)

    def hasScan(self, scan):
        return self.scans.has_key(scan)

    def getScans(self):
        return self.scans

    def getScan(self, scan):
        try:
            return self.scans[scan]
        except KeyError:
            return None