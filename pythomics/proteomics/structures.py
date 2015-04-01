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
import pandas as pd
import traceback

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
        self.rawId = None

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


class WideTable(object):
    """
    Breaks up a massive table into individual chunks and stores these chunks as separate tables.
    Extends the select and append method to be pandas-like and handle this indexing to mimic a single
    DataFrame
    """
    def __init__(self, store, chunk_size=None, min_itemsize=255):
        self.store = store
        self.index_name = '_store_index'
        self.chunk_size_name = '_chunk_size'
        self.index_itemsize = min_itemsize
        try:
            self.store_index = self.store.select(self.index_name)
        except KeyError:
            self.store_index = pd.DataFrame(columns=['chunk'])
            self.store.put(self.index_name, self.store_index)
        finally:
            self.tables = [i[0] for i in self.store.iteritems()]
        try:
            # stored as a series
            self.chunk_size = self.store.select(self.chunk_size_name)
        except KeyError:
            self.chunk_size = pd.Series(chunk_size if chunk_size is not None else 100)
            self.store.append(self.chunk_size_name, self.chunk_size)
        self.chunk_size = self.chunk_size[0]

    def get_chunk_name(self, name, chunk):
        return '_{0}_{1}'.format(name, chunk)

    def append(self, name, data, destroy=True):
        if isinstance(data, pd.Series):
            self.store.append(name, data)
        else:
            # get the existing chunks from the data and insert them to those blocks
            chunks = self.store_index['chunk']
            # prefix column name with name of store for multi file storage
            existing_map = chunks[data.columns]
            new = existing_map[existing_map.isnull()] # Series of new columns
            existing = existing_map.dropna().to_dict() # dict of {'column': 'chunk'}, we want {'chunk': columns set}
            existing_rev = {}
            for key, value in existing.iteritems():
                try:
                    existing_rev[value].append(key)
                except KeyError:
                    existing_rev[value] = [key]

            # update existing entries
            columns_added = []
            for chunk_index, chunk_columns in existing_rev.iteritems():
                chunk_name = self.get_chunk_name(name, chunk_index)
                cd = self.store.select(chunk_name)
                cd = cd.merge(data.loc[:, chunk_columns], how='outer')
                self.store.append(chunk_name, cd)
                del cd
                columns_added.extend(chunk_columns)

            new_cols = list(set(data.columns)-set(columns_added))
            if new_cols:
                if destroy:
                    data = data[new_cols]
                else:
                    data = data.copy()[new_cols]
                current_chunk = chunks.astype('int').max()
                if pd.isnull(current_chunk):
                    current_chunk = 0
                space_taken = len(chunks.where(chunks == str(current_chunk)).dropna())
                end = self.chunk_size-space_taken
                if end == 0:
                    # we have no more room
                    current_chunk += 1
                    if len(new_cols) > self.chunk_size:
                        end = self.chunk_size
                    else:
                        end = len(new_cols)
                elif end > len(new_cols):
                    end = len(new_cols)
                start = 0
                while start < len(new_cols):
                    chunk_name = self.get_chunk_name(name, current_chunk)
                    td = data.iloc[:,start:end]
                    self.store.append(chunk_name, td)
                    for col_name in td.columns:
                        self.store_index.loc[col_name,'chunk'] = str(current_chunk)
                    del td
                    self.tables.append(chunk_name)
                    start = end
                    end += self.chunk_size
                    if end > len(new_cols):
                        end = len(new_cols)
                self.store.put(self.index_name, self.store_index)
                del data

    def get_matching_tables(self, name):
        return [i for i in self.tables if i.startswith('/_{0}_'.format(name))]

    def getColumnChunk(self, name, column):
        chunk = self.store_index.loc[str(column), 'chunk']
        return self.get_chunk_name(name, chunk)

    def selectColumn(self, name, column=None):
        chunk_name = self.getColumnChunk(name, column)
        return self.store.select(chunk_name, where="columns=='%s'" % column).dropna()

    def select(self, name, query_columns=None, sort=True, sparse=True, where=None, **kwargs):
        t = pd.DataFrame()
        queries = {}
        if query_columns is not None:
            if not isinstance(query_columns, (list, tuple)):
                query_columns = [query_columns]
            for i in query_columns:
                chunk_name = self.getColumnChunk(name, i)
                col = str(i)
                try:
                    queries[chunk_name].append(col)
                except KeyError:
                    queries[chunk_name] = [col]
        else:
            for i in self.get_matching_tables(name):
                queries[i] = None # None is a key for all
        for chunk_name, query_columns in queries.iteritems():
            query_kw = {}
            if where is not None:
                query_kw.update({'where': where})
            if query_columns is not None:
                query_kw.update({'columns': query_columns})
            values = self.store.select(chunk_name, **query_kw).dropna(axis=[0, 1], how='all')
            if values.empty:
                continue
            t = pd.concat([t,values])
        if sort:
            t.sort_index(inplace=True)
        if sparse:
            return t.to_sparse()
        else:
            return t

    def select_mz(self, name, start_mz=None, end_mz=None, columns=None):
        if isinstance(start_mz, (list,tuple)):
            query = ' | '.join(['(index >= {0} & index <= {1})'.format(start, end) for start, end in zip(start_mz, end_mz)])
        else:
            query = 'index > {0} & index < {1}'.format(start_mz, end_mz)
        return self.select(name, where=query, query_columns=columns)

    def create_table_index(self, name, *args, **kwargs):
        for i in self.get_matching_tables(name):
            try:
                self.store.create_table_index(i, *args, **kwargs)
            except:
                sys.stderr.write('Traceback on create table index: {0}'.format(traceback.format_exc()))
                continue