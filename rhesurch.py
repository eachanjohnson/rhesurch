#!/usr/bin/env python

#
# Rhesearch
#
# (c) 2016 Eachan Johnson
#

# Import things we need
import __main__
__main__.pymol_argv = ['pymol','-qc']  # Pymol: quiet and no GUI
import csv  # For reading and writing CSV files
#import rdkit  # For cheminformatics
import urllib, urllib2, requests  # For communicating with search servers
import json  # For parsing result data
import time  # For checking for results at intervals
import xml.etree.ElementTree as ET  # For parsing XML results
import pymol
pymol.finish_launching()
import sys, os  # For moving and renaming files
import shutil
import docopt


__author__ = 'Eachan Johnson'

__doc__ = '''
Usage: rhesearch
'''

# Define classes
## Primitive classes
class Searcher(object):

    def __init__(self, type):
        self.type = type

class Converter(object):

    def __init__(self, base_url):
        self.base_url = base_url


class LocalConverter(Converter):

    def __init__(self, filename):
        """
        :param filename: CSV file providing keys in the first column, and values in the second column
        :return: Abstraction of a dictionary built from a CSV file.
        """
        Converter.__init__(self, base_url=filename)
        self.filename = filename
        self.dict = self.generate_dict()

    def generate_dict(self):
        d = {}
        with open(self.filename, 'rU') as f:
            c = csv.reader(f)
            for row in c:
                key = row[0]
                value = row[1]
                d[key] = value
        return d

    def convert(self, key):
        """
        :param key: Item to be converted
        :return: Product of conversion
        """
        if type(key) != list:
            return self.dict[key]
        else:
            l = []
            for number in key:
                try:
                    l.append(self.dict[number])
                except KeyError:
                    l.append('C')
            return l


class XMLConverter(Converter):

    def __init__(self, base_url):
        Converter.__init__(self, base_url=base_url)

    def get_xml(self, parameters, contact):
        #print 'Getting XML'
        data = urllib.urlencode(parameters)
        request = urllib2.Request(self.base_url, data)
        request.add_header('User-Agent', 'Python {}'.format(contact))
        try:
            response = urllib2.urlopen(request)
        except urllib2.HTTPError:
            time.sleep(10)
            try:
                response = urllib2.urlopen(request)
            except urllib2.HTTPError:
                return None
        response_string = response.read()
        return response_string

    def get_data_from_xml(self, xml_string, data_to_retrieve):
        #print 'Getting data from {}'.format(xml_string[0:500])
        d = {}
        try:
            response_tree = ET.fromstring(xml_string)
        except ET.ParseError:
            try:
                response_tree = ET.fromstring(xml_string)
            except ET.ParseError:
                return None
        if type(data_to_retrieve) != list:
            data_to_retrieve = [data_to_retrieve]
        else:
            pass
        for key in data_to_retrieve:
            if type(key) != list:
                for n, entry in enumerate(response_tree.iter(key)):
                    if n == 0:
                        d[key] = entry.text
                    else:
                        break
            else:
                for n, entry1 in enumerate(response_tree.iter(key[0])):
                    if n == 0:
                        try:
                            d[key[1]] = entry1.find(key[1]).text
                        except:
                            d[key[1]] = ''
                        #print d
                    else:
                        break
        return d


class JsonResult(object):

    def __init__(self, json_data):
        """
        :param json_data:
        :return:
        """
        self.json_data = json_data
        self.__dict__ = json.load(self.json_data)
        self.dict = self.__dict__


class PymolSession(object):

    def __init__(self, pdb):
        pymol.cmd.reinitialize()
        if type(pdb) != list:
            self.pdb = [pdb]
        else:
            self.pdb = pdb
        self.filename = '{}.pse'.format('-'.join([code for code in self.pdb]))
        for pdb in self.pdb:
            pymol.cmd.fetch(pdb)

    def save(self):
        print('Saving PyMOL session: {}'.format(self.filename))
        pymol.cmd.save(self.filename)
        return self.filename

    def cealign(self, mol_list):
        print('Getting RMSD of {} and {} using PyMOL.'.format(mol_list[0], mol_list[1]))
        try:
            rmsd = pymol.cmd.cealign('{} and chain A'.format(mol_list[0]), '{} and chain A'.format(mol_list[1]))['RMSD']
        except:
            try:
                rmsd = pymol.cmd.cealign(mol_list[0], mol_list[1])['RMSD']
            except:
                rmsd = 'NA'
        pymol.cmd.orient()
        return rmsd

    def cartoon(self):
        pymol.cmd.hide('everything')
        pymol.cmd.show('cartoon')
        return None

    def ray(self, w=2000, h=2000, space='cmyk', bg='white'):
        pymol.cmd.space(space)
        pymol.cmd.bg_color(bg)
        pymol.cmd.ray(str(w), str(h))
        filename = '{}.png'.format(self.filename.split('.pse')[0])
        pymol.cmd.save(filename)
        return filename

    def reinit(self):
        pymol.cmd.reinitialize()
        return None

    def quit(self):
        pymol.cmd.quit()
        return None

## Convenience classes
class NcbiGiToUniprotConverter(XMLConverter):

    def __init__(self):
        """
        :return: Abstraction of a search engine which takes NCBI GI numbers and returns UniProt accession codes and \
        names.
        """
        XMLConverter.__init__(self, base_url='http://www.uniprot.org/mapping/')

    def convert(self, gi):
        d = {'accession': '', 'fullName': ''}
        #print 'Getting UniProt ID for gi|{}'.format(gi)
        params = {
            'from':'P_GI',
            'to': 'ACC',
            'format': 'xml',
            'query': '{}'.format(gi)
        }
        #print 'Constructed params'
        response_string = self.get_xml(parameters=params, contact='ejohnson@broadinstitute.org')
        #print 'Got response'
        if response_string is not None:
            #print 'Cleaning Uniprot XML'
            head = response_string.split('<uniprot')[0]
            tail = '">'.join(response_string.split('">')[1:])
            response_string = head + '<uniprot>' + tail
            d = self.get_data_from_xml(xml_string=response_string,
                                       data_to_retrieve=['accession', ['recommendedName', 'fullName']])
        else:
            pass
        if d is not None:
            return d
        else:
            return {'accession': '', 'fullName': ''}


class RvToUniprotConverter(XMLConverter):

    def __init__(self):
        """
        :return: Abstraction of a search engine which takes NCBI GI numbers and returns UniProt accession codes and \
        names.
        """
        XMLConverter.__init__(self, base_url='http://www.uniprot.org/mapping/')

    def convert(self, rv_number):
        d = {'accession': '', 'fullName': ''}
        #print 'Getting UniProt ID for gi|{}'.format(gi)
        params = {
            'from':'TUBERCULIST_ID',
            'to': 'ACC',
            'format': 'xml',
            'query': '{}'.format(rv_number)
        }
        #print 'Constructed params'
        response_string = self.get_xml(parameters=params, contact='ejohnson@broadinstitute.org')
        #print 'Got response'
        if response_string is not None:
            #print 'Cleaning Uniprot XML'
            head = response_string.split('<uniprot')[0]
            tail = '">'.join(response_string.split('">')[1:])
            response_string = head + '<uniprot>' + tail
            d = self.get_data_from_xml(xml_string=response_string,
                                       data_to_retrieve=['accession', ['recommendedName', 'fullName']])
        else:
            pass
        if d is not None:
            print d['fullName']
            return d
        else:
            return {'accession': '', 'fullName': ''}


class PubChemSearch(Searcher):

    def __init__(self):
        """
        :return: Abstraction of a search engine for PubChem.
        """
        Searcher.__init__(self, type = 'PubChem')
        self.base_url = 'http://pubchem.ncbi.nlm.nih.gov/rest/pug'

    def cleanup_smiles(self, smiles):
        clean_smiles = '%5B'.join(smiles.split('['))
        clean_smiles = '%5D'.join(clean_smiles.split(']'))
        clean_smiles = '%23'.join(clean_smiles.split('#'))
        clean_smiles = '%2F'.join(clean_smiles.split('/'))
        clean_smiles = '%3D'.join(clean_smiles.split('='))
        clean_smiles = '%2B'.join(clean_smiles.split('+'))
        clean_smiles = '%5C'.join(clean_smiles.split('\\'))
        #print 'Cleaned smiles is {}'.format(clean_smiles)
        return clean_smiles

    def bioassay_from_smiles(self, smiles):
        result = ''
        print('Getting PubChem BioAssays for compound with SMILES {}'.format(smiles))
        clean_smiles = self.cleanup_smiles(smiles)
        search_url = '{}/compound/smiles/{}/assaysummary/CSV'.format(self.base_url, clean_smiles)
        try:
            search_csv = urllib2.urlopen(url=search_url)
        except urllib2.HTTPError:
            time.sleep(10)
            try:
                search_csv = urllib2.urlopen(url=search_url)
            except urllib2.HTTPError:
                try:
                   time.sleep(10)
                   search_csv = urllib2.urlopen(url=search_url)
                except urllib2.HTTPError:
                    search_csv = None
        if search_csv is not None:
            result = search_csv.read()
        return result

    def name_from_smiles(self, smiles):
        result = ''
        #print 'Getting name for compound with SMILES {}:'.format(smiles)
        clean_smiles = self.cleanup_smiles(smiles)
        search_url = '{}/compound/smiles/{}/synonyms/TXT'.format(self.base_url, clean_smiles)
        try:
            search_result = urllib2.urlopen(url=search_url)
        except urllib2.HTTPError:
            time.sleep(1)
            try:
                search_result = urllib2.urlopen(url=search_url)
            except urllib2.HTTPError:
                try:
                   time.sleep(1)
                   search_result = urllib2.urlopen(url=search_url)
                except urllib2.HTTPError:
                    search_result = None
        if search_result is not None:
            result = search_result.read().split('\n')[0]
            #print result
        return result


class BlastSearch(Searcher):

    def __init__(self, organism='txid83332%20%5BORGN%5D'):
        Searcher.__init__(self, type = 'BLAST')
        self.organism = organism
        self.base_url = 'http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?'

    def create_request(self, gi_list):
        print('Searching {} for gi|{}...'.format(self.type, ', gi|'.join(gi_list)))
        search_parameters = {
            'program': 'blastp',
            'service': 'plain',
            'database': 'nr',
            'entrez_query': self.organism,
            'query': '%20AND%20'.join(gi_list),
            'alignment_view': 'pairwise',
            'format_object': 'Alignment',
            'format_type': 'XML',
            'expect': '1e-4',
            'email': 'ejohnson@broadinstitute.org'
        }
        parameter_url = '&'.join(['='.join([param.upper(), search_parameters[param]]) for param in search_parameters])
        search_url = '{}CMD=Put&{}'.format(self.base_url, parameter_url)
        print('Retrieving URL: {}'.format(search_url))
        search_html = urllib2.urlopen(url=search_url)
        qblast_info = search_html.read().split('QBlastInfoBegin')[1].split('QBlastInfoEnd')[0]
        rid = qblast_info.split('RID = ')[1].split('RTOE')[0].rstrip()
        return rid

    def retrieve_search(self, rid):
        print('RID is {}'.format(rid))
        rid_url = '{}CMD=Get&RID={}&FORMAT_TYPE=XML'.format(self.base_url, rid)
        print('Retrieving URL: {}'.format(rid_url))
        rid_xml = urllib2.urlopen(url=rid_url).read()
        try:
            status = rid_xml.split('QBlastInfoBegin')[1].split('QBlastInfoEnd')[0].split('Status=')[1].rstrip()
        except IndexError:
            status = 'READY'
        else:
            print(status)
        counter = 0
        if status == 'WAITING':
            while status == 'WAITING':
                #print counter
                if counter >= 10 and status != 'READY':
                    return None
                counter += 1
                time.sleep(60)
                rid_xml = urllib2.urlopen(url=rid_url).read()
                try:
                    status = rid_xml.split('QBlastInfoBegin')[1].split('QBlastInfoEnd')[0].split('Status=')[1].rstrip()
                except IndexError:
                    status = 'READY'
                    rid_xml = ET.fromstring(rid_xml)
                print(status)
        elif status == 'READY':
            print(status)
            rid_xml = ET.fromstring(rid_xml)
        else:
            self.retrieve_search(rid=rid)
        blast_dict = {}
        for n, hit in enumerate(rid_xml.iter('Hit')):
            name = ''
            for id in hit.iter('Hit_id'):
                hit_id = id.text
            for evalue in hit.iter('Hsp_evalue'):
                signif = float(evalue.text)
            for description in hit.iter('Hit_def'):
                name = description.text.split('>')[0]
            if 'pdb' in hit_id:
                pdb = hit_id.split('pdb|')[1].split('|')[0]
            else:
                pdb = 'unknown_{}'.format(n)
            if signif <= 1e-4:
                blast_dict[hit_id.split('gi|')[1].split('|')[0]] = {'E value': signif, 'name': name, 'pdb': pdb,
                                                                    'gi': hit_id.split('gi|')[1].split('|')[0]}
        time.sleep(3)
        return blast_dict

    def do_search(self, gi_list, rid_list=[]):
        d = {}
        try:
            this_rid = self.create_request(gi_list=gi_list)
        except Exception as e:
            print('Encountered error: {}'.format(e))
            this_rid = self.create_request(gi_list=gi_list)
        d['rid'] = this_rid
        d['pdbs'] = self.retrieve_search(rid=d['rid'])
        if d['pdbs'] is None:
            print('Trying BLAST again...')
            if len(rid_list) != 0:
                for rid in rid_list + [this_rid]:
                    print('Checking back on RID {}'.format(rid))
                    if d['pdbs'] is None:
                        d['rid'] = rid
                        d['pdbs'] = self.retrieve_search(rid=d['rid'])
                    else:
                        pass
                if d['pdbs'] is None:
                    new_rid_list = rid_list
                    new_rid_list.append(this_rid)
                    d = self.do_search(gi_list=gi_list, rid_list=new_rid_list)
                else:
                    pass
            else:
                d = self.do_search(gi_list=gi_list, rid_list=[this_rid])
        else:
            pass
        return d


class PdbSearch(Searcher):

    def __init__(self):
        Searcher.__init__(self, type = 'PDB')
        self.base_url = 'http://www.rcsb.org/pdb/rest/search'

    def name_from_pdb(self, pdb_code):
        pass

    def pdb_from_uniprot(self, uniprot):
        #print 'Getting PDB code for UniProt ID {}'.format(uniprot)
        params = '''
        <orgPdbQuery>
        <queryType>org.pdb.query.simple.UpAccessionIdQuery</queryType>
        <description>{}</description>
        <accessionIdList>{}</accessionIdList>
        </orgPdbQuery>
        '''.format(uniprot, uniprot)
        request = urllib2.Request(self.base_url, data=params)
        response = urllib2.urlopen(request)
        response_string = response.read()
        if response_string:
            pdb = response_string[0:4]
            if pdb != 'Prob':
                return pdb
            else:
                return ''
        else:
            return ''

class StringSearch(Searcher):

    def __init__(self):
        """
        :return: Abstraction of a search engine for the STRING database.
        """
        Searcher.__init__(self, type = 'STRING')
        self.base_url = 'http://string-db.org/api/tsv-no-header/'

    def get_interactors(self, gi):
        search_url = '{}/interactorsList?identifiers={}&species=83332'.format(self.base_url, gi)
        response = urllib2.urlopen(url=search_url)
        interactors = set([line.split('.')[1] for line in response])
        return interactors


class Protein(object):

    def __init__(self):
        self.gi = ''
        self.pdb = ''
        self.name = ''
        self.uniprot = ''
        self.locus = ''
        self.similar = []
        self.similar_rmsd = {}

    def get_name_from_gi(self):
        converter = NcbiGiToUniprotConverter()
        acc_name = converter.convert(gi=self.gi)
        self.uniprot = acc_name['accession']
        try:
            self.name = acc_name['fullName']
        except:
            self.name = ''
        return self.name

    def get_pdb_from_uniprot(self):
        pdb_search = PdbSearch()
        self.pdb = pdb_search.pdb_from_uniprot(uniprot=self.uniprot)
        return self.pdb

    def reinitialize(self):
        if gi != '':
            self.get_name_from_gi()
            if uniprot != '':
                self.get_pdb_from_uniprot()
            else:
                pass
        else:
            pass
        return None

    def update_similar(self, blast_dict):
        for result in blast_dict:
            protein_result = Protein()
            protein_result.gi = result
            protein_result.reinitialize()
            self.similar_pdb.append(protein_result)
        return None

    def update_pdb_rmsd(self, rmsd_dict):
        for rmsd in rmsd_dict:
            self.similar_rmsd[rmsd] = rmsd_dict[rmsd]
        return None


class Compound(object):

    def __init__(self, id, conc, screen_target):
        self.id = id
        smiles_converter = LocalConverter(filename='brd2smiles.csv')
        self.smiles = smiles_converter.convert([self.id])[0]
        pubchem = PubChemSearch()
        self.name = pubchem.name_from_smiles(self.smiles)
        self.concentration = conc
        self.screen_target = screen_target
        self.assay_targets = []

    def update_assay_targets(self, gi_list):
        self.assay_targets += [Protein(gi=gi) for gi in gi_list]
        return None


# Define functions
def main():

    return None

# Boilerplate
if __name__ == '__main__':
    try:
        import __main__
        __main__.pymol_argv = ['pymol','-qc']  # Pymol: quiet and no GUI
        import pymol
        main()
    except KeyboardInterrupt:
        print('Goodbye!')
        quit()
else:
    pass
