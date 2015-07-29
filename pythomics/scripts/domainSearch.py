#!/usr/bin/env python
from __future__ import division, absolute_import

__author__ = 'chris'

description = """
This will search NCBI for domains in a protein fasta file.
"""

import sys
import time
import requests

from pythomics.templates import CustomParser

parser = CustomParser(description=description)
parser.add_fasta()
parser.add_out()
parser.add_argument('--db', help='The database to search', default='cdd', choices=['cdd', 'pfam', 'smart',
                                                                                   'tigrfam', 'cog', 'kog'])

def main():
    args = parser.parse_args()
    files = {'queries': args.fasta}
    nbci_url = 'http://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi?'
    response = requests.post('{nbci}tdata=hits&dmode=std&db={database}&compbasedadj=0&filter=true&evalue=0.0001&cddefl=true'.format(nbci=nbci_url, database=args.db), files=files)
    if response.status_code != 200:
        sys.stderr.write('Error interfacing with NCBI: {}'.format(response.text))
        return 1
    #sys.stderr.write('Submission response is {}\n'.format(response.text))
    header, search_id, _, status = response.text.split('\n')[:4]
    search_id = search_id.split('\t')[1]
    status_code = int(status.split('\t')[1])
    backoff = 1
    while status_code != 0:
        if status_code != 3:
            sys.stderr.write('Error interfacing with NCBI: {}'.format(response.text))
            return 1
        response = requests.get('{nbci}cdsid={job})'.format(nbci=nbci_url, job=search_id))
        header, _, _, status = response.text.split('\n')[:4]
        try:
            status_code = int(status.split('\t')[1])
        except:
            sys.stderr.write('status fail with response\n{}\n'.format(response.text))
            return 1
        if status_code != 0:
            if backoff < 60:
                backoff = int(float(backoff)*1.2+1.0)
            sys.stderr.write('Still running, sleeping for {}\n'.format(backoff))
            time.sleep(backoff)
    with args.out as o:
        o.write(response.text)


if __name__ == "__main__":
    sys.exit(main())
