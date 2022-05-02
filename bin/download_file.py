#!/usr/bin/env python

import requests
import wget

url = 'http://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/240/505/GCF_000240505.1_ASM24050v2/GCF_000240505.1_ASM24050v2_genomic.gbff.gz'
r = requests.get(url, allow_redirects=True)

with open('GCF_000240505.1_ASM24050v2_genomic.gbff.gz', 'wb') as out:
    out.write(r.content)
#
#
# import requests
#
# print('Beginning file download with requests')
#
# url = 'http://i3.ytimg.com/vi/J---aiyznGQ/mqdefault.jpg'
# r = requests.get(url)
#
# with open('/Users/scott/Downloads/cat3.jpg', 'wb') as f:
#     f.write(r.content)
#
# # Retrieve HTTP meta-data
# print(r.status_code)
# print(r.headers['content-type'])
# print(r.encoding)

#wget.download(url, 'GCF_000240505.1_ASM24050v2/GCF_000240505.1_ASM24050v2_genomic.gbff.gz')
