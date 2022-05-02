#!/usr/bin/env python

from bs4 import BeautifulSoup as bs
# import simplejson as json
import requests
from sys import argv, exit
from time import sleep

# taxid = argv[1]
taxid = '1764'
# eutils link for genome ids when querying on taxid
tax_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=genome&term=txid{}&retmode=json"
genome_url = "https://www.ncbi.nlm.nih.gov/genome/?term={}"

def fetch(url, req_type):
    timer = 0
    try:
        r = requests.get(url)
        while r.status_code != 200:
            r = requests.get(url)
            sleep(10)
            timer += 10
            if timer > 100:
                raise Exception
        return r
    except requests.exceptions.RequestException as e:
        print(e)
        exit(1)
    except:
        print('Too many tries for request: {}\n{}'.format(req_type, url))
        exit(1)


# request results for taxid, request object is returned when successful
tax_req = fetch(tax_url.format(taxid), 'genome id to tax url')

if tax_req:
    req_json = tax_req.json()
    genome_ids = req_json['esearchresult']['idlist']
    if len(genome_ids) > 0:
        genome_id = genome_ids[0]
    else:
        print('no genome ids linked to taxid')
        exit(1)
else:
    print('request is empty from taxid')
    exit(1)

genome_req = fetch(genome_url.format(genome_id), 'genbank link from genome page')

contents = genome_req.text
if contents:
    soup = bs(contents, 'html.parser')
    divs = soup.findAll("div", {"class": "refgenome_sensor"})

    success = False
    if len(divs) > 0:
        for i in divs[0].find_all('a'):
            if i.contents[0] == 'GenBank' and i['href'].endswith('gbff.gz'):
                print(i['href'])
                success = True
        if not success:
            print('No reference sequence on Genome page')
    else:
        print('No reference sequence on Genome page')
else:
    print('no contents in Genome xml file')
