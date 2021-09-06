import pandas as pd
import json
import requests
from bs4 import BeautifulSoup, UnicodeDammit
import gzip
import argparse
import logging
import logging.config
import sys


'''
This script retrieves TEP (target enabling package) from the structural genomics consoritum.

Output example:

ENSG00000118007: {
                id: 'ENSG00000118007',
                symbol: 'STAG1',
                link: 'https://www.thesgc.org/tep/stag1'
}

'''
def id_lookup(ensembl_ids):

    headers={ "Content-Type" : "application/json", "Accept" : "application/json"}
    r = requests.post('http://rest.ensembl.org/lookup/id', headers=headers, data=json.dumps({ "ids" : ensembl_ids }))

    decoded = r.json()

    # Parse response:
    parsed = []
    for gene_id, data in decoded.items():
        parsed.append({
            'gene_id': gene_id,
            'symbol': data['display_name']
            })

    return pd.DataFrame(parsed)


def uniprot_lookup(uniprot_id):
    url = f'http://rest.ensembl.org/xrefs/symbol/homo_sapiens/{uniprot_id}?content-type=application/json'
    r = requests.get(url)
    data = r.json()

    # Parse gene id:
    for item in data:
        if item['type']:
            return item['id']

    # If gene id is not found:
    logging.info(f'Failed to retrieve Ensembl id for: {uniprot_id}')
    return None


def retrieve_tep_list():

    def get_gene_name(row):
        return row.findAll('td')[1].text

    def get_url(row):
        gene_cell = row.findAll('td')[0]
        tep_url = gene_cell.find('a').get('href')
        
        if not tep_url.startswith('http'):
            tep_url = 'https://www.thesgc.org' + tep_url
            
        return tep_url
        
    def get_therapeutic_area(row):
        therapeutic_area = row.findAll('td')[2].text
        return therapeutic_area


    url = 'https://www.thesgc.org/tep'

    response = requests.get(url)

    html = response.text
    uhtml = UnicodeDammit(html)

    soup = BeautifulSoup(uhtml.unicode_markup, features="html.parser")

    TEP_table = soup.findAll('table')[-1]

    TEP_raw_data = []

    for row in TEP_table.find('tbody').findAll('tr'):
        TEP_raw_data.append({
            "TEP_url": get_url(row),
            "disease": get_therapeutic_area(row),
            "description": get_gene_name(row)
        })
        
    return pd.DataFrame(TEP_raw_data)


def tep_lookup(link):
    r = requests.get(link)
    html = r.text
    uhtml = UnicodeDammit(html)

    soup = BeautifulSoup(uhtml.unicode_markup, features="html.parser")
    uniprot_ids = []

    for a in soup.findAll('a'):
        url = a.get('href')
        try:
            if 'uniprot' in url:
                uniprot_ids.append(url.split('/')[-1])
        except TypeError:
            continue

    if uniprot_ids == []:
        logging.info(f'Failed to retrieve uniprot ids from this TEP: {link}')

    return uniprot_ids

def write_output(tep_df, outputFile, jsonl=True):
    if jsonl:
        tep_df.to_json(outputFile, orient='records', lines=True, compression='gzip')
    else:
        TEPs = {row['gene_id']:{
                'id': row['gene_id'],
                'symbol': row['symbol'],
                'link': row['TEP_url'],
                'disease': row['disease'],
                'uniprot_id': row['uniprot_id']}
                for i, row in tep_df.iterrows()}
        with gzip.open(outputFile, 'wt', encoding="ascii") as f:
            json.dump(TEPs, f)

def main():

    # Reading output file name from the command line:
    parser = argparse.ArgumentParser(description="This script fetches TEP data from the Structural Genomics Consortium.")
    parser.add_argument('--output', '-o', type=str, help='Output file. gzipped JSON', required = True)
    parser.add_argument('--lines', required=False, action='store_false', help='Flag to output the results in JSON multiline.')  
    parser.add_argument('--logFile', type=str, help='File into which the logs are saved', required=False)
    args = parser.parse_args()

    outputFile = args.output

    # If no logfile is specified, logs are written to the standard error:
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)s %(module)s - %(funcName)s: %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
    )
    if args.logFile:
        logging.config.fileConfig(filename=args.logFile)
    else:
        logging.StreamHandler(sys.stderr)

    # The TEP list compiled as dataframe:
    tep_list = retrieve_tep_list()

    logging.info(f'Number of TEPs retrieved: {len(tep_list)}')

    # Adding uniprot IDs:
    logging.info(f'Retrieving uniprot IDs.')
    tep_list['uniprot_id'] = tep_list.TEP_url.apply(tep_lookup)

    # Explode uniprot IDs:
    tep_list = tep_list.explode('uniprot_id')
    logging.info(f'After exploding data by uniprot IDs, number of rows: {len(tep_list)}')

    # Fetching corresponding ensembl id:
    logging.info(f'Fetching Ensembl IDs for uniprot ids.')
    tep_list['gene_id'] = tep_list.uniprot_id.apply(uniprot_lookup)

    # Get the symbols for the genes:
    logging.info(f'Retrieving gene symbols for each Ensembl ID.')
    gene_ids = tep_list.gene_id.unique().tolist()

    symbol_mapping_df = id_lookup(gene_ids)

    tep_merged = tep_list.merge(symbol_mapping_df, on='gene_id', how='outer')

    # Printing:
    logging.info(f'Saving data to {outputFile}.')
    
    write_output(tep_merged, outputFile)

if __name__ == '__main__':
    main()
