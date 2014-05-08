import os, argparse, sys, re
import random
import unittest
import process_seqs
from Bio import Entrez
import urllib2
import gzip
from subprocess32 import check_output, CalledProcessError
import shutil

class TestDownloadFastaFunctions(unittest.TestCase):
    
    def setUp(self):
        Entrez.email = 'dar326@cornell.edu'
        
    def test_get_taxonomy(self):
        process_seqs.get_taxonomy(os.getcwd(), '169')
        self.assertTrue(os.path.isfile('taxonomy.txt'))
        d = {}
        with open("taxonomy.txt") as f:
            for line in f:
                (key, val) = line.split()
                d[key] = val
        self.assertEqual(len(d), 12)
        self.assertEqual(d['superkingdom'], 'Eukaryota')
        self.assertEqual(d['kingdom'], 'Metazoa')
        self.assertEqual(d['phylum'], 'Chordata')
        self.assertEqual(d['subphylum'], 'Craniata')
        self.assertEqual(d['class'], 'Mammalia')
        self.assertEqual(d['superorder'], 'Euarchontoglires')
        self.assertEqual(d['order'], 'Rodentia')
        self.assertEqual(d['suborder'], 'Sciurognathi')
        self.assertEqual(d['family'], 'Muridae')
        self.assertEqual(d['subfamily'], 'Murinae')
        self.assertEqual(d['genus'], 'Mus')
        self.assertEqual(d['subgenus'], 'Mus')
        os.remove('taxonomy.txt')
        self.assertFalse(os.path.isfile('taxonomy.txt'))

    def test_invalid_bioproject_id_raises_error_in_get_taxonomy(self):
        handle = Entrez.elink(dbfrom="bioproject", db="taxonomy", id='48361')
        link_record = Entrez.read(handle)
        tax_id = link_record[0]['LinkSetDb'][0]['Link'][0]['Id']
        handle = Entrez.efetch(db="taxonomy", id=tax_id)
        tax_record = Entrez.read(handle)
        handle.close()
        for taxon in tax_record[0]['LineageEx']:
            if('Rank' in taxon and taxon['Rank'] != 'no rank'):
                print "%s\t%s" % (taxon['Rank'], taxon['ScientificName'])
                            
        
        self.assertRaises(IndexError, process_seqs.get_taxonomy(os.getcwd(), '48361'))

    def test_genome_to_taxonomy(self):
        handle = Entrez.esearch(db="genome", term='(txid29258[Organism:exp] OR txid35237[Organism:exp]) AND complete[Status] AND "RefSeq" AND genome_pubmed[FILT]', usehistory="y")
        results = Entrez.read(handle)
        webenv = results["WebEnv"]
        query_key = results["QueryKey"]
        count = int(results["Count"])
        fetch_handle = Entrez.esummary(db="genome", retstart=0, retmax=200, webenv=webenv, query_key=query_key)
        genome_ids = Entrez.read(fetch_handle)
        fetch_handle.close()
        fetch_handle = Entrez.elink(dbfrom="genome", db="taxonomy", webenv=webenv, query_key=query_key, usehistory="y")
        link_results = Entrez.read(fetch_handle)
        fetch_handle.close()
        self.assertEqual(len(link_results[0]['IdList']), count)
        tax_ids = link_results[0]['LinkSetDb'][0]['Link']
        self.assertEqual(len(tax_ids), count)
        handle = Entrez.epost(db='taxonomy', id=",".join(tax_ids))
        results = Entrez.read(handle)
        handle.close()
        webenv = results["WebEnv"]
        query_key = results["QueryKey"]
        taxonomy_count = int(results["Count"])
        self.assertEqual(taxonomy_count, count)
        fetch_handle = Entrez.efetch(db='taxonomy', restart=0, retmax=200, webenv=webenv, query_key=query_key)
        taxonomy_records = Entrez.read(fetch_handle)
        d = {}
        for taxon in taxonomy_records:
            if('Rank' in taxon and (taxon['Rank'] != 'no rank')):
                d[taxon['LineageEx']['Rank']] = taxon['LineageEx']['ScientificName']
        
        
if __name__ == '__main__':
    unittest.main()
