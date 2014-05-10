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

    @unittest.skip("skipping")
    def test_fetch_non_virus_bioproject_ids(self):
        pass

    @unittest.skip("skipping")
    def test_fetch_virus_bioproject_ids(self):
        pass

    @unittest.skip("skipping")
    def test_fetch_taxonomy_details(self):
        pass

    @unittest.skip("skipping")
    def test_get_taxonomy(self):
        process_seqs.get_taxonomy(os.getcwd(), '169')
        self.assertTrue(os.path.isfile('taxonomy.txt'))
        d = {}
        with open("taxonomy.txt") as f:
            for line in f:
                (key, val) = line.strip().split("\t")
                d[key] = val
        self.assertEqual(len(d), 13)
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
        self.assertEqual(d['species'], 'Mus musculus')
        os.remove('taxonomy.txt')
        self.assertFalse(os.path.isfile('taxonomy.txt'))

    @unittest.skip("skipping")
    def test_invalid_bioproject_id_raises_error_in_get_taxonomy(self):        
        self.assertRaises(IndexError, process_seqs.get_taxonomy, os.getcwd(), '209158')

    #@unittest.skip("skipping")
    def test_genome_to_taxonomy(self):
        bioproject_ids = ['168','14003','162087','196786', '47493','59067','14013','162089','196787','47507','59071','14014','162091','196788','47509',
                          '59073']

        d = process_seqs.fetch_classifications(bioproject_ids)

        for bp, taxonomy in d.iteritems():
            self.assertTrue('species' in taxonomy)
        
if __name__ == '__main__':
    unittest.main()
