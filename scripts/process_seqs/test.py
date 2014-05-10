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
        ids = process_seqs.fetch_non_virus_bp_ids()
        self.assertGreaterEqual(len(ids), 2292) # number of results as of May 10, 2014                                                 
        self.assertIn('128', ids) #yeast
        self.assertIn('168', ids) #human
        self.assertIn('116', ids) #arabidopsis
        self.assertIn('164', ids) #drosophila
        self.assertIn('158', ids) #c. elegans

    @unittest.skip("skipping")
    def test_fetch_virus_bioproject_ids(self):
        ids = process_seqs.fetch_virus_bp_ids()

        self.assertGreaterEqual(len(ids), 761) # number of results as of May 10, 2014
        self.assertIn('218024', ids)
        self.assertIn('15423', ids)
        self.assertIn('14331', ids)

    @unittest.skip("skipping")
    def test_fetch_nucleotide_sequence_ids(self):
        pjids = ['169','14003','162087','196786', '47493','59067','14013','162089','196787','47507','59071','14014','162091','196788','47509',
                                  '59073', '86861', '215233', '88071', '86645', '89395', '213395', '225', '128']
        project_sequence_ids = process_seqs.fetch_link_ids(pjids, "bioproject", "nuccore")
        self.assertEqual(len(project_sequence_ids), len(pjids))
        self.assertIn('225', project_sequence_ids) # e. coli
        self.assertIn('545778205', project_sequence_ids['225'])
        self.assertEqual(len(project_sequence_ids['225']), 1)
        self.assertIn('128', project_sequence_ids) # yeast
        self.assertEqual(len(project_sequence_ids['128']), 17)
        self.assertIn('330443753', project_sequence_ids['128'])
        self.assertIn('330443743', project_sequence_ids['128'])
        self.assertIn('330443715', project_sequence_ids['128'])
        self.assertIn('330443688', project_sequence_ids['128'])
        self.assertIn('330443681', project_sequence_ids['128'])
        self.assertIn('330443667', project_sequence_ids['128'])
        self.assertIn('330443638', project_sequence_ids['128'])
        self.assertIn('330443595', project_sequence_ids['128'])
        self.assertIn('330443590', project_sequence_ids['128'])
        self.assertIn('330443578', project_sequence_ids['128'])
        self.assertIn('330443543', project_sequence_ids['128'])
        self.assertIn('330443531', project_sequence_ids['128'])
        self.assertIn('330443520', project_sequence_ids['128'])
        self.assertIn('330443489', project_sequence_ids['128'])
        self.assertIn('330443482', project_sequence_ids['128'])
        self.assertIn('330443391', project_sequence_ids['128'])
        self.assertIn('6226515', project_sequence_ids['128'])
        
    @unittest.skip("skipping")
    def test_get_taxonomy(self):
        bioproject_ids = ['169','14003','162087','196786', '47493','59067','14013','162089','196787','47507','59071','14014','162091','196788','47509',
                                                    '59073', '86861', '215233', '88071', '86645', '89395', '213395']
        classifications = process_seqs.fetch_classifications(bioproject_ids)
        process_seqs.get_taxonomy(os.getcwd(), '169', classifications)
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
                          '59073', '86861', '215233', '88071', '86645', '89395', '213395']

        d = process_seqs.fetch_classifications(bioproject_ids)
        self.assertEqual(len(bioproject_ids), len(d))
        for bp, taxonomy in d.iteritems():
            self.assertTrue('species' in taxonomy)
        
if __name__ == '__main__':
    unittest.main()
