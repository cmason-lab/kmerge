import os, argparse, sys, re
import random
import unittest
import process_seqs
from Bio import Entrez
import urllib2
import gzip
from subprocess32 import check_output, CalledProcessError
import shutil
import threadpool
import time

class TestDownloadFastaFunctions(unittest.TestCase):
    
    def setUp(self):
        Entrez.email = 'dar326@cornell.edu'

    #@unittest.skip("skipping")
    def test_fetch_assembly_ids(self):
        ids = process_seqs.fetch_assembly_ids()
        self.assertGreaterEqual(len(ids), 6530) # number of results as of June 21, 2015                                                 
        self.assertIn('285498', ids) #yeast
        #self.assertIn('320101', ids) #human
        self.assertIn('237408', ids) #arabidopsis
        self.assertIn('202931', ids) #drosophila
        self.assertIn('253201', ids) #hpv

    #@unittest.skip("skipping")
    def test_fetch_nucleotide_sequence_ids(self):
        asids = ['253171', '253601', '315421', '285498', '48671', '79781']
        assembly_sequence_ids = process_seqs.fetch_link_ids(asids, "assembly", "nuccore", "assembly_nuccore_refseq")
        self.assertEqual(len(assembly_sequence_ids), len(asids))
        self.assertIn('79781', assembly_sequence_ids) # e. coli
        self.assertIn('556503834', assembly_sequence_ids['79781'])
        self.assertEqual(len(assembly_sequence_ids['79781']), 1)
        self.assertIn('285498', assembly_sequence_ids) # yeast
        self.assertEqual(len(assembly_sequence_ids['285498']), 17)
        self.assertIn('330443753', assembly_sequence_ids['285498'])
        self.assertIn('330443743', assembly_sequence_ids['285498'])
        self.assertIn('330443715', assembly_sequence_ids['285498'])
        self.assertIn('330443688', assembly_sequence_ids['285498'])
        self.assertIn('330443681', assembly_sequence_ids['285498'])
        self.assertIn('330443667', assembly_sequence_ids['285498'])
        self.assertIn('330443638', assembly_sequence_ids['285498'])
        self.assertIn('330443595', assembly_sequence_ids['285498'])
        self.assertIn('330443590', assembly_sequence_ids['285498'])
        self.assertIn('330443578', assembly_sequence_ids['285498'])
        self.assertIn('330443543', assembly_sequence_ids['285498'])
        self.assertIn('330443531', assembly_sequence_ids['285498'])
        self.assertIn('330443520', assembly_sequence_ids['285498'])
        self.assertIn('330443489', assembly_sequence_ids['285498'])
        self.assertIn('330443482', assembly_sequence_ids['285498'])
        self.assertIn('330443391', assembly_sequence_ids['285498'])
        self.assertIn('6226515', assembly_sequence_ids['285498'])

    #@unittest.skip("skipping")
    def test_get_sequence_from_refseq(self):
        d = os.getcwd()
        fasta_handle = gzip.open("%s/sequence.fa.gz" % d, 'w+b')
        acc = 'NC_000913.3'
        db_dir = '/zenodotus/dat01/mason_lab_scratch_reference/cmlab/GENOMES/BLAST/'
        self.assertTrue(process_seqs.get_sequence_from_refseq(fasta_handle, acc, db_dir))
        fasta_handle.close()
        os.remove("%s/sequence.fa.gz" % d)


    #@unittest.skip("skipping")
    def test_full_pipeline(self):
        batch_size = 200
        base = os.getcwd() + "/"
        seq_format = "fasta"
        db_dir = "/zenodotus/dat01/mason_lab_scratch_reference/cmlab/GENOMES/BLAST/"

        pool = threadpool.ThreadPool(1)

        ids = []
        try:
            ids = process_seqs.fetch_assembly_ids(batch_size)
        except Exception as inst:
            sys.stderr.write("Error converting non-virus assembly ids to bioproject ids\n")


        classifications = process_seqs.fetch_classifications(ids[:30])

        self.assertEqual(len(classifications), 30)
        # get the asids that have classifications and ignore rest
        asids = classifications.keys()[:6]
    
        for ncbi_asid in asids:
            request = threadpool.WorkRequest(process_seqs.process_genomes, args=[base, ncbi_asid, classifications, seq_format, db_dir, True, 0, 2])
            pool.putRequest(request)

        pool.wait()
        for asid in asids:
            shutil.rmtree('./%s' % asid)

    #@unittest.skip("skipping")
    def test_get_taxonomy(self):
        assembly_ids = ['315421','284548','375868','587628', '196698','39508','253601','376268','587928','110908','44548','253171','380438',
                        '587648','110928','52131', '31388', '407318', '48671', '360488', '445688', '46691', '253201']
        classifications = process_seqs.fetch_classifications(assembly_ids, 30)
        process_seqs.get_taxonomy(os.getcwd(), '315421', classifications)
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

        process_seqs.get_taxonomy(os.getcwd(), '52131', classifications)
        self.assertTrue(os.path.isfile('taxonomy.txt'))
        d = {}
        with open("taxonomy.txt") as f:
            for line in f:
                (key, val) = line.strip().split("\t")
                d[key] = val
                
        self.assertEqual(d['strain'], 'Mycobacterium tuberculosis PanR0707')                                                                     
        os.remove('taxonomy.txt')
        self.assertFalse(os.path.isfile('taxonomy.txt'))
        
        process_seqs.get_taxonomy(os.getcwd(), '253201', classifications)
        self.assertTrue(os.path.isfile('taxonomy.txt'))
        d = {}
        with open("taxonomy.txt") as f:
            for line in f:
                (key, val) = line.strip().split("\t")
                d[key] = val
    
        self.assertEqual(d['strain'], 'Human papillomavirus type 60')
        os.remove('taxonomy.txt')
        self.assertFalse(os.path.isfile('taxonomy.txt'))
                                                                                                            
        
    #@unittest.skip("skipping")
    def test_assembly_to_taxonomy(self):
        assembly_ids = ['315421','284548','375868','587628', '196698','39508','253601','376268','587928','110908','44548','253171','380438',
                        '587648','110928', '31388', '407318', '48671', '360488', '445688', '46691']

        d = process_seqs.fetch_classifications(assembly_ids)
        missing = set(assembly_ids) - set(d.keys())
        self.assertEqual(len(missing), 0)
        for bp, taxonomy in d.iteritems():
            self.assertTrue('species' in taxonomy)

            
if __name__ == '__main__':
    unittest.main()
