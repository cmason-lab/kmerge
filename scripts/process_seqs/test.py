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

class TestDownloadFastaFunctions(unittest.TestCase):
    
    def setUp(self):
        Entrez.email = 'dar326@cornell.edu'

    #@unittest.skip("skipping")
    def test_fetch_bioproject_ids(self):
        ids_and_stats = process_seqs.fetch_assembly_ids_and_stats()
        self.assertGreaterEqual(len(ids_and_stats), 6530) # number of results as of June 21, 2015                                                 
        self.assertIn('285498', ids_and_stats) #yeast
        self.assertEqual(ids_and_stats['285498']['num_sequences'], 17)
        self.assertEqual(ids_and_stats['285498']['sequence_length'], 12157105) 
        self.assertIn('320101', ids_and_stats) #human
        self.assertEqual(ids_and_stats['320101']['num_sequences'], 25)
        self.assertEqual(ids_and_stats['320101']['sequence_length'], 3226010022)
        self.assertIn('237408', ids_and_stats) #arabidopsis
        self.assertEqual(ids_and_stats['237408']['num_sequences'], 7)
        self.assertEqual(ids_and_stats['237408']['sequence_length'], 119667750)
        self.assertIn('202931', ids_and_stats) #drosophila
        self.assertEqual(ids_and_stats['202931']['num_sequences'], 8)
        self.assertEqual(ids_and_stats['202931']['sequence_length'], 143726002)
        self.assertIn('266321', ids_and_stats) #hiv
        self.assertEqual(ids_and_stats['266321']['num_sequences'], 1)
        self.assertEqual(ids_and_stats['266321']['sequence_length'], 9181)

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

        ids_and_stats = []
        try:
            ids_and_stats = process_seqs.fetch_assembly_ids_and_stats(batch_size)
        except Exception as inst:
            sys.stderr.write("Error converting non-virus assembly ids to bioproject ids\n")


        classifications = process_seqs.fetch_classifications(ids_and_stats.keys()[:6])
        # get the asids that have classifications and ignore rest
        asids = classifications.keys()
    
        for ncbi_asid in asids:
            request = threadpool.WorkRequest(process_seqs.process_genomes, args=[base, ncbi_asid, ids_and_stats[ncbi_asid], classifications, seq_format, db_dir, True, 0, 2])
            pool.putRequest(request)

        pool.wait()
        for asid in asids:
            shutil.rmtree('./%s' % asid)

    #@unittest.skip("skipping")
    def test_get_taxonomy(self):
        assembly_ids = ['315421','284548','375868','587628', '196698','39508','253601','376268','587928','110908','44548','253171','380438','587648','110928',
                                                    '31388', '407318', '48671', '360488', '445688', '46691']
        classifications = process_seqs.fetch_classifications(assembly_ids)
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

    #@unittest.skip("skipping")
    def test_genome_to_taxonomy(self):
        assembly_ids = ['315421','284548','375868','587628', '196698','39508','253601','376268','587928','110908','44548','253171','380438','587648','110928',
                          '31388', '407318', '48671', '360488', '445688', '46691']

        d = process_seqs.fetch_classifications(assembly_ids)
        self.assertEqual(len(assembly_ids), len(d))
        for bp, taxonomy in d.iteritems():
            self.assertTrue('species' in taxonomy)
        
if __name__ == '__main__':
    unittest.main()
