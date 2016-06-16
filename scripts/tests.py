import pytest
import sys
sys.path.append("/home/darryl/Development/kmerge/scripts/loader")
import loader
import numpy as np 
from os import listdir
import os
import os.path
from os.path import isfile, join
import unittest
import cPickle as pickle
from gensim.models.tfidfmodel import TfidfModel
from gensim import models
from gensim.corpora.dictionary import Dictionary
import logging
from operator import itemgetter
import gzip
from sklearn.preprocessing import LabelEncoder
from gensim.corpora.dictionary import Dictionary
import itertools
from Bio.Seq import Seq
#import pyhash
import pandas as pd
from scipy.sparse import csr_matrix
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import graphlab as gl
from graphlab import SArray
from pyspark import SparkContext, SparkConf

logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO)

class KMergeCorpus(object):
    def __init__(self, groups, reference, hashf, ks, ke, filter_file=None):
        self.groups = groups
        self.reference = reference
        self.hash = hashf
        self.ks = ks
        self.ke = ke
        self.selection = None
        self.filter = None
        if filter_file or isfile("hashes.k%s.npz" % ks):
            self.filter = LabelEncoder()
            self.filter.fit(np.load(filter_file)["hashes"])
            pickle.dump(self.filter, open("encoder.k%s.pkl" % self.ks, "wb"))
            #self.filter = {h:i for i,h in enumerate(np.load(filter_file)["hashes"])}

    def __len__(self):
        return len(self.selection)

    def set_selection(self, s):
        self.selection = s


    def store(self):
        # build as sparse matrix first
        self.sa = SArray(dtype=dict)
        for i, group in enumerate(self.selection):
            #if self.selection and i not in self.selection:
            #if self.selection and group not in self.selection:
            #    continue
            counts_file = "%s/%s_%s_%s/%s.counts.bin" % (self.reference, self.hash, self.ks, self.ke, group)
            hashes_file = "%s/%s_%s_%s/%s.hashes.bin" % (self.reference, self.hash, self.ks, self.ke, group)
            c = loader.create_list(counts_file,False)
            h = loader.create_list(hashes_file,True)
            print "%s, %s, %s" % (group, len(h), len(c))
            values = {}
            if self.filter != None:
                mask = np.in1d(h, self.filter.classes_)
                values = dict(zip(self.filter.transform(np.array(h, dtype=np.uint32)[mask]), np.array(c, dtype=np.uint32)[mask]))
            else:
                values = dict(zip(h, c))
            if len(values) == 0:
                print "No data for %s" % group
                continue
            sa = SArray([values])
            self.sa = self.sa.append(sa)
        
    def __iter__(self):
        for i, group in enumerate(self.selection):
            #if self.selection and i not in self.selection:
            #if self.selection and group not in self.selection:
            #    continue
            counts_file = "%s/%s_%s_%s/%s.counts.bin" % (self.reference, self.hash, self.ks, self.ke, group)
            hashes_file = "%s/%s_%s_%s/%s.hashes.bin" % (self.reference, self.hash, self.ks, self.ke, group)
            c = loader.create_list(counts_file,False)
            h = loader.create_list(hashes_file,True)
            values = []
            if self.filter != None:
                mask = np.in1d(h, self.filter.classes_)
                values = zip(self.filter.transform(np.array(h)[mask]), np.array(c)[mask])
            else:
                values = zip(h, c)
            if len(values) == 0:
                print "No data for %s" % group
            yield values



class OnlineIDFFeatureSelectionTest(unittest.TestCase):
    
    def setup_method(self, method):
        self.reference = pytest.config.getoption('reference')
        self.hash = pytest.config.getoption('hash')
        self.ks = pytest.config.getoption('ks')
        self.ke = pytest.config.getoption('ke')
        self.t = pytest.config.getoption('t')
        self.filter_file = pytest.config.getoption('ffh')
        self.groups = list(set([ int(f.split(".")[0]) for f in listdir("%s/%s_%s_%s" % (self.reference, self.hash, self.ks, self.ke)) if isfile(join("%s/%s_%s_%s" % (self.reference,self.hash, self.ks, self.ke),f)) ]))
        self.k_corpus = KMergeCorpus(self.groups, self.reference, self.hash, self.ks, self.ke, None)
        sel = [int(group) for group in open("bacteria_groups.txt").readlines()]
        self.k_corpus.set_selection(sel)
        
    def test_create_idf_values(self):
        tt = models.TfidfModel(self.k_corpus)
        sorted_vals = sorted(tt.idfs.items(), key=itemgetter(1), reverse=True)
        pickle.dump(sorted_vals, gzip.open('idf_vals.%s.%s.pklz' % (self.hash, self.ks), "wb"))
        np.savez_compressed("hashes.k%s.npz" % self.ks, hashes=np.array([t[0] for t in sorted_vals]))
        #hist, bin_edges = np.histogram([sorted_vals[i][1] for i in xrange(0,len(sorted_vals))], bins=100)
        #hashes=[sorted_vals[i][0] for i in xrange(0,len(sorted_vals)) if sorted_vals[i][1] < bin_edges[1] or sorted_vals[i][1] > bin_edges[-2]]

        
class FunctionalGroupTest(unittest.TestCase):

    def setup_method(self, method):
        self.reference = pytest.config.getoption('reference')
        self.hash = pytest.config.getoption('hash')
        self.ks = pytest.config.getoption('ks')
        self.ke = pytest.config.getoption('ke')
        self.t = pytest.config.getoption('t')
        self.filter_file = pytest.config.getoption('ffh')
        self.groups = list(set([ int(f.split(".")[0]) for f in listdir("%s/%s_%s_%s" %
                                 (self.reference, self.hash, self.ks, self.ke)) if isfile(join("%s/%s_%s_%s" % (self.reference,self.hash, self.ks, self.ke),f)) ]))
        self.k_corpus = KMergeCorpus(self.groups, self.reference, self.hash, self.ks, self.ke, self.filter_file)
        sel = [int(group) for group in open("bacteria_groups.txt").readlines()]
        self.k_corpus.set_selection(sel)
                

        
    def evaluate(this, num_topics=None):
        print "Trying num_topics = %s" % int(num_topics)
        # evaluate on either log_perplexity or UMass topic coherence (top_topics function)
        model = models.ldamulticore.LdaMulticore(this.k_corpus, num_topics=int(num_topics), workers=this.t)
        return model.top_topics(this.k_corpus, num_words=20)
        
    def test_hdp_model(self):
        print "Building model"
        d = Dictionary.from_corpus(self.k_corpus)
        model = models.hdpmodel.HdpModel(self.k_corpus, d, T=500)
        #model = models.ldamulticore.LdaMulticore(self.k_corpus, num_topics=479, workers=self.t-1)
        #model = models.ldamodel.LdaModel(self.k_corpus, num_topics=479, passes=100)
        model.save("kmer_hdp.k%s" % self.ks, ignore=['corpus'])
        
        print "Done"

    def test_distributed_lda_model(self):
        print "Building model"
        conf = SparkConf().setAppName("test").setMaster("local[4]")
        self.k_corpus.store()
        # 484 reference pathway maps as of 03/05/2016
        #model = gl.topic_model.create(self.k_corpus.sa, num_topics=484, num_iterations=100)
        model = gl.topic_model.create(self.k_corpus.sa, num_topics=484, num_iterations=10)
        model.save("kmer_lda.gl.k%s" % self.ks)
        print "Done"


    def test_distributed_lda_model_alias(self):
        print "Building model with alias lda"
        self.k_corpus.store()
        # 484 reference pathway maps as of 03/05/2016
        model = gl.topic_model.create(self.k_corpus.sa, num_topics=484, num_iterations=10, method='alias')
        model.save("kmer_lda.gl.alias.k%s" % self.ks)
        print "Done"
                                                        
        
    def test_lda_model(self):
        print "Building model"
        # 482 reference pathway maps as of 12/17/2015
        model = models.ldamodel.LdaModel(self.k_corpus, num_topics=482, passes=5)
        del model.corpus
        model.save("kmer_lda.k%s" % self.ks)
        print "Done"

        
        
class PrepareKMerTopicsForAssemblyTest(unittest.TestCase):

    def setup_method(self, method):
        self.reference = pytest.config.getoption('reference')
        self.hash = pytest.config.getoption('hash')
        self.ks = pytest.config.getoption('ks')
        self.ke = pytest.config.getoption('ke')
        self.t = pytest.config.getoption('t')
        
    def test_functionality(self):
        model = models.hdpmodel.HdpModel.load("kmer_hdp.k11")
        encoder = pickle.load(open("encoder.k%s.pkl" % self.ks, "rb"))
        bases=['A','T','G','C']
        kmers=[''.join(p) for p in itertools.product(bases, repeat=self.ks)]
        
        hasher = pyhash.murmur3_32()
        kmer_map = {}
        for kmer in kmers:
            rc_kmer = str(Seq(kmer).reverse_complement())
            combine = ""
            if kmer < rc_kmer:
                combine = kmer + rc_kmer
            else:
                combine = rc_kmer + kmer
            kmer_map.setdefault(hasher(combine), set())
            kmer_map[hasher(combine)].add(kmer if kmer < rc_kmer else rc_kmer)

            
            
        # determine probability cutoff for including kmer in topic
        df = pd.DataFrame([{k:v for k,v in topic[1]} for topic in model.show_topics(topics=-1, topn=-1, formatted=False)])
        # get top 25 percent of kmers
        qt = df.quantile(q=0.25, axis=1)
        # replace this with whatever cutoff is chosen
        df = df.apply(lambda row: row.where(row >= qt[int(row.name)]), axis=1)
        df.fillna(0, inplace=True)
        topics = csr_matrix(df)

        for row in xrange(topics.shape[0]):
            sequences = []
            output_handle = open("assemble/hdp_k11/topic_%s.fasta" % row, "w")
            for col in topics[row,:].nonzero()[1]:
                k_idx = encoder.inverse_transform(int(df.columns.values[col]))
                for i, kmer in enumerate(kmer_map[k_idx]):
                    # create k-mer twice for velvet
                    sequences.append(SeqRecord(Seq(kmer),id='%s_%s' % (k_idx,i), name='', description=''))
            SeqIO.write(sequences, output_handle, "fasta")
            output_handle.close()
        
