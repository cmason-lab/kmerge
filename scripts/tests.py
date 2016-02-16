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
import optunity
from operator import itemgetter
import gzip
from sklearn.preprocessing import LabelEncoder
from gensim.corpora.dictionary import Dictionary
import itertools
from Bio.Seq import Seq
import pyhash
import pandas as pd
from scipy.sparse import csr_matrix
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

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
        if filter_file:
            self.filter = LabelEncoder()
            self.filter.fit(np.load(filter_file)["hashes"])
            pickle.dump(self.filter, open("encoder.k%s.pkl" % self.ks, "wb"))
            #self.filter = {h:i for i,h in enumerate(np.load(filter_file)["hashes"])}

    def __len__(self):
        return len(self.selection)

    def set_selection(self, s):
        self.selection = s
    
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
        self.filter_file = pytest.config.getoption('ff')
        self.groups = list(set([ int(f.split(".")[0]) for f in listdir("%s/%s_%s_%s" % (self.reference, self.hash, self.ks, self.ke)) if isfile(join("%s/%s_%s_%s" % (self.reference,self.hash, self.ks, self.ke),f)) ]))
        self.k_corpus = KMergeCorpus(self.groups, self.reference, self.hash, self.ks, self.ke, None)
        sel = [int(group) for group in open("bacteria_groups.txt").readlines()]
        self.k_corpus.set_selection(sel)
        
    def test_create_idf_values(self):
        tt = models.TfidfModel(self.k_corpus)
        sorted_vals = sorted(tt.idfs.items(), key=itemgetter(1), reverse=True)
        pickle.dump(sorted_vals, gzip.open('idf_vals.%s.%s.pklz' % (self.hash, self.ks), "wb"))
        #np.savez_compressed("hashes.k%s.npz" % self.ks, hashes=np.array([t[0] for t in sorted_vals]))
        
class FunctionalGroupTest(unittest.TestCase):

    def setup_method(self, method):
        self.reference = pytest.config.getoption('reference')
        self.hash = pytest.config.getoption('hash')
        self.ks = pytest.config.getoption('ks')
        self.ke = pytest.config.getoption('ke')
        self.t = pytest.config.getoption('t')
        self.filter_file = pytest.config.getoption('ff')
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
        #params = {"num_topics": [50,500]}
        #cv_f = optunity.cross_validated(x=np.array(range(len(self.groups))), y=np.array(range(len(self.groups))), num_folds=10, regenerate_folds=True)(self.evaluate)
        #optimal_configuration, info, _ = optunity.maximize(self.evaluate, 30, **params)
        #print optimal_configuration
        #print info
        #print "Saving model"
        #model = models.ldamulticore.LdaMulticore(self.k_corpus, num_topics=int(optimal_configuration["num_topics"]), workers=self.t)
        print "Building model"
        d = Dictionary.from_corpus(self.k_corpus)
        model = models.hdpmodel.HdpModel(self.k_corpus, d, T=500)
        #model = models.ldamulticore.LdaMulticore(self.k_corpus, num_topics=479, workers=self.t-1)
        #model = models.ldamodel.LdaModel(self.k_corpus, num_topics=479, passes=100)
        model.save("kmer_hdp.k%s" % self.ks, ignore=['corpus'])
        
        print "Done"

    def test_distributed_lda_model(self):
        print "Building model"
        # 482 reference pathway maps as of 12/17/2015
        model = models.ldamodel.LdaModel(self.k_corpus, num_topics=482, passes=100, distributed=True)
        model.save("kmer_lda.k%s" % self.ks, ignore=['corpus'])
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
        
