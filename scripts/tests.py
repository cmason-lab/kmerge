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

logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO)

class KMergeCorpus(object):
    def __init__(self, groups, reference, hashf, ks, ke):
        self.groups = groups
        self.reference = reference
        self.hash = hashf
        self.ks = ks
        self.ke = ke
        self.selection = None

    def __len__(self):
        return len(self.groups)

    def set_selection(self, s):
        self.selection = s
    
    def __iter__(self):
        for i, group in enumerate(self.groups):
            if self.selection and i not in self.selection:
                continue
            counts_file = "%s/%s_%s_%s/%s.counts.bin" % (self.reference, self.hash, self.ks, self.ke, group)
            hashes_file = "%s/%s_%s_%s/%s.hashes.bin" % (self.reference, self.hash, self.ks, self.ke, group)
            c = loader.create_list(counts_file,False)
            h = loader.create_list(hashes_file,True)
            values = []
            for i in range(0, len(h)):
                values.append((h[i], c[i]))
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
        self.groups = list(set([ int(f.split(".")[0]) for f in listdir("%s/%s_%s_%s" % (self.reference, self.hash, self.ks, self.ke)) if isfile(join("%s/%s_%s_%s" % (self.reference,self.hash, self.ks, self.ke),f)) ]))
        self.k_corpus = KMergeCorpus(self.groups, self.reference, self.hash, self.ks, self.ke)

    def test_create_idf_values(self):
        tt = models.TfidfModel(self.k_corpus)
        sorted_vals = sorted(tt.idfs.items(), key=itemgetter(1), reverse=True)
        pickle.dump(sorted_vals, gzip.open('idf_vals.%s.%s.pklz' % (self.hash, self.ks), "wb"))

class FunctionalGroupTest(self):

    def setup_method(self, method):
        self.reference = pytest.config.getoption('reference')
        self.hash = pytest.config.getoption('hash')
        self.ks = pytest.config.getoption('ks')
        self.ke = pytest.config.getoption('ke')
        self.t = pytest.config.getoption('t')
        self.groups = list(set([ int(f.split(".")[0]) for f in listdir("%s/%s_%s_%s" %
                                 (self.reference, self.hash, self.ks, self.ke)) if isfile(join("%s/%s_%s_%s" % (self.reference,self.hash, self.ks, self.ke),f)) ]))
        self.k_corpus = KMergeCorpus(self.groups, self.reference, self.hash, self.ks, self.ke)

    def evaluate(this, x_train, y_train, x_test, y_test, num_topics=None):
        # evaluate on either log_perplexity or UMass topic coherence (top_topics function)
        corpus = KMergeCorpus(self.groups, self.reference, self.hash, self.ks, self.ke)
        corpus.set_selection(x_train.to_list())
        model = models.ldamulticore.LdaMulticore(corpus, num_topics=num_topics, workers=self.t)
        # should it be evaluated on the test documents or entire corpus?
        corpus.set_selection(x_test.to_list())
        return model.top_topics(corpus)
        
    def test_lda_model(self):
        params = {"num_topics": [50,500]}
        cv_f = optunity.cross_validated(x=np.array(range(len(self.groups))), y=np.array(range(len(self.groups))), num_folds=10, regenerate_folds=True)(self.evaluate)
        optimal_configuration, info, _ = optunity.maximize(cv_f, 300, pmap=optunity.parallel.create_pmap(self.t-1), **params)
                                                
