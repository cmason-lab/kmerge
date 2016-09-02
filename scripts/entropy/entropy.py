import sys
sys.path.append("/home/darryl/Development/kmerge/scripts/loader")
import loader
import subprocess
import os
import numpy as np
from scipy.stats import entropy
from scipy.sparse import csr_matrix
from numpy.linalg import norm


def JSD(P, Q):
    _P = P / norm(P, ord=1)
    _Q = Q / norm(Q, ord=1)
    _M = 0.5 * (_P + _Q)
    return 0.5 * (entropy(_P, _M, base=2) + entropy(_Q, _M, base=2))


for dir in os.listdir("."):
    fp = open("%s/entropy.txt" % dir, "wb")
    fp.write("type\tk\tentropy\n")
    for k in range(3,33,2):
        # determine entropy of unhashed k-mers
        c = []
        for l in subprocess.check_output(['jellyfish', 'dump', '-c', '/zenodotus/masonlab/darryl_scratch/information_theory/%s/raw/k%s.jf' % (dir, k)]).split("\n"):
            if len(l):
                c.append(np.uint32(l.split()[1]))
        max_size = 4**k/2
        pad_length = max_size - len(c)
        ent = entropy(c, base=2)
        fp.write("~hashed\t%d\t%f\n" % (k, ent))
            
        c = np.array((loader.create_list("/zenodotus/masonlab/darryl_scratch/information_theory/%s/murmur/k%s/%s.counts.bin" % (dir,k,dir),False)), dtype=np.uint32)
        for d_name, d_type in dict({"8bit":np.uint8, "16bit":np.uint16, "32bit":np.uint32}).items():
            h = np.array((loader.create_list("/zenodotus/masonlab/darryl_scratch/information_theory/%s/murmur/k%s/%s.hashes.bin" % (dir,k,dir),True)), dtype=d_type)
            m = csr_matrix((c, ([0]*len(h),h)), shape=(1,np.iinfo(d_type).max+1))
            #ent = entropy(m.toarray()[0], base=2)
            ent = entropy(m.data, base=2)
            fp.write("%s\t%d\t%f\n" % (d_name, k, ent))

    fp.close()
             
