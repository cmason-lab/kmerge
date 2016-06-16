import sys
sys.path.append("/home/darryl/Development/kmerge/scripts/loader")
import loader
import subprocess
import os
import numpy as np
from scipy.stats import entropy
from scipy.sparse import csr_matrix

for dir in os.listdir("."):
    fp = open("%s/entropy.txt" % dir, "wb")
    fp.write("type\tk\tentropy\tKL divergence\tredundancy\n")
    i = 0
    for d_name, d_type in dict({"8bit":np.uint8, "16bit":np.uint16, "32bit":np.uint32}).items():
        for k in range(3,31,2):
            if i == 0:
                if k <= 17:
                    # determine entropy of unhashed k-mers
                    c = []
                    for l in subprocess.check_output(['jellyfish', 'dump', '-c', '/zenodotus/masonlab/darryl_scratch/information_theory/%s/raw/k%s.jf' % (dir, k)]).split("\n"):
                        if len(l):
                            c.append(np.uint32(l.split()[1]))
                    max_size = 4**k/2
                    pad_length = max_size - len(c)
                    ent = entropy(c if pad_length else c + [0]*pad_length, base=2)
                    if k > 3:
                        kl = entropy(c if pad_length else c + [0]*pad_length, c_last, base=2)
                    else:
                        kl = np.nan
                    fp.write("~hashed\t%d\t%f\t%f\t%f\n" % (k, ent, kl, (1-(ent/np.log2(max_size)))))
                    c_last = np.array(c if pad_length else c + [0]*pad_length)
                else:
                    fp.write("~hashed\t%s\tnp.nan\tnp.nan\tnp.nan\n" % k)
            
            c = np.array((loader.create_list("/zenodotus/masonlab/darryl_scratch/information_theory/%s/murmur/k%s/%s.counts.bin" % (dir,k,dir),False)), dtype=np.uint32)
            h = np.array((loader.create_list("/zenodotus/masonlab/darryl_scratch/information_theory/%s/murmur/k%s/%s.hashes.bin" % (dir,k,dir),True)), dtype=d_type)
            m = csr_matrix((c, ([0]*len(h),h)), shape=(1,np.iinfo(d_type).max+1))
            ent = entropy(m.toarray()[0], base=2)
            if k > 3:
                kl = entropy(m.toarray()[0], m_last.toarray()[0], base=2)
            else:
                kl = np.nan
            fp.write("%s\t%d\t%f\t%f\t%f\n" % (np.typename(np.sctype2char(d_type)), k, ent, kl, (1-ent/np.log2(np.iinfo(d_type).max+1))))
            m_last = csr_matrix(m)
        i += 1
    fp.close()
             
