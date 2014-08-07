import unittest
#import reorg_data
import tables
import numpy as np
import os
import subprocess as sbp

class TestReorganizeData(unittest.TestCase):
    
    def setUp(self):
        pass
   
    def test_convert_column_counts_to_matrix_counts(self):
        shape = (5,0)
        atom = tables.UInt32Atom()
        filters = tables.Filters(complevel=9, complib='zlib')
        h5fh = tables.open_file("earray1.h5", mode='a', filters=filters)
        ea = h5fh.create_earray(h5fh.root, 'counts', atom, shape, "counts matrix", filters, 2**32-1)
        self.assertEqual(ea.shape[1], 0)
        col1 = np.array([1,6,11,16,21], dtype=np.uint32, ndmin=2).transpose()
        ea.append(col1)
        self.assertEqual(ea.shape[1], 1)
        col2 = np.array([2,7,12,17,22], dtype=np.uint32, ndmin=2).transpose()
        ea.append(col2)
        self.assertEqual(ea.shape[1], 2)
        col3 = np.array([3,8,13,18,23], dtype=np.uint32, ndmin=2).transpose()
        ea.append(col3)
        self.assertEqual(ea.shape[1], 3)
        col4 = np.array([4,9,14,19,24], dtype=np.uint32, ndmin=2).transpose()
        ea.append(col4)
        self.assertEqual(ea.shape[1], 4)
        col5 = np.array([5,10,15,20,25], dtype=np.uint32, ndmin=2).transpose()
        ea.append(col5)
        self.assertEqual(ea.shape[1], 5)
        print h5fh.root.counts[:]
        self.assertTrue(np.array_equal(h5fh.root.counts[:,0], col1[:,0]))
        self.assertTrue(np.array_equal(h5fh.root.counts[:,1], col2[:,0]))
        self.assertTrue(np.array_equal(h5fh.root.counts[:,2], col3[:,0]))
        self.assertTrue(np.array_equal(h5fh.root.counts[:,3], col4[:,0]))
        self.assertTrue(np.array_equal(h5fh.root.counts[:,4], col5[:,0]))
        h5fh.close()
        os.remove("earray1.h5")

    def test_append_to_matrix_counts(self):
        shape = (5,0)
        atom = tables.UInt32Atom()
        filters = tables.Filters(complevel=9, complib='zlib')
        h5fh = tables.open_file("earray_append.h5", mode='a', filters=filters)
        ea = h5fh.create_earray(h5fh.root, 'counts', atom, shape, "counts matrix", filters, 2**32-1)
        self.assertEqual(ea.shape[1], 0)
        col1 = np.array([1,6,11,16,21], dtype=np.uint32, ndmin=2).transpose()
        ea.append(col1)
        h5fh.close()

        # append to file
        h5fh = tables.open_file("earray_append.h5", mode='a', filters=filters)
        self.assertIn('counts', h5fh.root)
        counts =  h5fh.root.counts
        self.assertEqual(counts.shape[1], 1)
        col2 = np.array([2,7,12,17,22], dtype=np.uint32, ndmin=2).transpose()
        counts.append(col2)
        self.assertEqual(counts.shape[1], 2)
        self.assertTrue(np.array_equal(h5fh.root.counts[:,0], col1[:,0]))
        self.assertTrue(np.array_equal(h5fh.root.counts[:,1], col2[:,0]))
        h5fh.close()
        os.remove("earray_append.h5")


    def test_group_counts_can_be_added_to_matrix_counts(self):
        sbp.check_call(["python", "reorg_data.py", "--master_file", "reference.h5", "--contrib_file", "/home/darryl/Development/kmerge/tests/sandbox/reference.h5"])

        h5fh = tables.open_file("reference.h5", mode='a')
        grouph5fh = tables.open_file("/home/darryl/Development/kmerge/tests/sandbox/reference.h5", mode='r')
        for group_id in list(h5fh.root._v_groups):
            group_num, col_num = group_id.split("|")
            path = "/%s" % group_num
            col = np.resize(grouph5fh.getNode(path, 'count').read(), (2**32-1,1))
            self.assertTrue(np.array_equal(h5fh.root.counts[:,int(col_num)], col[:,0]))	    	
        
        grouph5fh.close()
        h5fh.close()
        os.remove("reference.h5")

if __name__ == '__main__':
    unittest.main()
