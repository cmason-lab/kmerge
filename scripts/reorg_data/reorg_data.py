import tables
import numpy as np
import argparse

def groups_to_matrix(m_file, c_file):
    filters = tables.Filters(complevel=9, complib='zlib')
    h5fh = tables.open_file(m_file, mode='a', filters=filters)
    if not 'counts' in h5fh.root:
        atom = tables.UInt32Atom()
        shape = (2**32-1,0)
        h5fh.create_earray(h5fh.root, 'counts', atom, shape, "counts matrix", expectedrows=2**32-1)
    counts =  h5fh.root.counts

    grouph5fh = tables.open_file(c_file, mode='r')
    for group_num in list(grouph5fh.root._v_groups):
        path = "/%s" % group_num
        print "Processing counts for %s" % path
        counts.append(np.resize(grouph5fh.getNode(path, 'count').read(), (2**32-1,1)))
        col_num = counts.shape[1]-1
        new_group = h5fh.create_group(h5fh.root, "%s|%s" % (group_num, col_num))
        print "Adding taxonomy data for %s" % path
        h5fh.copy_node(grouph5fh.getNode(path, 'taxonomy'), new_group, 'taxonomy', recursive=True)
        print "Finished processing %s" % path
    grouph5fh.close()
    h5fh.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Test of RIQO method for species identification')
    parser.add_argument('--master_file', dest='m_file', help='file of kmer hash counts in matrix form', default='reference_1.h5')
    parser.add_argument('--contrib_file', dest='c_file', help='file of kmer hash counts organized by group', default='reference.h5')
    args = parser.parse_args()

    groups_to_matrix(args.m_file, args.c_file)
