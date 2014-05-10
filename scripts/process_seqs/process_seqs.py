#!/usr/bin/python
import os, argparse, sys, re
from Bio import Entrez
import threadpool
import urllib2
import gzip
from subprocess32 import check_output, CalledProcessError
import shutil

class FileType(argparse.FileType):
    def __call__(self, string):
        if string.endswith('.gz'):
            return gzip.open(string, 'ab')
        return super(FileType, self).__call__(string)

    
def remove_ambiguous_bases(file_handle, output_filename, seq_format, pat):
    id = ""
    first = True
    sequences = []
    seq_num = 0
       
    out_file = gzip.open(output_filename, 'wb')

    for line in file_handle:
        line = line.strip()
        if not line:
            # some fasta files contain empty lines
            return
        if line[0] == ">":
            id = line
            seq_num = 0
            out_file.write( ("" if first else "\n") + line + "/" + str(seq_num) + "\n")
            seq_num = seq_num+1
            first = False
        else:
            sequences = re.split(pat, line)
            if sequences[-1] == '':
                #algorithm requires that split only return empty strings at front of array
                sequences.pop()
            if sequences: # whole line could be N
                out_file.write("%s" % sequences.pop(0))
            for sequence in sequences:
                out_file.write("\n" + id + "/" + str(seq_num) + "\n" + sequence)
                seq_num = seq_num + 1
    out_file.write("\n")
    out_file.close()

def get_taxonomy(d, ncbi_pjid, classifications): # will raise KeyError if key not present
    tax_record = classifications[ncbi_pjid]
    tax_file = open("%s/taxonomy.txt" % d, 'w')
    for rank, name in tax_record.iteritems():
        tax_file.write("%s\t%s\n" % (rank, name))
    tax_file.close()

def fetch_classifications(bioproject_ids, batch_size=20):
    if batch_size < 20:
        batch_size = 20 # minimum size for NCBI records
    lookup = fetch_link_ids(bioproject_ids, "bioproject", "taxonomy", batch_size)
    tax_ids = [tax_id for sublist in lookup.values() for tax_id in sublist]
    bp_ids = [bp_id for bp_id in lookup.keys()]
    # reverse lookup keys and values
    lookup = dict(zip(tax_ids, bp_ids))
    post_handle = Entrez.epost(db='taxonomy', id=",".join(tax_ids), usehistory="y")
    results = Entrez.read(post_handle)
    post_handle.close()
    webenv = results["WebEnv"]
    query_key = results["QueryKey"]
    taxonomy_records = []
    for start in range(0, len(tax_ids), batch_size):
        fetch_handle = Entrez.efetch(db='taxonomy', retstart=start, retmax=batch_size, webenv=webenv, query_key=query_key)
        taxonomy_records.extend(Entrez.read(fetch_handle))
    d = dict.fromkeys(bp_ids)    
    for record in taxonomy_records:
        bp_taxonomy = {}
        species_found = False
        for taxon in record['LineageEx']:
            if('Rank' in taxon and (taxon['Rank'] != 'no rank')):
                if taxon['Rank'] == 'species':
                    species_found = True
                bp_taxonomy[taxon['Rank']] = taxon['ScientificName']
        if not species_found:
            bp_taxonomy[record['Rank']] = record['ScientificName']
        d[lookup[record['TaxId']]] = bp_taxonomy
    return d

def fetch_non_virus_bp_ids(batch_size=20):
    if batch_size < 20:
        batch_size = 20
    pjids = []
    handle = Entrez.esearch(db="assembly", term='(chromosome[ASLV] OR "Gapless Chromosome"[ASLV] OR "Chromosome with gaps"[ASLV]) AND full-genome-representation[Property] AND assembly_nuccore_refseq[FILT] AND assembly_pubmed[FILT] AND (latest[Property] OR latest_refseq[Property] OR latest_genbank[Property]) NOT suppressed_refseq[Property]', usehistory="y")
    results = Entrez.read(handle)
    webenv = results["WebEnv"]
    query_key = results["QueryKey"]
    count = int(results["Count"])
    for start in range(0,count,batch_size):
        sum_handle = Entrez.esummary(db="assembly", retstart=start, retmax=batch_size, webenv=webenv, query_key=query_key)
        summaries = Entrez.read(sum_handle,validate=False)
        pjids.extend([ summary['RS_BioProjects'][0]['BioprojectId'] for summary in summaries['DocumentSummarySet']['DocumentSummary']])
    return pjids

def fetch_virus_bp_ids(batch_size=20):
    if batch_size < 20:
        batch_size = 20
    handle = Entrez.esearch(db="genome", term='(txid29258[Organism:exp] OR txid35237[Organism:exp]) AND complete[Status] AND RefSeq[All Fields] AND genome_pubmed[FILT]', usehistory="y")
    results = Entrez.read(handle)
    webenv = results["WebEnv"]
    query_key = results["QueryKey"]
    count = int(results["Count"])
    summaries = []
    for start in range(0,count,batch_size):
        sum_handle = Entrez.esummary(db="genome", retstart=start, retmax=batch_size, webenv=webenv, query_key=query_key)
        summaries.extend(Entrez.read(sum_handle))
    
    pjids = [ summary['ProjectID'] for summary in summaries]
    return pjids

def fetch_link_ids(lookup_ids, from_db, to_db, split_size=20):
    if split_size < 20:
        split_size = 20
    results = []
    for split_ids in [ lookup_ids[i:i+split_size] for i in range(0, len(lookup_ids), split_size) ]:
        handle = Entrez.elink(dbfrom=from_db, db=to_db, id=split_ids)
        results.extend(Entrez.read(handle))
        handle.close()

    id_map = {}
    for elem in results:
        id_map[elem['IdList'][0]] = [ rec['Id'] for rec in elem['LinkSetDb'][0]['Link']]

    return id_map
                                                                                                    

def process_genomes(base, record, seq_format, pat, db_dir, retry=0, max_retry=2):
    ncbi_pjid = record['IdList'][0]
    d = base + ncbi_pjid    
    if retry > max_retry:
        sys.stderr.write("!%s!: Could not obtain sequences\n" % ncbi_pjid)
        shutil.rmtree(d)
        return
    try:
        os.makedirs(d)
    except OSError:
        if not os.path.isdir(d):
            raise
    seq_count = None
    results = None
    try:
        handle = Entrez.esearch(db="nucleotide", term='%s[BioProject] AND biomol_genomic[PROP] AND srcdb_refseq[PROP] NOT ALTERNATE_LOCUS[Keyword] NOT FIX_PATCH[Keyword] NOT NOVEL_PATCH[Keyword] NOT CONTIG[Title] NOT SCAFFOLD[Title]' % ncbi_pjid, usehistory='y')
        results = Entrez.read(handle)
    except Exception as inst:
        sys.stderr.write("!%s!\n" % ncbi_pjid)
        sys.stderr.write("%s\n" % type(inst))
        return
    webenv = results["WebEnv"]
    query_key = results["QueryKey"]
    count = int(results["Count"])
    batch_size = 100                                                                                                    
    if not retry:
        # write out classifications to taxonomy file
        try:
            get_taxonomy(d, ncbi_pjid)
        except Exception as inst:
            sys.stderr.write("!%s!: Error getting taxonomy information\n" % ncbi_pjid)
            # delet directory if can't get taxonomy
            shutil.rmtree(d)
            return
        results = None
        fasta_handle = open("%s/sequence.fa" % d, 'w+')
        seq_count = 0
        for start in range(0,count,batch_size):
            results = None
            try:
                handle = Entrez.efetch(db="nuccore", retstart=start, retmax=batch_size, webenv=webenv, query_key=query_key, rettype="acc")
                results = handle.read()
                handle.close()
            except Exception as inst:
                sys.stderr.write("!%s!\n" % ncbi_pjid)
                sys.stderr.write("%s\n" % type(inst))
            nul_f = open(os.devnull, 'w')
            for acc in results.split("\n"):
                try:
                    output = check_output(["blastdbcmd", "-db", "%s/refseq_genomic" % db_dir, "-dbtype", "nucl", "-entry", "%s" % acc, ],
                                universal_newlines=True, stderr=nul_f)
                    fasta_handle.write(output)
                    seq_count = seq_count + 1
                except Exception as inst:
                    pass
            nul_f.close()
        if seq_count:
            fasta_handle.seek(0)
            remove_ambiguous_bases(fasta_handle, "%s/%s.fasta.gz" % (d, ncbi_pjid), seq_format, pat)
            fasta_handle.close()
            os.remove("%s/sequence.fa" % d)
        else:
            fasta_handle.close()
            os.remove("%s/sequence.fa" % d)
    # proceed if retrying or no sequences found from local database
    if retry or not seq_count:
        try:
            fasta_handle = open("%s/sequence.fa" % d, 'w+')
            for start in range(0,count,batch_size):
                handle = Entrez.efetch(db="nuccore", retstart=start, retmax=batch_size, webenv=webenv, query_key=query_key, rettype=seq_format)
                fasta_handle.write(handle.read())
                handle.close()
            fasta_handle.seek(0)
            remove_ambiguous_bases(fasta_handle, "%s/%s.fasta.gz" % (d, ncbi_pjid), seq_format, pat)
            fasta_handle.close()
            os.remove("%s/sequence.fa" % d)
        except urllib2.HTTPError:
            fasta_handle.close()
            os.remove("%s/sequence.fa" % d)
            #if we have a HTTPError, retry
            fasta_file = "%s/%s.fasta.gz" % (d, ncbi_pjid)
            if os.path.isfile(fasta_file):
                os.remove(fasta_file)
            sys.stderr.write("Retrying %s\n" % ncbi_pjid)    
            process_genomes(base, record, seq_format, pat, db_dir, retry+1, max_retry)
            
def main():
    parser = argparse.ArgumentParser(description='Process sequences for KMerge')
    parser.add_argument('mode', metavar='MODE', help='run mode', choices=['download_fasta', 'remove_N'])
    parser.add_argument('i', metavar='INPUT', type=FileType(), default=sys.stdin, help='location of input file')
    parser.add_argument('-d', '--dir', metavar='OUTPUT', default=os.getcwd(), help='location to write output files')
    parser.add_argument('-b', '--db_dir', metavar='DB_DIR', default=os.getcwd(), help='location of refseq genomic database')
    parser.add_argument('-f', '--format', metavar='FORMAT', default='fastq', help='format of input file', choices=['fasta', 'fastq', 'fastq-solexa'])
    parser.add_argument('-t', '--num_threads', metavar='NUM THREADS', default=1, type=int, help='specify number of threads to use')
    parser.add_argument('-r', '--max_retry', metavar='MAX_RETRY', default=2, help='max number of times to retry genome download on http error') 
    args = parser.parse_args()

    mode = args.mode
    pat = re.compile(r'N+', re.IGNORECASE)
    Entrez.email = 'dar326@cornell.edu'
    split_size = 100
    
    
    if mode == 'download_fasta':
        batch_size = 200
        base = args.dir + "/"
        seq_format = "fasta"

        pool = threadpool.ThreadPool(args.num_threads)
        
        nv_pjids = []
        #non-viruses
        try:
            nv_pjids = fetch_non_virus_bp_ids(batch_size)
        except Exception as inst:
            sys.stderr.write("Error converting non-virus assembly ids to bioproject ids\n")
            sys.exit()

        v_pjids = []
        #viruses
        try:
            v_pjids = fetch_virus_bp_ids(batch_size)    
        except Exception as inst:
            sys.stderr.write("Error converting non-virus assembly ids to bioproject ids\n")
            sys.exit()

        
        pjids = list(set(nv_pjids + v_pjids))

        classifications = fetch_classifications(pjids)
        # get the pjids that have classifications and ignore rest
        pjids = classifications.keys()
        
        records = []
        try:
            records = process_seqs.fetch_link_ids(pjids, "bioproject", "nuccore", split_size)        
        except Exception as inst:
            sys.stderr.write("Error linking to sequences in nuccore\n")
            sys.exit()

        for record in records:
            request = threadpool.WorkRequest(process_genomes, args=[base, record, seq_format, pat, args.db_dir, 0, args.max_retry])
            pool.putRequest(request) 
            
        pool.wait()
        
    elif mode == 'remove_N':
        handle = args.i
        seq_format = args.format
        base = args.dir

        remove_ambiguous_bases(handle, base + "/sample.fasta.gz", seq_format, pat)


            
if __name__ == '__main__':
    main()
