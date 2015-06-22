#!/usr/bin/python
import os, argparse, sys, re
from Bio import Entrez
import threadpool
import urllib2
import gzip
from subprocess32 import check_output, CalledProcessError
import shutil
from bs4 import BeautifulSoup
from Bio import SeqIO

def get_taxonomy(d, ncbi_asid, classifications): # will raise KeyError if key not present
    tax_record = classifications[ncbi_asid]
    tax_file = open("%s/taxonomy.txt" % d, 'w')
    for rank, name in tax_record.iteritems():
        tax_file.write("%s\t%s\n" % (rank, name))
    tax_file.close()

def fetch_classifications(assembly_ids, batch_size=20):
    if batch_size < 20:
        batch_size = 20 # minimum size for NCBI records
    lookup = fetch_link_ids(assembly_ids, "assembly", "taxonomy", None, batch_size)
    tax_ids = [tax_id for sublist in lookup.values() for tax_id in sublist]
    as_ids = [as_id for as_id in lookup.keys()]
    # reverse lookup keys and values
    lookup = dict(zip(tax_ids, as_ids))
    post_handle = Entrez.epost(db='taxonomy', id=",".join(tax_ids), usehistory="y")
    results = Entrez.read(post_handle)
    post_handle.close()
    webenv = results["WebEnv"]
    query_key = results["QueryKey"]
    taxonomy_records = []
    for start in range(0, len(tax_ids), batch_size):
        fetch_handle = Entrez.efetch(db='taxonomy', retstart=start, retmax=batch_size, webenv=webenv, query_key=query_key)
        taxonomy_records.extend(Entrez.read(fetch_handle))
    d = dict.fromkeys(as_ids)
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


def fetch_assembly_ids_and_stats(batch_size=20):
    if batch_size < 20:
        batch_size = 20
    assembly_ids = []
    handle = Entrez.esearch(db="assembly", term='("Complete Genome"[ASLV] OR chromosome[ASLV] OR "Gapless Chromosome"[ASLV] OR "Chromosome with gaps"[ASLV]) AND full-genome-representation[Property] AND assembly_nuccore_refseq[FILT] AND (latest[Property] OR latest_refseq[Property] OR latest_genbank[Property]) NOT replaced[Property] NOT suppressed_refseq[Property]', usehistory="y")
    results = Entrez.read(handle)
    handle.close()
    webenv = results["WebEnv"]
    query_key = results["QueryKey"]
    count = int(results["Count"])
    d = {}
    for start in range(0,count,batch_size):
        handle = Entrez.esummary(db="assembly", retstart=start, retmax=batch_size, webenv=webenv, query_key=query_key)
        data = Entrez.read(handle,validate=False)
        handle.close()
        for assembly in data['DocumentSummarySet']['DocumentSummary']:
            assembly_uid = assembly.attributes['uid']
            if not assembly_uid in d:
                d[assembly_uid] = {}
                soup = BeautifulSoup(assembly['Meta'])
                d[assembly_uid]["sequence_length"] = int(soup.find(category="total_length").text)
                d[assembly_uid]["num_sequences"] = int(soup.find(category="replicon_count").text)
    return d


def fetch_link_ids(lookup_ids, from_db, to_db, linkname=None, split_size=20):
    if split_size < 20:
        split_size = 20
    results = []
    for split_ids in [ lookup_ids[i:i+split_size] for i in range(0, len(lookup_ids), split_size) ]:
        if linkname:
            handle = Entrez.elink(dbfrom=from_db, db=to_db, linkname=linkname, id=split_ids)
        else:
            handle = Entrez.elink(dbfrom=from_db, db=to_db, id=split_ids)
        results.extend(Entrez.read(handle))
        handle.close()

    id_map = {}
    for elem in results:
        from_id = elem['IdList'][0]
        if len(elem['LinkSetDb']) == 0: # no taxonomy record
            continue
        to_list = elem['LinkSetDb'][0]['Link']
        to_ids = [ rec['Id'] for rec in to_list]
        id_map[from_id] = to_ids

    return id_map

def get_sequence_from_refseq(file_handle, accession, db_dir):
    try:
        output = check_output(["blastdbcmd", "-db", "%s/refseq_genomic" % db_dir, "-dbtype", "nucl", "-entry", "%s" % accession, ],
                              universal_newlines=True, stderr=open(os.devnull, 'w'))
        file_handle.write(output)
        return True
    except Exception as inst:
        return False

def verify_genome(d, ncbi_asid, stats):
    count = 0
    seq_length = 0
    fasta_handle = gzip.open("%s/%s.fasta.gz" % (ncbi_asid, ncbi_asid), 'rb')
    for record in SeqIO.parse(fasta_handle, "fasta"):
        count = count+1
        seq_length = seq_length + len(record)
    
    if count != stats['num_sequences']:
        print "Sequence count does not match for BioProjectID %s (my count: %d, actual: %d)" % (ncbi_asid, count, stats['num_sequences'])

    if abs(seq_length - stats['sequence_length']) / float(stats['sequence_length']) > 0.1:
        print "Sequence length difference is above 10 pct threshold (%d, %d) for BioProjectID %s" % (seq_length, stats['sequence_length'], ncbi_asid)
        
    fasta_handle.close()

def process_genomes(base, ncbi_asid, stats, classifications, seq_format, db_dir, check_refseq=True, retry=0, max_retry=2):
    d = base + ncbi_asid    
    if retry > max_retry:
        sys.stderr.write("!%s!: Could not obtain sequences\n" % ncbi_asid)
        shutil.rmtree(d)
        return

    if retry > 0:
        sys.stderr.write("Retrying %s\n" % ncbi_asid)
        
    try:
        os.makedirs(d)
    except OSError:
        if not os.path.isdir(d):
            raise

    try:
        handle = Entrez.elink(dbfrom="assembly", db="nucleotide", id=[ncbi_asid], linkname="assembly_nuccore_refseq")
        results = Entrez.read(handle)
        handle.close()
        nuccore_ids = []
        for link in results[0]['LinkSetDb'][0]['Link']:
            nuccore_ids.append(link['Id'])
        # post the ids so that we don't run into any trouble pulling down the sequences
        post_handle = Entrez.epost(db="nuccore", id=",".join(nuccore_ids), usehistory="y")
        results = Entrez.read(post_handle)
        post_handle.close()
        webenv = results["WebEnv"]
        query_key = results["QueryKey"]
        batch_size = 100
        results = ""
        for start in range(0, len(nuccore_ids), batch_size):
            fetch_handle = Entrez.efetch(db='nuccore', retstart=start, retmax=batch_size, webenv=webenv, query_key=query_key, rettype="acc")
            results = results + fetch_handle.read()
            fetch_handle.close()
        get_taxonomy(d, ncbi_asid, classifications)
        fasta_handle = gzip.open("%s/%s.fasta.gz" % (d, ncbi_asid), 'wb')
        seq_count = 0
        if check_refseq:
            for acc in [line for line in results.split("\n") if line != '']:
                if get_sequence_from_refseq(fasta_handle, acc, db_dir):
                    seq_count = seq_count + 1
        check_refseq = False
        # try to download sequences if no sequences found from local database
        if not seq_count:
            for start in range(0, len(nuccore_ids), batch_size):
                handle = Entrez.efetch(db="nuccore", retstart=start, retmax=batch_size, webenv=webenv, query_key=query_key, rettype=seq_format)
                fasta_handle.write(handle.read())
                handle.close()
        fasta_handle.close()
        verify_genome(d, ncbi_asid, stats)
    except (urllib2.HTTPError, AttributeError):
        #if we have a HTTPError or AttributeError, retry
        fasta_file = "%s/%s.fasta.gz" % (d, ncbi_asid)
        if os.path.isfile(fasta_file):
            fasta_handle.close()
            os.remove(fasta_file)
 
        process_genomes(base, ncbi_asid, classifications, seq_format, db_dir, check_refseq, retry+1, max_retry)
    #except Exception as inst:
    #    sys.stderr.write("!%s!:%s\n" % (ncbi_asid, type(inst)))
        
            
def main():
    parser = argparse.ArgumentParser(description='Process sequences for KMerge')
    parser.add_argument('-d', '--dir', metavar='OUTPUT', default=os.getcwd(), help='location to write output files')
    parser.add_argument('-b', '--db_dir', metavar='DB_DIR', default=os.getcwd(), help='location of refseq genomic database')
    parser.add_argument('-f', '--format', metavar='FORMAT', default='fastq', help='format of input file', choices=['fasta', 'fastq', 'fastq-solexa'])
    parser.add_argument('-t', '--num_threads', metavar='NUM THREADS', default=1, type=int, help='specify number of threads to use')
    parser.add_argument('-r', '--max_retry', metavar='MAX_RETRY', default=2, help='max number of times to retry genome download on http error')
    parser.add_argument('-e', '--email', metavar='EMAIL', help='email address to use for Entrez API queries')
    args = parser.parse_args()

    Entrez.email = args.email
    split_size = 100
    
    
    batch_size = 200
    base = args.dir + "/"
    seq_format = "fasta"

    pool = threadpool.ThreadPool(args.num_threads)
        
    asids = []
    try:
        asids = fetch_non_assembly_ids(batch_size)
    except Exception as inst:
        sys.stderr.write("Error converting non-virus assembly ids to bioproject ids\n")
        sys.exit()


        
    asids = list(set(asids))
    
    classifications = fetch_classifications(asids)
    stats = fetch_sequence_stats(asids)
    # get the pjids that have classifications and ignore rest
    asids = [id for id, c in classifications.iteritems() if c != None]

    for ncbi_asid in asids:
        request = threadpool.WorkRequest(process_genomes, args=[base, ncbi_asid, classifications, seq_format, args.db_dir, True, 0, args.max_retry])
        pool.putRequest(request) 
            
    pool.wait()

            
if __name__ == '__main__':
    main()
