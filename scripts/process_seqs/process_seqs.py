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
    lookup = fetch_link_ids(assembly_ids, "assembly", "taxonomy", "assembly_taxonomy", batch_size)
    tax_ids = [tax_id for sublist in lookup.values() for tax_id in sublist]
    as_ids = [as_id for as_id in lookup.keys()]
    # match assembly ids and taxonomy ids
    lookup = dict(zip(as_ids, tax_ids))
    post_handle = Entrez.epost(db='taxonomy', id=",".join(tax_ids), usehistory="y")
    results = Entrez.read(post_handle)
    post_handle.close()
    webenv = results["WebEnv"]
    query_key = results["QueryKey"]
    taxonomy_records = []
    for start in range(0, len(tax_ids), batch_size):
        fetch_handle = Entrez.efetch(db='taxonomy', retstart=start, retmax=batch_size, webenv=webenv, query_key=query_key)
        taxonomy_records.extend(Entrez.read(fetch_handle))
    taxa = {}
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
        taxa[record['TaxId']] = bp_taxonomy
    d = {}
    for as_id, taxa_id in lookup.iteritems():
        d[as_id] = taxa[taxa_id]
    return d


def fetch_assembly_ids(batch_size=20):
    if batch_size < 20:
        batch_size = 20
    assembly_ids = []
    handle = Entrez.esearch(db="assembly", term='("Complete Genome"[ASLV] OR "Chromosome"[ASLV] OR "Gapless Chromosome"[ASLV] OR "Chromosome with gaps"[ASLV]) AND full-genome-representation[Property] AND assembly_nuccore_refseq[FILT] AND (latest[Property] OR latest_refseq[Property] OR latest_genbank[Property]) NOT replaced[Property] NOT suppressed_refseq[Property] AND (txid2[Organism:exp] OR txid2759[Organism:exp] OR txid2157[Organism:exp] OR txid35237[Organism:exp] OR txid29258[Organism:exp])', usehistory="y")
    results = Entrez.read(handle)
    handle.close()
    webenv = results["WebEnv"]
    query_key = results["QueryKey"]
    count = int(results["Count"])
    data = []
    for start in range(0,count,batch_size):
        handle = Entrez.efetch(db="assembly", retstart=start, retmax=batch_size, webenv=webenv, query_key=query_key)
        data.extend(Entrez.read(handle,validate=False))
        handle.close()

    return data


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

def verify_genome(d, ncbi_asid, num_sequences):
    count = 0
    seq_length = 0
    fasta_handle = gzip.open("%s/%s.fasta.gz" % (ncbi_asid, ncbi_asid), 'rb')
    for record in SeqIO.parse(fasta_handle, "fasta"):
        count = count+1
    
    if count != num_sequences:
        sys.stderr.write("!%s!: Sequence count does not match for Assembly ID %s (my count: %d, actual: %d)\n" % (ncbi_asid, ncbi_asid, count, num_sequences))
        
    fasta_handle.close()

def process_genomes(base, ncbi_asid, classifications, seq_format, db_dir, check_refseq=True, retry=0, max_retry=2):
    d = base + ncbi_asid    
    if retry > max_retry:
        sys.stderr.write("!%s!: Could not obtain sequences\n" % ncbi_asid)
        shutil.rmtree(d)
        return

    if retry > 0:
        sys.stdout.write("Retrying %s\n" % ncbi_asid)
    else:
        sys.stdout.write("Processing %s\n" % ncbi_asid)    
        os.makedirs(d)

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
        seqs_not_found = []
        if check_refseq:
            for acc in [line for line in results.split("\n") if line != '']:
                if get_sequence_from_refseq(fasta_handle, acc, db_dir):
                    seq_count = seq_count + 1
        check_refseq = False
        # try to download sequences if not all sequences found from local database
        if seq_count != len(nuccore_ids):
            # reduce batch size
            batch_size = 10
            # move back to the beginning of the file for writing
            fasta_handle.seek(0)
            sys.stdout.write("Attempting to download sequences for %s\n" % ncbi_asid)
            for start in range(0, len(nuccore_ids), batch_size):
                handle = Entrez.efetch(db="nuccore", retstart=start, retmax=batch_size, webenv=webenv, query_key=query_key, rettype=seq_format)
                fasta_handle.write(handle.read())
                handle.close()
            sys.stdout.write("Finished downloading sequences for %s\n" % ncbi_asid)
        fasta_handle.close()
        verify_genome(d, ncbi_asid, len(nuccore_ids))
    except (urllib2.HTTPError, AttributeError, IndexError) as e:
        exc_type = sys.exc_info()[0]
        exc_obj = sys.exc_info()[1]
        exc_tb = sys.exc_info()[2]
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        sys.stdout.write("Retrying %s due to error: %s, %s, %s" % (ncbi_asid, exc_type, fname, exc_tb.tb_lineno))
        #if we have a HTTPError or AttributeError or IndexError, retry
        fasta_file = "%s/%s.fasta.gz" % (d, ncbi_asid)
        if os.path.isfile(fasta_file):
            fasta_handle.close()
            os.remove(fasta_file)
 
        process_genomes(base, ncbi_asid, classifications, seq_format, db_dir, check_refseq, retry+1, max_retry)
    except Exception as inst:
        exc_type = sys.exc_info()[0]
        exc_obj = sys.exc_info()[1]
        exc_tb = sys.exc_info()[2]
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
                        
        sys.stderr.write("!%s!: %s, %s, %s\n" % (ncbi_asid, exc_type, fname, exc_tb.tb_lineno))
        
            
def main():
    parser = argparse.ArgumentParser(description='Process sequences for KMerge')
    parser.add_argument('-d', '--dir', metavar='OUTPUT', default=os.getcwd(), help='location to write output files')
    parser.add_argument('-b', '--db_dir', metavar='DB_DIR', default=os.getcwd(), help='location of refseq genomic database')
    parser.add_argument('-f', '--format', metavar='FORMAT', default='fasta', help='format of input file', choices=['fasta', 'fastq', 'fastq-solexa'])
    parser.add_argument('-t', '--num_threads', metavar='NUM THREADS', default=1, type=int, help='specify number of threads to use')
    parser.add_argument('-r', '--max_retry', metavar='MAX_RETRY', default=2, help='max number of times to retry genome download on http error')
    parser.add_argument('-e', '--email', metavar='EMAIL', help='email address to use for Entrez API queries')
    parser.add_argument('-i', '--input_file', metavar='INPUT_FILE', default=None, help='location of file of NCBI assembly ids to process')
    args = parser.parse_args()

    Entrez.email = args.email
    split_size = 100
    
    
    batch_size = 200
    base = args.dir + "/"
    seq_format = "fasta"

    pool = threadpool.ThreadPool(args.num_threads)

    ids = []
    try:
        if not args.input_file:
            ids = fetch_assembly_ids(batch_size)
        else:
            for line in open(args.input_file).readlines():
                ids.append(line.strip())
    except Exception as inst:
        sys.stderr.write("Error obtaining assembly ids\n")
        

    classifications = fetch_classifications(ids)
    # get the asids that have classifications and ignore rest
    asids = classifications.keys()

    sys.stdout.write("%s genomes will be processed\n" % len(asids))
    
    for ncbi_asid in asids:
        request = threadpool.WorkRequest(process_genomes, args=[base, ncbi_asid, classifications, seq_format, args.db_dir, True, 0, args.max_retry])
        pool.putRequest(request)
                                                                                                                     
            
    pool.wait()

            
if __name__ == '__main__':
    main()
