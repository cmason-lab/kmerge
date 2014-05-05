#!/usr/bin/python
import os, argparse, sys, re
from Bio import SeqIO
from Bio import Entrez
from StringIO import StringIO
import pandas as pd
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import DNAAlphabet
from Bio.Seq import Seq
from Bio import bgzf
import threadpool
from bs4 import BeautifulSoup
from ftplib import FTP
from urlparse import urlparse
import urllib2
import gzip
import time
from subprocess32 import check_output, CalledProcessError
import shutil

class FileType(argparse.FileType):
    def __call__(self, string):
        if string.endswith('.gz'):
            return bgzf.open(string, 'ab')
        return super(FileType, self).__call__(string)

    
def remove_ambiguous_bases(file_handle, output_filename, seq_format, pat):
    out_file = bgzf.open(output_filename, 'ab')
    output = StringIO()
    for record in SeqIO.parse(file_handle, seq_format):
        sequences = re.split(pat, str(record.seq))
        seq_num = 0
        for sequence in sequences:
            new_rec = SeqRecord(Seq(sequence, DNAAlphabet), description="", id="%s" % (record.id + ("" if len(sequences) == 1 else "/%d" % (seq_num))))
            seq_num = seq_num+1
            SeqIO.write(new_rec, output, "fasta")
    out_file.write(output.getvalue())
    output.close()
    out_file.close()

def get_taxonomy(d, ncbi_pjid):
    handle = Entrez.elink(dbfrom="bioproject", db="taxonomy", id=ncbi_pjid)
    link_record = Entrez.read(handle)
    tax_id = link_record[0]['LinkSetDb'][0]['Link'][0]['Id']
    handle = Entrez.efetch(db="taxonomy", id=tax_id)
    tax_record = Entrez.read(handle)
    handle.close()
    tax_file = open("%s/taxonomy.txt" % d, 'w')
    for taxon in tax_record[0]['LineageEx']:
        if('Rank' in taxon and taxon['Rank'] != 'no rank'):
            tax_file.write("%s\t%s\n" % (taxon['Rank'], taxon['ScientificName']))
                                                                                                                                                                                                                                                                
            
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
    if not retry:
        # write out classifications to taxonomy file
        try:
            get_taxonomy(d, ncbi_pjid)
        except Exception as inst:
            sys.stderr.write("!%s!: Error getting taxonomy information for %s\n" % ncbi_pjid)
            # delet directory if can't get taxonomy
            shutil.rmtree(d)
            return
        results = None
        try:
            handle = Entrez.esearch(db="nucleotide", term='%s[BioProject] AND biomol_genomic[PROP] AND srcdb_refseq[PROP] NOT ALTERNATE_LOCUS[Keyword] NOT FIX_PATCH[Keyword] NOT NOVEL_PATCH[Keyword] NOT CONTIG[Title] NOT SCAFFOLD[Title]' % ncbi_pjid, usehistory='y')
            results = Entrez.read(handle)
        except Exception as inst:
            sys.stderr.write("!%s!\n" % ncbi_pjid)
            sys.stderr.write("%s\n" % type(inst))
        webenv = results["WebEnv"]
        query_key = results["QueryKey"]
        count = int(results["Count"])
        batch_size = 100
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
            for start in range(0,count,batch_size):
                handle = Entrez.efetch(db="nuccore", retstart=start, retmax=batch_size, webenv=webenv, query_key=query_key, rettype=seq_format)
                remove_ambiguous_bases(handle, "%s/%s.fasta.gz" % (d, ncbi_pjid), seq_format, pat)
                handle.close()
        except urllib2.HTTPError:
            #if we have a HTTPError, retry
            fasta_file = "%s/%s.fasta.gz" % (d, ncbi_pjid)
            if os.path.isfile(fasta_file):
                os.remove(fasta_file)
            time.sleep(2) # sleep for 2 seconds to let server clear up any issues
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
    parser.add_argument('-g', '--no_gold', default=False, action='store_true', help='do not use GOLD input file ("gold.csv")')
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
        
        
        #non-viruses
        nv_pjids = []
        try:
            nv_handle = Entrez.esearch(db="assembly", term='(chromosome[ASLV] OR "Gapless Chromosome"[ASLV]) AND full-genome-representation[Property] AND assembly_nuccore_refseq[FILT] AND assembly_pubmed[FILT] AND (latest[Property] OR latest_refseq[Property] OR latest_genbank[Property])', usehistory="y")
            nv_results = Entrez.read(nv_handle)
            nv_webenv = nv_results["WebEnv"]
            nv_query_key = nv_results["QueryKey"]
            count = int(nv_results["Count"])
            asm_ids = []
            for start in range(0,count,batch_size):
                fetch_handle = Entrez.esummary(db="assembly", retstart=start, retmax=batch_size, webenv=nv_webenv, query_key=nv_query_key)
                nv_summaries = Entrez.read(fetch_handle)                    
                nv_pjids.extend([ summary['RS_BioProjects'][0]['BioprojectId'] for summary in nv_summaries['DocumentSummarySet']['DocumentSummary']])
        except Exception as inst:
            sys.stderr.write("Error converting non-virus assembly ids to bioproject ids\n")
            sys.exit()

            
        #viruses
        v_handle = Entrez.esearch(db="genome", term='(txid29258[Organism:exp] OR txid35237[Organism:exp]) AND complete[Status] AND "RefSeq" AND genome_pubmed[FILT]', usehistory="y")
        v_results = Entrez.read(v_handle)
        v_webenv = v_results["WebEnv"]
        v_query_key = v_results["QueryKey"]
        count = int(v_results["Count"])
        v_summaries = []
        for start in range(0,count,batch_size):
            fetch_handle = Entrez.esummary(db="genome", retstart=start, retmax=batch_size, webenv=v_webenv, query_key=v_query_key)
            v_summaries.extend(Entrez.read(fetch_handle))

        v_pjids = [ summary['ProjectID'] for summary in v_summaries]

        
        pjids = list(set(nv_pjids + v_pjids))

        records = []
        for split_ids in [ pjids[i:i+split_size] for i in range(0, len(pjids), split_size) ]:
            try:
                handle = Entrez.elink(dbfrom="bioproject", db="nuccore", id=split_ids)
                records.extend(Entrez.read(handle))
                handle.close()
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
