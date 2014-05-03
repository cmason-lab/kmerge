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
import gzip

class FileType(argparse.FileType):
    def __call__(self, string):
        if string.endswith('.gz'):
            return bgzf.open(string, 'ab')
        return super(FileType, self).__call__(string)

def download_genome_sequences(base, asm_id):
    summary = None
    ncbi_pjid = None
    try:
        handle = Entrez.esummary(db='assembly', id=asm_id)
        summary = Entrez.read(handle)
        handle.close()
        ncbi_pjid = summary['DocumentSummarySet']['DocumentSummary'][0]['RS_BioProjects'][0]['BioprojectId']
    except Exception as inst:
        sys.stderr.write("Summary data for linked assembly %s to bioproject %s inaccessible\n" % (asm_id, ncbi_pjid))
        return

    os.mkdir(base + ncbi_pjid)
    get_taxonomy(base + ncbi_pjid, ncbi_pjid)
    
    soup = None
    try:
        soup = BeautifulSoup(summary['DocumentSummarySet']['DocumentSummary'][0]['Meta'])
    except Exception as inst:
        sys.stderr.write("Meta data for linked assembly %s to bioproject %s inacessible\n" % (asm_id, ncbi_pjid))
        return

    chr_count = int(soup.find(category="chromosome_count").text)
    seq_count = int(soup.find(category="replicon_count").text)
    length = int(soup.find(category="total_length").text)

    ftp_url = None
    try:
        ftp_url = soup.ftpsites.find(type="GenBank").text
        print ftp_url
        url_split = urlparse(ftp_url)
        ftp = FTP(url_split.netloc)
        ftp.login()
        ftp.cwd(url_split.path + "Primary_Assembly/assembled_chromosomes/FASTA")
        files = ftp.nlst("*.fa.gz")
        assert len(files) == chr_count
        gzf = gzip.GzipFile("%s/%s.full.fa.gz" % (base + ncbi_pjid, ncbi_pjid), 'wb')
        sio = StringIO()        
        for filename in files:
            try:
                # writing to StringIO object
                ftp.retrbinary('RETR ' + filename, sio.write)
            except Exception as inst:
                sys.stderr.write("Unable to download %s from linked assembly %s to bioproject %s inacessible\n" % (filename, asm_id, ncbi_pjid))
                return
            f_handle.close()
        # uncompress data
        zippy = gzip.GzipFile(fileobj=sio)
        gzf.write(zippy.read())
        sio.close()
        zippy.close()
        gzf.close()
    except Exception as inst:
        sys.stderr.write("GenBank url for linked assembly %s to bioproject %s inacessible\n" % (asm_id, ncbi_pjid))
        return
    
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
    try:
        handle = Entrez.elink(dbfrom="bioproject", db="taxonomy", id=ncbi_pjid)
        link_record = Entrez.read(handle)
        tax_id = link_record[0]['LinkSetDb'][0]['Link'][0]['Id']
        handle = Entrez.efetch(db="taxonomy", id=tax_id)
        tax_record = Entrez.read(handle)
        handle.close()
        tax_file = open("%s/taxonomy.txt" % d, 'w')
        for taxon in tax_record[0]['LineageEx']:
            if(taxon['Rank'] != 'no rank'):
                tax_file.write("%s\t%s\n" % (taxon['Rank'], taxon['ScientificName']))
    except Exception as inst:
        sys.stderr.write("Error getting taxonomy information for %s\n" % ncbi_pjid)

def convert_virus_accessions(ids, split_size=100):
    # get virus NCBI project ids
    virus_ncbi_pjids = set()
    records = []
    for split_ids in [ ids[i:i+split_size] for i in range(0, len(ids), split_size) ]:
        handle = Entrez.elink(dbfrom="nuccore", db="bioproject", id=split_ids)
        records.extend(Entrez.read(handle))
        handle.close()
    for record in records:
        try:
            ncbi_pjid = record['LinkSetDb'][0]['Link'][0]['Id']
            virus_ncbi_pjids.add(ncbi_pjid)
        except Exception as inst:
            sys.stderr.write("!%s!\n" % record['IdList'][0])
        
    for ncbi_pjid in virus_ncbi_pjids:
        sys.stdout.write("%s\n" % ncbi_pjid)
                                                                                                                                                                                                                                                                
            
def find_refseq_records(ncbi_pjid_list, split_size=100):
    records = []
    for split_ids in [ ncbi_pjid_list[i:i+split_size] for i in range(0, len(ncbi_pjid_list), split_size) ]:
        handle = Entrez.elink(dbfrom="bioproject", db="bioproject", id=split_ids)
        records.extend(Entrez.read(handle))
        handle.close()
    query_ids = set()
    skip_ids = set()
    # loop over each record of project ids and linked project ids
    for record in records:
        base_id = record['IdList'][0]
        if base_id in query_ids or base_id in skip_ids:
            continue
        summary_record = None
        try:
            handle = Entrez.esummary(db="bioproject", id=base_id)
            summary_record = Entrez.read(handle)
            handle.close()
        except:
            sys.stderr.write("!%s!\n" % base_id)
        if summary_record['DocumentSummarySet']['DocumentSummary'][0]['Project_Data_Type'] == "RefSeq Genome":
            # add id to query set
            query_ids.add(base_id)
            for other in record['LinkSetDb']:
                # add related ids to skip set
                skip_ids.add(other['Link'][0]['Id'])
            continue
        skip_rest = False
        for other in record['LinkSetDb']:
            id = other['Link'][0]['Id']
            if id in query_ids or id in skip_ids:
                continue
            if skip_rest:
                skip_ids.add(id)
                continue
            # get summary
            summary_record = None
            try:
                handle = Entrez.esummary(db="bioproject", id=id)
                summary_record = Entrez.read(handle)
                handle.close()
            except:
                sys.stderr.write("!%s!\n" % id)
            if summary_record['DocumentSummarySet']['DocumentSummary'][0]['Project_Data_Type'] == "RefSeq Genome":
                query_ids.add(id)
                skip_rest = True
        # if no "RefSeq Genome" record found, print organism name for manual search
        if not skip_rest: # didn't find refseq genome
            sys.stderr.write("Further Review: %s\n" % (summary_record['DocumentSummarySet']['DocumentSummary'][0]['Project_Name']))
    # print pjid records for genomes
    for query_id in query_ids:
        sys.stdout.write("%s\n" % query_id)
            
def process_genomes(base, record, seq_format, pat):
    ncbi_pjid = record['IdList'][0]
    d = base + ncbi_pjid
    os.mkdir(d)
    try:
        # write out classifications to taxonomy.txt
        get_taxonomy(d, ncbi_pjid)
    
        ncid_list = ",".join([ link['Id'] for link in record['LinkSetDb'][0]['Link'] ])
        handle = Entrez.esearch(db="nucleotide", term='%s[BioProject] AND biomol_genomic[PROP] AND srcdb_refseq[PROP] NOT ALTERNATE_LOCUS[Keyword] NOT FIX_PATCH[Keyword] NOT NOVEL_PATCH[Keyword] NOT CONTIG[Title] NOT SCAFFOLD[Title]' % ncbi_pjid, usehistory='y')
        results = Entrez.read(handle)
        webenv = results["WebEnv"]
        query_key = results["QueryKey"]
        count = int(results["Count"])
        batch_size = 100
        #gzf = gzip.GzipFile("%s/%s.fasta.gz" % (d, ncbi_pjid), 'wb')
        for start in range(0,count,batch_size):
            handle = Entrez.efetch(db="nuccore", retstart=start, retmax=batch_size, webenv=webenv, query_key=query_key, rettype=seq_format)
            remove_ambiguous_bases(handle, "%s/%s.fasta.gz" % (d, ncbi_pjid), seq_format, pat)
            #gzf.write(handle.read())
            handle.close()
        #gzf.close()
    except Exception as inst:
        sys.stderr.write("!%s!\n" % ncbi_pjid)
        sys.stderr.write("%s\n" % type(inst))
            
def main():
    parser = argparse.ArgumentParser(description='Process sequences for KMerge')
    parser.add_argument('mode', metavar='MODE', help='run mode', choices=['download_fasta', 'remove_N', 'virusTOprjna', 'prjnaTOrefseq'])
    parser.add_argument('i', metavar='INPUT', type=FileType(), default=sys.stdin, help='location of input file')
    parser.add_argument('-d', '--dir', metavar='OUTPUT', default=os.getcwd(), help='location to write output files')
    parser.add_argument('-f', '--format', metavar='FORMAT', default='fastq', help='format of input file', choices=['fasta', 'fastq', 'fastq-solexa'])
    parser.add_argument('-t', '--num_threads', metavar='NUM THREADS', default=1, type=int, help='specify number of threads to use')
    parser.add_argument('-g', '--no_gold', default=False, action='store_true', help='do not use GOLD input file ("gold.csv")')
    args = parser.parse_args()

    mode = args.mode
    pat = re.compile(r'N+', re.IGNORECASE)
    Entrez.email = 'dar326@cornell.edu'
    split_size = 100
    
    if mode == 'virusTOprjna':
        ids = []
        for line in args.i:
            ids.append(line.strip())
            
        convert_virus_accessions(ids)
    elif mode == 'download_fasta':
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
        pjids = nv_pjids + v_pjids
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
            request = threadpool.WorkRequest(process_genomes, args=[base, record, seq_format, pat])
            pool.putRequest(request) 
            
        pool.wait()
        
    elif mode == 'remove_N':
        handle = args.i
        seq_format = args.format
        base = args.dir

        remove_ambiguous_bases(handle, base + "/sample.fasta.gz", seq_format, pat)


            
if __name__ == '__main__':
    main()
