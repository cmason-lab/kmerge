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
    tax_file.write("%s\t%s\n" % (tax_record[0]['Rank'], tax_record[0]['ScientificName']))
    tax_file.close()

def fetch_classifications(bioproject_ids, batch_size):
    handle = Entrez.elink(dbfrom="bioproject", db="taxonomy", id=bioproject_ids)
    link_results = Entrez.read(handle)
    handle.close()
    tax_ids = [elem['LinkSetDb'][0]['Link'][0]['Id'] for elem in link_results]
    bp_ids = [elem['IdList'][0] for elem in link_results]
    lookup = dict(zip(tax_ids, bp_ids))
    post_handle = Entrez.epost(db='taxonomy', id=",".join(tax_ids), usehistory="y")
    results = Entrez.read(post_handle)
    post_handle.close()
    webenv = results["WebEnv"]
    query_key = results["QueryKey"]
    taxonomy_records = []
    fetch_handle = Entrez.efetch(db='taxonomy', webenv=webenv, query_key=query_key)
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
        
        
        #non-viruses
        nv_pjids = []
        try:
            nv_handle = Entrez.esearch(db="assembly", term='(chromosome[ASLV] OR "Gapless Chromosome"[ASLV] OR "Chromosome with gaps"[ASLV]) AND full-genome-representation[Property] AND assembly_nuccore_refseq[FILT] AND assembly_pubmed[FILT] AND (latest[Property] OR latest_refseq[Property] OR latest_genbank[Property]) NOT suppressed_refseq[Property]', usehistory="y")
            nv_results = Entrez.read(nv_handle)
            nv_webenv = nv_results["WebEnv"]
            nv_query_key = nv_results["QueryKey"]
            count = int(nv_results["Count"])
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
