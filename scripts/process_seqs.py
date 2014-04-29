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
    tax_file = open("%s/taxonomy.txt" % d, 'w')
    for taxon in tax_record[0]['LineageEx']:
        if(taxon['Rank'] != 'no rank'):
            tax_file.write("%s\t%s\n" % (taxon['Rank'], taxon['ScientificName']))

def find_refseq_records(ncbi_pjid_list):
    handle = Entrez.elink(dbfrom="bioproject", db="bioproject", id=ncbi_pjid_list)
    records = Entrez.read(handle)
    query_ids = set()
    skip_ids = set()
    # loop over each record of project ids and linked project ids
    for record in records:
        base_id = record['IdList'][0]
        if base_id in query_ids or base_id in skip_ids:
            continue
        handle = Entrez.esummary(db="bioproject", id=base_id)
        summary_record = Entrez.read(handle)
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
            handle = Entrez.esummary(db="bioproject", id=id)
            summary_record = Entrez.read(handle)
            if summary_record['DocumentSummarySet']['DocumentSummary'][0]['Project_Data_Type'] == "RefSeq Genome":
                query_ids.add(id)
                skip_rest = True
        # if no "RefSeq Genome" record found, print organism name for manual search
        if not skip_rest: # didn't find refseq genome
            sys.stderr.write("Further Review: %s\n" % (summary_record['DocumentSummarySet']['DocumentSummary'][0]['Project_Name']))
    # print pjid records for genomes
    for query_id in query_ids:
        sys.stdout.write("%s\n" % query_id)
            
#def process_gold_genomes(base, genomes, ncbi_pjid, seq_format, pat):
def process_genomes(base, ncbi_pjid, seq_format, pat):
    d = base + ncbi_pjid
    os.mkdir(d)
    try:
        #row = genomes[genomes['NCBI PROJECT ID'] == ncbi_pjid]
        # write out classifications to taxonomy.txt
        get_taxonomy(d, ncbi_pjid)
        handle = Entrez.elink(dbfrom="bioproject", db="nuccore", id=ncbi_pjid)
        record = Entrez.read(handle)
        links = record[0]['LinkSetDb'][0]['Link']
        records = []
        for link in links:
            ncid = link['Id']
            handle = Entrez.efetch(db="nuccore", id=ncid, rettype=seq_format)
            remove_ambiguous_bases(handle, "%s/%s.fasta.gz" % (d, ncbi_pjid), seq_format, pat)
    except Exception as inst:
        print "!%s!" % ncbi_pjid
        print type(inst)
        print inst.args
        print inst

def process_entrez_genomes(base, accession_id, seq_format, pat):
    try:
        handle = Entrez.esearch(db="nucleotide", term=accession_id)
        search_record = Entrez.read(handle)
        ncid = search_record['IdList'][0]
        handle = Entrez.elink(dbfrom="nuccore", db="bioproject", id=ncid)
        link_record = Entrez.read(handle)
        ncbi_pjid = link_record[0]['LinkSetDb'][0]['Link'][0]['Id']
        print "%s <-> %s" % (accession_id, ncbi_pjid)
        d = base + ncbi_pjid
        os.mkdir(d)
        handle = Entrez.efetch(db="nuccore", id=ncid, rettype=seq_format)
        remove_ambiguous_bases(handle, "%s/%s.fasta.gz" % (d, ncbi_pjid), seq_format, pat)
        get_taxonomy(d, ncbi_pjid)
    except Exception as inst:
        print accession_id
        print type(inst)
        print inst.args
        print inst
            
def main():
    parser = argparse.ArgumentParser(description='Download genomes')
    parser.add_argument('mode', metavar='MODE', help='run mode', choices=['download_fasta', 'remove_N', 'virusTOprjna', 'prjnaTOrefseq'])
    parser.add_argument('-i', '--in', metavar='INPUT', type=FileType(), default=sys.stdin, help='location of input file')
    parser.add_argument('-d', '--dir', metavar='OUTPUT', default=os.getcwd(), help='location to write output files')
    parser.add_argument('-f', '--format', metavar='FORMAT', default='fastq', help='format of input file', choices=['fasta', 'fastq', 'fastq-solexa'])
    parser.add_argument('-t', '--num_threads', metavar='NUM THREADS', default=1, type=int, help='specify number of threads to use')
    args = parser.parse_args()

    mode = args.mode
    pat = re.compile(r'N+', re.IGNORECASE)
    Entrez.email = 'dar326@cornell.edu'

    
    if mode == 'virusTOprjna':
        pass
    elif mode == 'prjnaTOrefseq':
        ids = []
        for line in open("prjna_list.txt", "r"):
            ids.append(line.strip())

        find_refseq_records(ids)
    elif mode == 'download_fasta':
        Entrez.email = 'dar326@cornell.edu'

        base = args.dir + "/"
        seq_format = "fasta"
        
        df = pd.read_csv("gold.csv", dtype=object)
        genomes = df[(df['PROJECT TYPE'] == 'Whole Genome Sequencing') & (df['PROJECT STATUS'] == 'Complete and Published') & (df['AVAILABILITY'] == 'Public') & (df['SEQUENCING STATUS'] == 'Complete')].dropna(subset=['NCBI PROJECT ID'])

        # get virus NCBI project ids (refactor as virusTOprjna
        virus_ncbi_pjids = set()
        for accession_id in [ line.strip() for line in open("entrez.dsdna.txt").readlines() + open("entrez.ssdna.txt").readlines()]:
            try:
                handle = Entrez.esearch(db="nucleotide", term=accession_id)
                search_record = Entrez.read(handle)
                ncid = search_record['IdList'][0]
                handle = Entrez.elink(dbfrom="nuccore", db="bioproject", id=ncid)
                link_record = Entrez.read(handle)
                ncbi_pjid = link_record[0]['LinkSetDb'][0]['Link'][0]['Id']
                virus_ncbi_pjids.add(ncbi_pjid)
                print "!%s , %s!" % (accession_id, ncbi_pjid)
            except Exception as inst:
                print "!%s!" % accession_id
                print type(inst)
                print inst.args
                print inst

        pool = threadpool.ThreadPool(args.num_threads)
        
        # process genomes
        for ncbi_pjid in list(genomes['NCBI PROJECT ID']) + list(virus_ncbi_pjids):
            #request = threadpool.WorkRequest(process_gold_genomes, args=[base, genomes, ncbi_pjid, seq_format, pat])
            request = threadpool.WorkRequest(process_genomes, args=[base, ncbi_pjid, seq_format, pat])
            pool.putRequest(request)
            #process_gold_genomes(base, genomes, ncbi_pjid, seq_format, pat)

        #process Entrez genomes
        #for accession_id in [ line.strip() for line in open("entrez.dsdna.txt").readlines() + open("entrez.ssdna.txt").readlines()]:
            #request = threadpool.WorkRequest(process_entrez_genomes, args=[base, accession_id, seq_format, pat])
            #pool.putRequest(request)
            #process_entrez_genomes(base, accession_id, seq_format, pat)

        pool.wait()
        
    elif mode == 'remove_N':
        handle = args.i
        seq_format = args.format
        base = args.dir

        remove_ambiguous_bases(handle, base + "/sample.fasta.gz", seq_format, pat)


            
if __name__ == '__main__':
    main()
