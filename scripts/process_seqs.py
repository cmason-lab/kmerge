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
    handle.close()
    tax_file = open("%s/taxonomy.txt" % d, 'w')
    for taxon in tax_record[0]['LineageEx']:
        if(taxon['Rank'] != 'no rank'):
            tax_file.write("%s\t%s\n" % (taxon['Rank'], taxon['ScientificName']))

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
        handle = Entrez.efetch(db="nuccore", id=ncid_list, rettype=seq_format)
        remove_ambiguous_bases(handle, "%s/%s.fasta.gz" % (d, ncbi_pjid), seq_format, pat)
        handle.close()
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
    
    elif mode == 'prjnaTOrefseq':
        ids = []
        if not args.no_gold:
            df = pd.read_csv("gold.csv", dtype=object)
            genomes = df[(df['PROJECT TYPE'] == 'Whole Genome Sequencing') & (df['PROJECT STATUS'] == 'Complete and Published') & (df['AVAILABILITY'] == 'Public') & (df['SEQUENCING STATUS'] == 'Complete')].dropna(subset=['NCBI PROJECT ID'])
            ids.extend(list(genomes['NCBI PROJECT ID']))
        if args.i:        
            for line in args.i:
                ids.append(line.strip())

        find_refseq_records(ids)
    elif mode == 'download_fasta':

        base = args.dir + "/"
        seq_format = "fasta"
        

        pool = threadpool.ThreadPool(args.num_threads)
        
        # process genomes
        ids = []
        for line in args.i:
            ids.append(line.strip())
        records = []
        for split_ids in [ ids[i:i+split_size] for i in range(0, len(ids), split_size) ]:
            handle = Entrez.elink(dbfrom="bioproject", db="nuccore", id=split_ids)
            records.extend(Entrez.read(handle))
            handle.close()
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
