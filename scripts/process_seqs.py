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

def process_gold_genomes(base, genomes, ncbi_pjid, seq_format, pat):
    d = base + ncbi_pjid
    os.mkdir(d)
    try:
        row = genomes[genomes['NCBI PROJECT ID'] == ncbi_pjid]
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
        print ncbi_pjid
        print type(inst)
        print inst.args
        print inst

    os.chdir(base)

def process_entrez_genomes(base, accession_id, seq_format, pat):
    try:
        handle = Entrez.esearch(db="nucleotide", term=accession_id)
        search_record = Entrez.read(handle)
        ncid = search_record['IdList'][0]
        handle = Entrez.elink(dbfrom="nuccore", db="bioproject", id=ncid)
        link_record = Entrez.read(handle)
        ncbi_pjid = link_record[0]['LinkSetDb'][0]['Link'][0]['Id']
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

    os.chdir(base)
            
def main():
    parser = argparse.ArgumentParser(description='Download genomes')
    parser.add_argument('mode', metavar='MODE', help='run mode', choices=['download', 'remove_N'])
    parser.add_argument('-i', '--in', metavar='INPUT', type=FileType(), default=sys.stdin, help='location of input file')
    parser.add_argument('-d', '--dir', metavar='OUTPUT', default=os.getcwd(), help='location to write output files')
    parser.add_argument('-f', '--format', metavar='FORMAT', default='fastq', help='format of input file', choices=['fasta', 'fastq', 'fastq-solexa'])
    parser.add_argument('-t', '--num_threads', metavar='NUM THREADS', default=1, type=int, help='specify number of threads to use')
    args = parser.parse_args()

    mode = args.mode
    pat = re.compile(r'N+', re.IGNORECASE)
    
    if mode == 'download':
        Entrez.email = 'dar326@cornell.edu'

        base = args.dir + "/"
        seq_format = "fasta"
        
        df = pd.read_csv("gold.csv", dtype=object)
        genomes = df[(df['PROJECT TYPE'] == 'Whole Genome Sequencing') & (df['PROJECT STATUS'] == 'Complete and Published') & (df['AVAILABILITY'] == 'Public') & (df['SEQUENCING STATUS'] == 'Complete')].dropna(subset=['NCBI PROJECT ID'])

        pool = threadpool.ThreadPool(args.num_threads)
        
        # process GOLD genomes
        for ncbi_pjid in genomes['NCBI PROJECT ID']:
            request = threadpool.WorkRequest(process_gold_genomes, args=[base, genomes, ncbi_pjid, seq_format, pat])
            pool.putRequest(request)
            #process_gold_genomes(base, genomes, ncbi_pjid, seq_format, pat)

        #process Entrez genomes
        for accession_id in [ line.strip() for line in open("entrez.dsdna.txt").readlines() + open("entrez.ssdna.txt").readlines()]:
            request = threadpool.WorkRequest(process_entrez_genomes, args=[base, accession_id, seq_format, pat])
            pool.putRequest(request)
            #process_entrez_genomes(base, accession_id, seq_format, pat)

        pool.wait()
    elif mode == 'remove_N':
        handle = args.i
        seq_format = args.format
        base = args.dir

        remove_ambiguous_bases(handle, base + "/sample.fasta.gz", seq_format, pat)


            
if __name__ == '__main__':
    main()
