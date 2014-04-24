#!/usr/bin/python
import os, argparse
from Bio import SeqIO
from Bio import Entrez
from StringIO import StringIO
import gzip
import pandas as pd

base = '/zenodotus/dat01/mason_lab_scratch/cmlab_scratch008/darryl/k-mer/data/'

parser = argparse.ArgumentParser(description='Download genomes')
parser.add_argument('splitNum', metavar='s', help='split number')
args = parser.parse_args()

splitNum = args.splitNum

df = pd.read_csv("gold.csv", dtype=object))
genomes = df[(df['PROJECT TYPE'] == 'Whole Genome Sequencing') & (df['PROJECT STATUS'] == 'Complete and Published') & (df['AVAILABILITY'] == 'Public') & (df['SEQUENCING STATUS'] == 'Complete')].dropna(subset=['NCBI PROJECT ID'])

Entrez.email = 'dar326@cornell.edu'

for ncbi_pjid in genomes['NCBI PROJECT ID']:
    os.mkdir(base + ncbi_pjid)
    os.chdir(base + ncbi_pjid)
    try:
        row = genomes[genomes['NCBI PROJECT ID'] == ncbi_pjid]
        # write out classifications to taxonomy.txt
        handle = Entrez.elink(dbfrom="bioproject", db="nuccore", id=ncbi_pjid)
        record = Entrez.read(handle)
        links = record[0]['LinkSetDb'][0]['Link']
        records = []
        for link in links:
            ncid = link['Id']
            handle = Entrez.efetch(db="nuccore", id=ncid, rettype="fasta")
            records = records + list(SeqIO.parse(handle, "fasta"))
            # Great place to split records on Ns!!!
            outHandle = StringIO()
            SeqIO.write(records, outHandle, "fasta")
        SeqIO.write(records, open("%s.fasta" % ncid, 'w'), 'fasta')
        f_in = open("%s.fasta" % ncbi_pjid, 'rb')
        f_out = gzip.open("%s.fa.gz" % ncbi_pjid, 'wb')
        # split on Ns here and write to f_in
        f_out.writelines(f_in)
        f_out.close()
        f_in.close()
    except Exception as inst:
        print ncbi_pjid
        print type(inst)
        print inst.args
        print inst
    
    os.chdir(base)

# for Entrez (virus) genomes

#handle = Entrez.esearch(db="nucleotide", term="NC_015932")
#record = Entrez.read(handle)
#ncid=record['IdList'][0]
#handle = Entrez.elink(dbfrom="nuccore", db="bioproject", id=ncid)
#record = Entrez.read(handle)
#ncbi_pjid = record[0]['LinkSetDb'][0]['Link'][0]['Id']
#handle = Entrez.efetch(db="nuccore", id=ncid, rettype="fasta")
#records = []
#records = records.extend(list(SeqIO.parse(handle, "fasta")))
#f_in = open("%s.fasta" % ncbi_pjid, 'rb')
#f_out = gzip.open("%s.fa.gz" % ncbi_pjid, 'wb')
#f_out.writelines(f_in)
#f_out.close()
#f_in.close()


# get taxonomy
#handle = Entrez.elink(dbfrom="bioproject", db="taxonomy", id=ncbi_pjid)
#record = Entrez.read(handle)
#tax_id = record[0]['LinkSetDb'][0]['Link'][0]['Id']
#handle = Entrez.efetch(db="taxonomy", id=tax_id)
#record = Entrez.read(handle)
#tax_file = open("taxonomy.txt", 'wb')
#for taxon in record[0]['LineageEx']:
    #if(taxon['Rank'] != 'no rank'):
        #tax_file.write("%s\t%s\n" % (taxon['Rank'], taxon['ScientificName']))
