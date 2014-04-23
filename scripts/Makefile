# get_genomes.py
# perl ../scripts/contiguous_fasta.pl <(zcat genome.test.fa.gz) | ~/bin/fastx_toolkit/bin/fasta_formatter -w 80 | gzip > genome.test.contig.fa.gz

# for FASTQ

# perl ../scripts/contiguous_fasta.pl <(zcat sample.fastq.gz | ~/fastx_toolkit/bin/fastq_to_fasta -Q33) | ~/bin/fastx_toolkit/bin/fasta_formatter -w 80 | gzip > sample.contig.fa.gz

#1) for dir in $(ls -d */); do org="${dir%?}"; rename .fasta.gz .fa.gz $org/*.fasta.gz; perl ~/Development/kmerge/scripts/contiguous_fasta.pl <(zcat $org/*.fa.gz) | gzip > $org/$org.fasta.gz;  done
#2) rm */*.fa.gz

all: reference.h5

# download gold.xls 
# convert gold.xls to gold.csv
# download genomes
# remove Ns from genome files
# build reference.h5