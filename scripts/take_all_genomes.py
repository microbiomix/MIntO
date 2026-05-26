"""
The script accept in input:
    - cluster file from vamb
    - contigs file
Ouput: 
    the ouput folder
Parameters:
    Minimum length of the contigs (default 500000)
"""

import os
from os import listdir
from os.path import isfile, join
import warnings
import argparse as ap
import sys
import numpy as np
from itertools import groupby
import math
import pandas as pd
import time
import yaml

# Open a file, possibly gzipped

def open_file(filename, mode, compresslevel=2):
    import gzip

    if filename.endswith('.gz'):
        return gzip.open(filename, mode, compresslevel=compresslevel)
    else:
        return open(filename, mode)

def fasta_iter(infile):
    """
    Source: https://www.biostars.org/p/710/
    modified from Brent Pedersen
    Correct Way To Parse A Fasta File In Python
    given a fasta file. yield tuples of header, sequence
    """

    from itertools import groupby

    #first open the file outside
    fh = open_file(infile, mode='rt')

    # ditch the boolean (x[0]) and just keep the header or sequence since
    # we know they alternate.
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        # drop the ">"
        headerStr = header.__next__()[1:].strip()
        # join all sequence lines to one.
        seq = "".join(s.strip() for s in faiter.__next__())
        yield (headerStr, seq)

def read_params():
	p = ap.ArgumentParser( description = ( "Given a cluster table from vamb, it retrieves the genome considereing the minimum length"),
		formatter_class = ap.ArgumentDefaultsHelpFormatter )


	p.add_argument("--vamb_cluster_tsv", 
		help = "the table with the cluster obtained by vamb",
		)
	
	p.add_argument("--binsplit_char",
		help = "String used to concatenate sample-id with contig-id to make unique fasta headers",
		default = "_NODE"
		)

	p.add_argument("--binning_method_name", 
		help = "the binning method used (aaey, aaez, vae512, vaevae512, graphmb",
		default= "")

	p.add_argument("--contigs_file", 
		help = "File containing all the contigs file to parse")

	p.add_argument("--output_folder",
		help = "where to copy the genomes ")

	p.add_argument("--min_fasta_length", 
		help = "Minimum threshold for the fasta length", 
		default = 500000 )

	p.add_argument("--discarded_genomes_info", 
		help = "Path where to put info about discarded genomes")
		
	return(p.parse_args())

args = read_params()

cluster_tsv =  args.vamb_cluster_tsv
contigs_file = args.contigs_file
binning_method = args.binning_method_name
min_fasta = int(args.min_fasta_length)
output = args.output_folder
discarded_genomes = args.discarded_genomes_info


# create a dictionary with contigs and length

dictionary_of_contigs_length = {}

for header, sequence in fasta_iter(contigs_file):
	contig_name = str(header).split(" ")[0] # take only the name without the flag

	if not contig_name in dictionary_of_contigs_length:
		dictionary_of_contigs_length[contig_name] = sequence
	
# open the cluster.tsv table
cluster_tsv = pd.read_csv(cluster_tsv, sep = "\t", names = ["bin", "contig"])

# create a dictionary with bin as a key and the nodes as key, only if the sum of the scaffold are >MIN_FASTA_LENGTH
bins_dictionary = {}

for i in range(len(cluster_tsv)):
	contig = cluster_tsv["contig"][i].split(" ")[0]
	sample = cluster_tsv["contig"][i].split(args.binsplit_char)[0]
	bin_id = "{}_{}".format(sample, cluster_tsv["bin"][i])

	if not bin_id in bins_dictionary:
		bins_dictionary[bin_id] = [contig]
	else:
		bins_dictionary[bin_id].append(contig)


# file with discarder genomes that do not agree with the threshold
to_write = "# minimum threshold {}\nBin_id\tlength\n".format(min_fasta)

# only create bins with length > MIN_FASTA_LENGTH
c = 1
bin_with_length_lest_than_minfasta = 0
remaining_bins  = 0

for bins in bins_dictionary:
	#print("Processing {}".format(bins))
	renamed_bins = output + "/{}.{}.fna".format(binning_method, bins)
	bins_length = 0

	for contig in bins_dictionary[bins]:

		if contig in dictionary_of_contigs_length:
			contig_length = int(len((dictionary_of_contigs_length[contig])))
			bins_length =  bins_length + contig_length
			
		else:
			print("{} Strange! Contigs cannot be found in the dictionary".format(contig))

			for i in dictionary_of_contigs_length:
				if contig in i:
					print("Contig in FASTA file: {}".format(i))
					print("Contig in VAMB file: {}".format(contig))
				
					sys.exit()
			sys.exit()

	if bins_length >= min_fasta:
		remaining_bins  = remaining_bins  + 1


		with open(renamed_bins, "w") as fh:
			for contig in bins_dictionary[bins]:
				fh.write(">{}\n".format(contig))
				sequence = dictionary_of_contigs_length[contig]
				fh.write("{}\n".format(str(sequence)))
		
	else:
		bin_with_length_lest_than_minfasta = bin_with_length_lest_than_minfasta + 1
		to_write = to_write + "{}\t{}\n".format(bins, bins_length )
	
	c = c + 1

print("# Total Bins: {}".format(len(bins_dictionary)))
print("# Bins < {}: {}".format(min_fasta, bin_with_length_lest_than_minfasta))
print("# Remaining Bins: {}".format(remaining_bins))


with open(discarded_genomes, "w") as fh:
	fh.write(to_write)
