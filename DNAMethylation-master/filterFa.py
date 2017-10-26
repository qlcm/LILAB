#!/usr/bin/python
# -*- coding: utf-8 -*-
"""filter Fasta file, dump the random fasta

Usage:
  filterFa.py [options] FASTA
  filterFa.py --version

Arguments:
  FASTA  path of fasta

Options:
  -h --help         show this help message and exit
  -o, --outdir      directory out
  -l, --len N       the shortest length of chrom, [default: 1000]
  -r, --reg=<arg>   regex, [default: random|Un_]
  -v, --version     show version and exit
"""
import os
import os.path
import sys
import re
try:
	from docopt import docopt
except ImportError:
	exit('This script requires package `docopt`\n'
		'sudo pip install docopt\n'
		'https://github.com/docopt/docopt')
try:
	from Bio import SeqIO
except ImportError:
	exit('This script requires package `Biopython`\n'
		'sudo pip install Biopython\n')

if __name__ == '__main__':
		arguments = docopt(__doc__, version='0.01')
		filterlen = arguments["--len"]
		reg = arguments["--reg"]
		filtername = re.compile(reg)
		if arguments["--outdir"]:
			dir_out = arguments["--outdir"]
		else:
			dir_out = os.getcwd()
		file_in = arguments["FASTA"]

		f_out = os.path.join(dir_out,os.path.basename(file_in).split('.')[0] +'.filter.fasta')
		sequences = []
		for record in SeqIO.parse(open(file_in), "fasta"):
			if len(record.seq) > int(filterlen) and len(filtername.findall(record.id)) == 0:
				print(record.id + " is OK" ) 
				sequences.append(record)

		output_handle = open(f_out, "w")
		SeqIO.write(sequences, output_handle, "fasta")
		output_handle.close()
