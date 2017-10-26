#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
from __future__ import division

"""process methylation"""

__appname__ = "methylation"
__author__  = "dmulilab"
__version__ = "0.0pre0"
__license__ = "GNU GPL 3.0 or later"

import re
import os
import os.path
import argparse
import sys
import csv
import numpy as np
from scipy.interpolate import UnivariateSpline
from intervaltree import Interval, IntervalTree
from scipy.stats.kde import gaussian_kde
from scipy.optimize import fsolve

import logging

#####################################################################
# common functions
#####################################################################

# log fuction

def init_log(logfilename):
	logging.basicConfig(level = logging.DEBUG, 
		format = '%(asctime)s %(message)s', 
		datefmt = '%Y-%m-%d %H:%M',
		filename = logfilename,
		filemode = 'w')

	console = logging.StreamHandler()
	console.setLevel(logging.INFO)
	formatter = logging.Formatter('%(message)s')
	console.setFormatter(formatter)
	logging.getLogger('').addHandler(console)

	return(logging.getLogger(''))

# load Methylation 

def load_methy_csv(filename):
	dictMethy = {}
	with open(filename, 'r') as methyFile :
		lines = csv.reader(methyFile, delimiter = ',')
		next(lines, None)
		for line in lines:
			chrname = line[0].strip()
			if not chrname in dictMethy:
				dictMethy[chrname] = []
			dictMethy[chrname] += [(int(line[1]), int(line[2]) + int(line[3]), int(line[4]) + int(line[5]))]
	methyFile.close()
	return(dictMethy)

# load the fa chrom sizes 

def load_chromsizes(filename):
	dictchromsizes = {}
	with open(filename, 'r') as sizesFile :
		lines = csv.reader(sizesFile, delimiter = '\t')
		for line in lines:
			chrname = line[0].strip()
			dictchromsizes[chrname] = line[1]
	sizesFile.close()
	return(dictchromsizes)

# get methylation vector

def get_meth_score(dictMethy, dictchromsizes):
	dictscore = {}
	for chrname in dictMethy:
		seqlen = int(dictchromsizes[chrname])
		vscore = [0] * seqlen
		for m in dictMethy[chrname]:
			vscore[m[0]] = m[1] * 1.0 / (m[1] + m[2]) if (m[1] + m[2]) else 0.0 
		dictscore[chrname] = vscore
	return(dictscore)

# get Methylation density by chromsom

def _guassian(x, mu, sigma):
	return (1 / np.sqrt(2 * np.pi * sigma * sigma)) * np.exp(-np.power(x - mu, 2.0) / (2 * sigma * sigma))

def get_methy_density(scorev, winsize, func = "guassian"):

	# step1: calculate convolution

	if(func == 'rect'):
		winv = [1 / winsize] * winsize
	else:
		winv = [_guassian(x, 0, winsize / 2.355) for x in range(-int(winsize * 1.7), int(winsize * 1.7))]
	methydensity = np.convolve(scorev, winv, mode = "same")

	return(methydensity)

# write methylation density

def write_density_csv(dictDensity, filename):
	try:
		csvFile = open(filename, 'w')
	except IOError:
		print('error: write to csv file "' + filename + '" failed!')
		sys.exit(-1)
	
	for chrname in dictDensity:
		csvFile.write('\n'.join([format('%s\t%f' % (chrname, density)) for density in dictDensity[chrname]]) + '\n')		
	csvFile.close()

def write_density_wig(dictDensity, filename):
	try:
		wigFile = open(filename, 'w')
	except IOError:
		print('error: write to wig file "' + filename + '" failed!')
		sys.exit(-1)
	
	for chrname in dictDensity:
		wigFile.write('fixedStep chrom=' + chrname + ' start=1 step=1' + '\n')
		wigFile.write('\n'.join([format(x) for x in dictDensity[chrname]]) + '\n')		
	wigFile.close()


def main():

	# parse command line options

	parser = argparse.ArgumentParser(description = '')
	parser.add_argument('mtbrfile', metavar = 'MtbrFile', 
		type = str, 
		help='mtbr file of the methylation')
	parser.add_argument('chromsizesfile', metavar = 'chromsizesfile', 
		type = str, 
		help='the genome fasta chrom sizes')
	parser.add_argument('-F', '--convfunc', dest = 'convfunc',
		type = str, default = 'guassian', 
		help = 'convolution function')
	parser.add_argument('-W', '--winsize', dest = 'winsize',
		type = float, default = 655.0, 
		help = 'convolution window size')

	args = parser.parse_args()

	baseFileName = os.path.splitext(os.path.basename(args.mtbrfile))[0]
	global log
	log = init_log(baseFileName + '.log')

	if(not os.path.exists(args.mtbrfile)):
		log.info('error: Mtbr file "', args.mtbrfile, '"', ' doest not exist.')
		sys.exit(-1)
	if(not os.path.exists(args.chromsizesfile)):
		log.info('error: chrom sizes file "', args.chromsizesfile, '"', ' doest not exist.')
		sys.exit(-1)

	# load  methylation mtbr file

	log.info('[*] loading methlation file')
	dictMethy = load_methy_csv(args.mtbrfile)

	# get chrom sizes

	log.info('[*] getting chrom sizes')
	dictchromsizes = load_chromsizes(args.chromsizesfile)

	# get methylation score

	log.info('[*] getting methlation score')
	dictscore = get_meth_score(dictMethy, dictchromsizes)

	# get methylation densities

	dictMethylation = {}
	log.info('[*] calculating methylation level ...')
	for chrname in dictscore:
		log.info('    calculating methylation level for chromsome ' + chrname)
		methydensity = get_methy_density(dictscore[chrname], args.winsize, args.convfunc)
		dictMethylation[chrname] = methydensity

	# write the output file

	log.info('[*]writting methylation level csv file')
	write_density_csv(dictMethylation, baseFileName + '.methylation.csv')

	log.info('[*]writting methylation level wig file')
	write_density_wig(dictMethylation, baseFileName + '.methylation.wig')
	log.info('[*] done')

if __name__ == '__main__':
	main()


		