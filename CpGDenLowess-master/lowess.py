#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
from __future__ import division

"""process CpG density"""

__appname__ = "cgdensity"
__author__  = "dmulilab"
__version__ = "0.0pre0"
__license__ = "GNU GPL 3.0 or later"


import re
import os
import os.path
import argparse
import sys
import csv
import random
import gc
import math
import numpy as np
from math import ceil
from scipy import linalg
from scipy.interpolate import InterpolatedUnivariateSpline
import statsmodels.api as sm

import logging

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

def load_cg_csv(filename):
	cg = []
	with open(filename, 'r') as cgcsvFile :
		lines = csv.reader(cgcsvFile, delimiter = '\t')
		next(lines, None)
		for line in lines:
			cg += [[float(line[0]), float(line[1])]]
	cgcsvFile.close()
	return(cg)

def _rotate(x, y, x0, y0, theta):
	
	xt = x - x0
	yt = y - y0
	
	xr = xt * np.cos(theta) - yt * np.sin(theta)
	yr = xt * np.sin(theta) + yt * np.cos(theta)
	
	return(xr + x0, yr + y0)

def transform(x, y, x0, y0, delta, theta):
	
	xt, yt = _rotate(x + delta, y + delta, x0, y0, theta)

	return [xt, yt]

def write_csv(cgdata, filename):
	try:
		csvFile = open(filename, 'w')
	except IOError:
		log.info('error: write to csv file "' + filename + '" failed!')
		sys.exit(-1)
	
	csvFile.write('score\tdensity\n')	
	csvFile.write('\n'.join([format('%f\t%f' % (cg[0], cg[1])) for cg in cgdata]))		
	csvFile.close()

def write_cgtrans_wig(cgdata, filename, chrname, colindex):
	try:
		wigFile = open(filename, 'w')
	except IOError:
		log.info('error: write to wig file "' + filename + '" failed!')
		sys.exit(-1)

	wigFile.write('fixedStep chrom=' + chrname + ' start=1 step=1' + '\n')
	wigFile.write('\n'.join([format(cg[colindex]) for cg in cgdata]))
	wigFile.close()

def main():

	# parse command line options

	parser = argparse.ArgumentParser(description = '')
	parser.add_argument('cgcsvfile', metavar = 'CGCsvFile', 
		type = str, 
		help='CpG density and score csv file')
	parser.add_argument('-b', '--bin', dest = 'bin',
		type = int, default = 5000,
		help = 'bin count')
	parser.add_argument('-f', '--frac', dest = 'frac',
		type = float, default = 0.3333333,
		help = 'frac span')

	args = parser.parse_args()

	# set up logging system

	baseFileName = os.path.splitext(os.path.basename(args.cgcsvfile))[0]
	global log
	log = init_log(baseFileName + '.log')

	# check commandline varabile

	if(not os.path.exists(args.cgcsvfile)):
		log.info('error: CpG csv file "', args.cgcsvfile, '"', ' doest not exist.')
		sys.exit(-1)
	
	# load CpG csv file

	log.info('[*] loading CpG csv file')
	cg = np.array(load_cg_csv(args.cgcsvfile))
	density = cg[:, 0]
	score = cg[:, 1]

	# lowess fitting

	log.info('[*] performing lowess fitting')

	cglowess = sm.nonparametric.lowess(score, density, frac = args.frac)

	# build segements

	log.info('[*] building lowess segments')

	spl = InterpolatedUnivariateSpline(cglowess[:, 0], cglowess[:, 1])
	segscore = np.linspace(np.min(score), np.max(score), args.bin)
	segdensity = spl(segscore)
	seg = zip(segscore, segdensity)
	segcount = len(seg)
	for i in range(segcount - 1):
		segstart = seg[i]
		segend = seg[i + 1]
		a = (segend[1] - segstart[1]) * 1.0 / (segend[0] - segstart[0])
		delta = (1 - segstart[0] - segstart[1]) / 2.0
		theta = math.atan(a) - (135 * math.pi / 180)
		seg[i] += (delta, theta)

	print(seg)

	# tranform data points

	log.info('[*] transforming data points')
	
	cgtrans = []
	for score, density in cg:

		# locate segment 

		found = False
		for i  in range(segcount - 1):
			segstart = seg[i]
			segend = seg[i + 1]
			if(score >= segstart[0] and score <= segend[0]):
				found = True
				break;

		if( not found):
			log.info("error: data point (" + str(score) + "," + str(density), ") can not be located in segments")
			sys.exit(-1)

		# transform data point

		cgtrans += [transform(score, density, segstart[0], segstart[1], segstart[2], segstart[3])]

	log.info('[*] writting output files')

	log.info('    writting transformed csv file')
	write_csv(cgtrans, baseFileName + '.trans.csv')

	log.info('    writting wig file')
	write_cgtrans_wig(cgtrans, baseFileName + '.score.trans.wig', baseFileName, 0)
	write_cgtrans_wig(cgtrans, baseFileName + '.density.trans.wig', baseFileName, 1)

	log.info('[*] done')

if __name__ == '__main__':
	main()
