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
import numpy as np
from scipy.interpolate import UnivariateSpline
from intervaltree import Interval, IntervalTree
from scipy.stats.kde import gaussian_kde
from scipy.optimize import newton

import logging

reSplitCG = re.compile('[ATN]|CG|[CG]')
reCGPos = re.compile('(CG)')

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

#######################################################################
## library: peakdetect
######################################################################

def _datacheck_peakdetect(x_axis, y_axis):
    if x_axis is None:
        x_axis = range(len(y_axis))
    
    if len(y_axis) != len(x_axis):
        raise (ValueError, 
                'Input vectors y_axis and x_axis must have same length')
    
    #needs to be a numpy array
    y_axis = np.array(y_axis)
    x_axis = np.array(x_axis)
    return x_axis, y_axis

def peakdetect(y_axis, x_axis = None, lookahead = 300, delta=0):
    """
    Converted from/based on a MATLAB script at: 
    http://billauer.co.il/peakdet.html
    
    function for detecting local maximas and minmias in a signal.
    Discovers peaks by searching for values which are surrounded by lower
    or larger values for maximas and minimas respectively
    
    keyword arguments:
    y_axis -- A list containg the signal over which to find peaks
    x_axis -- (optional) A x-axis whose values correspond to the y_axis list
        and is used in the return to specify the postion of the peaks. If
        omitted an index of the y_axis is used. (default: None)
    lookahead -- (optional) distance to look ahead from a peak candidate to
        determine if it is the actual peak (default: 200) 
        '(sample / period) / f' where '4 >= f >= 1.25' might be a good value
    delta -- (optional) this specifies a minimum difference between a peak and
        the following points, before a peak may be considered a peak. Useful
        to hinder the function from picking up false peaks towards to end of
        the signal. To work well delta should be set to delta >= RMSnoise * 5.
        (default: 0)
            delta function causes a 20% decrease in speed, when omitted
            Correctly used it can double the speed of the function
    
    return -- two lists [max_peaks, min_peaks] containing the positive and
        negative peaks respectively. Each cell of the lists contains a tupple
        of: (position, peak_value) 
        to get the average peak value do: np.mean(max_peaks, 0)[1] on the
        results to unpack one of the lists into x, y coordinates do: 
        x, y = zip(*tab)
    """
    max_peaks = []
    min_peaks = []
    dump = []   #Used to pop the first hit which almost always is false
       
    # check input data
    x_axis, y_axis = _datacheck_peakdetect(x_axis, y_axis)
    # store data length for later use
    length = len(y_axis)
    
    
    #perform some checks
    if lookahead < 1:
        raise ValueError, "Lookahead must be '1' or above in value"
    if not (np.isscalar(delta) and delta >= 0):
        raise ValueError, "delta must be a positive number"
    
    #maxima and minima candidates are temporarily stored in
    #mx and mn respectively
    mn, mx = np.Inf, -np.Inf
    
    #Only detect peak if there is 'lookahead' amount of points after it
    for index, (x, y) in enumerate(zip(x_axis[:-lookahead], 
                                        y_axis[:-lookahead])):
        if y > mx:
            mx = y
            mxpos = x
        if y < mn:
            mn = y
            mnpos = x
        
        ####look for max####
        if y < mx-delta and mx != np.Inf:
            #Maxima peak candidate found
            #look ahead in signal to ensure that this is a peak and not jitter
            if y_axis[index:index+lookahead].max() < mx:
                max_peaks.append([mxpos, mx])
                dump.append(True)
                #set algorithm to only find minima now
                mx = np.Inf
                mn = np.Inf
                if index+lookahead >= length:
                    #end is within lookahead no more peaks can be found
                    break
                continue
            #else:  #slows shit down this does
            #    mx = ahead
            #    mxpos = x_axis[np.where(y_axis[index:index+lookahead]==mx)]
        
        ####look for min####
        if y > mn+delta and mn != -np.Inf:
            #Minima peak candidate found 
            #look ahead in signal to ensure that this is a peak and not jitter
            if y_axis[index:index+lookahead].min() > mn:
                min_peaks.append([mnpos, mn])
                dump.append(False)
                #set algorithm to only find maxima now
                mn = -np.Inf
                mx = -np.Inf
                if index+lookahead >= length:
                    #end is within lookahead no more peaks can be found
                    break
            #else:  #slows shit down this does
            #    mn = ahead
            #    mnpos = x_axis[np.where(y_axis[index:index+lookahead]==mn)]
    
    
    #Remove the false hit on the first value of the y_axis
    try:
        if dump[0]:
            max_peaks.pop(0)
        else:
            min_peaks.pop(0)
        del dump
    except IndexError:
        #no peaks were found, should the function return empty lists?
        pass
        
    return [max_peaks, min_peaks]

#####################################################################
# common functions
#####################################################################

def load_reference_fa(filename):
	dictRefSeq = {}
	with open(filename, 'r') as refSeqFile :
		chrname = ''
		seq = ''
		for line in refSeqFile:
			if(line[0] == '>'):
				
				# save current seq for current chr

				if(chrname != ''):
					dictRefSeq[chrname] = seq
				
				# new chrname & seq

				chrname = line[1:].strip()
				seq = ''
				log.info('    loading reference sequence: ' + chrname)
			else:
				seq += line.strip().upper()
		
		# write the last chr

		if(chrname != ''):
			dictRefSeq[chrname] = seq
	refSeqFile.close()
	return(dictRefSeq)

def load_cg_island(filename):
	dictCGI = {}
	with open(filename, 'r') as cgiFile :
		lines = csv.reader(cgiFile, delimiter = '\t')
		next(lines, None)
		for line in lines:
			chrname = line[1].strip()
			if not chrname in dictCGI:
				dictCGI[chrname] = []
			dictCGI[chrname] += [(int(line[2]), int(line[3]), int(line[5]))]
	cgiFile.close()
	return(dictCGI)

def cgi_avg_len(dictCGI):
	cgisum = 0.0
	cgicount = 0.0
	for chrname in dictCGI:
		cgicount += len(dictCGI[chrname])
		for cgi in dictCGI[chrname]:
			cgisum += int(cgi[2])
	return (cgisum / cgicount if cgicount > 0 else 0)

###################################################################
## functions
###################################################################

# get CpG position vector in [0, 1, ...] format

def get_cg_pos(seq):
	cgpos = [m.start() for m in reCGPos.finditer(seq)]
	cgvector = [0] * len(seq)
	for pos in cgpos:
		cgvector[pos] = 1

	return(cgvector)

# get background CpG density

def _cg_background_avg(cgv):
	return([sum(cgv) * 1.0 / len(cgv) if len(cgv) > 0 else 0.0] * len(cgv))

def _shuffle(vector):
	random.shuffle(vector)
	return vector

def _random_cg(cgv, N):
	if(N == 0):
		return(0.0)

	sumv = numpy.array([0]*len(cgv))
	result = []
	pool = multiprocessing.Pool()
	for i in xrange(N): 
		result.append(pool.apply_async(shuffle, (cgv, )))
	pool.close()
	pool.join()
	for res in result:
		sumv +=  numpy.array(res.get())
	return sumv

def _cg_background_random(cgv, N):
	rndcgv = _random_cg(cgv, N % 10)
	for i in range(N / 10):
		rndcgv += _random_cg(cgv, 10)
	return(rndcgv / N)

def get_cg_background(cgv, N):
	if(N):
		return(_cg_background_random(cgv, N))
	else:
		return(_cg_background_avg(cgv))

# get CpG density by chromsom

def _guassian(x, mu, sigma):
	return (1 / np.sqrt(2 * np.pi * sigma * sigma)) * np.exp(-np.power(x - mu, 2.0) / (2 * sigma * sigma))

def get_cg_density(refseq, winsize, func = "guassian"):

	# step1: get CpG position

	cgv = get_cg_pos(refseq)

	# step2: calculate convolution

	if(func == 'rect'):
		winv = [1 / winsize] * winsize
	else:
		# FWHM = 2.355 * sigma, range = (-4 * sigma, 4 * sigma) 
		winv = [_guassian(x, 0, winsize / 2.355) for x in range(-int(winsize * 1.7), int(winsize * 1.7))]
	cgdensity = np.convolve(cgv, winv, mode = "same")

	# step3: get background cg density

	#cgbackground = get_cg_background(cgv, N)

	return(cgdensity)
	

# detect peaks and find the FWHM

def _peak_boundary(peaki, peakindexes, valleyindexes, direction, v):
	maxindex = len(peakindexes) - 1
	if (direction == "left") and (peaki == 0):
		return 0

	if (direction == "right") and (peaki == maxindex):
		endvalleys = [valleyindex for valleyindex in valleyindexes if valleyindex >= peakindexes[peaki]]
		if(len(endvalleys)):
			return (min(endvalleys))
		return (len(v) - 1)

	nearestvalley = None
	currentpeak = peakindexes[peaki]
	
	if direction == "left":
		nextpeak = peakindexes[peaki - 1]
		for valleyindex in valleyindexes[::-1]:
			if (valleyindex > nextpeak) and (valleyindex < currentpeak):
				nearestvalley = valleyindex
				break
	else:
		nextpeak = peakindexes[peaki + 1]
		for valleyindex in valleyindexes:
			if (valleyindex < nextpeak) and (valleyindex > currentpeak):
				nearestvalley = valleyindex
				break
	
	if nearestvalley :
		index = nearestvalley
	else:
		index = int((currentpeak + nextpeak) / 2.0)

	return (index)

def get_peaks(v, winsize, delta):
	winlen = int(winsize)
	maxindex = len(v) - 1
	peakvalues, valleyvalues = peakdetect(np.array(v), delta = delta, lookahead = winlen / 4.0)
	peakindexes = [value[0] for value in peakvalues]
	valleyindexes = np.array([value[0] for value in valleyvalues])
	peaks = []
	for i, peakindex in enumerate(peakindexes):
		
		# prepare data for peak fitting

		leftboundary = _peak_boundary(i, peakindexes, valleyindexes, "left", v)
		rightboundary = _peak_boundary(i, peakindexes, valleyindexes, "right", v)
		halfheight = v[peakindex] / 2.0
		peakdata = v[leftboundary:rightboundary]
		peakmax = np.max(peakdata)
		peakmean = np.mean(peakdata)

		# UnivariateSpline(k = 3), peakdata must have 4 site 
		if len(peakdata) < 4:
			peaks += [[peakindex, leftboundary, rightboundary, peakmax, peakmean, "PEAK"]]
			continue

		peakdata = [x - halfheight for x in peakdata]
		indexes = range(leftboundary, rightboundary)

		## find FWHM

		spline = UnivariateSpline(indexes, peakdata, s = 0)
		root = spline.roots()

		rootcount = len(root)
		if(rootcount == 0):
			r1 = (leftboundary + peakindex) / 2.0
			r2 = (peakindex + rightboundary) / 2.0
		elif(rootcount == 1):
			r = root[0]
			if( r < peakindex):
				r1 = r
				r2 = (peakindex + rightboundary) / 2.0
			else:
				r1 = (leftboundary + peakindex) / 2.0
				r2 = r
		else:
			distanceleft = [peakindex - r for r in root if r < peakindex]
			distanceright = [r - peakindex for r in root if r > peakindex]
			if(len(distanceleft) == 0):
				r1 = (leftboundary + peakindex) / 2.0
			else:
				r1 = peakindex - min(distanceleft) 
			if(len(distanceright) == 0):
				r2 = (peakindex + rightboundary) / 2.0
			else:
				r2 = peakindex + min(distanceright)
		rangestart = int(r1) if r1 > 0 else 0
		rangeend = int(r2) if r2 < maxindex else maxindex
		peakmax = v[peakindex]
		peakmean = np.mean(v[rangestart:(rangeend + 1)])
		peaks += [[peakindex, rangestart, rangeend, peakmax, peakmean, "PEAK"]]

	return peaks

def get_valley(v, peaks):
	valley = []
	i = 0
	for i in range(len(peaks) - 1):
		start = peaks[i][2]
		end = peaks[i + 1][1]
		valleystart = min(start, end)
		valleyend = max(start, end)
		if valleystart == valleyend:
			valleymin = v[valleystart]
			valleymean = v[valleystart]
			valleyindex = valleystart
		else:
			values = v[valleystart:valleyend]
			valleymin, valleyindex = min((val, idx) for (idx, val) in enumerate(values))
			valleymean = np.mean(values)
			valleyindex += valleystart
		valley += [[valleyindex, valleystart, valleyend, valleymin, valleymean, "VALLEY"]]

	return(valley)

# find the overlap

def _cgi_overlap(cgis, regions):
	cgiInterval = IntervalTree(Interval(cg[0], cg[1]) for cg in cgis)

	vcgi = []
	vnoncgi = []
	vvalley = []
	for region in regions:
		if region[5] == "VALLEY":
			vvalley += [region[4]]
		else:
			if cgiInterval.overlaps(region[1], region[2]):
				vcgi += [region[4]]
			else:
				vnoncgi += [region[4]]

	return(vcgi, vnoncgi, vvalley)

def _kde_peak(func):
	xs = np.linspace(0, 1, 1e4)
	ys = func(xs)
	index = np.argmax(ys)
	return(xs[index], ys[index])

def get_density_threshold(dictRegion, dictCGI):
	cgi = []
	noncgi = []
	valley = []
	for chrname in dictRegion:
		if(not chrname in dictCGI):
			continue
		cgis = dictCGI[chrname]
		region = dictRegion[chrname]
		vcgi, vnoncgi, vvalley = _cgi_overlap(cgis, region)
		cgi += vcgi
		noncgi += vnoncgi
		valley += vvalley

	# kde the cgi and noncgi

	cgipdf = gaussian_kde(cgi)
	noncgipdf = gaussian_kde(noncgi)
	valleypdf = gaussian_kde(valley)
	cgimax, cgimaxpeak = _kde_peak(cgipdf)
 	noncgimax, noncgimaxpeak = _kde_peak(noncgipdf)
 	valleymax, valleymaxpeak = _kde_peak(valleypdf)
	HMintersetion = newton(lambda x : cgipdf(x) - noncgipdf(x), x0 = 0.0, tol = 1e-6, maxiter = 10000)
	MLintersetion = newton(lambda x : noncgipdf(x) - valleypdf(x), x0 = 0.0, tol = 1e-6, maxiter = 10000)
	HMthreshold = (noncgimax + cgimax) / 2.0
	MLthreshold = (noncgimax + valleymax) / 2.0

	for x in HMintersetion:
		if (x > noncgimax) and (x < cgimax) :
			HMthreshold = x
			break
	for x in MLintersetion:
		if (x < noncgimax) and (x > valleymax) :
			MLthreshold = x
			break

	return(HMthreshold, MLthreshold, cgipdf, cgimax, noncgipdf, noncgimax, valleypdf, valleymax)		

def classfy_regions(dictRegion, HMthreshold, MLthreshold):
	for chrname in dictRegion:
		regions = dictRegion[chrname]
		for region in regions:
			if region[4] > HMthreshold:
				region += ["H"]
			elif (region[4] > MLthreshold) and (region[4] < HMthreshold):
				region += ["M"] 
			else:
				region += ["L"]

	return(dictRegion)

def write_regions_csv(dictRegion, filename):
	try:
		csvfile = open(filename, 'w')
	except IOError:
		log.info('error: write to csv file "' + filename + '" failed!')
		sys.exit(-1)
	csvfile.write('chr\tpos\tstart\tend\tmax/min\tmean\ttype\tclass\n')
	for chrname in dictRegion:
		for region in dictRegion[chrname]:
			csvfile.write(format('%s\t%d\t%d\t%d\t%f\t%f\t%s\t%s\n') % 
				(chrname, region[0], region[1], region[2], region[3], region[4], region[5], region[6]))
	csvfile.close()

def write_density_csv(dictDensity, filename):
	try:
		csvFile = open(filename, 'w')
	except IOError:
		log.info('error: write to csv file "' + filename + '" failed!')
		sys.exit(-1)
	
	for chrname in dictDensity:
		csvFile.write('chr\tdensity\n')
		chrcode = chrname.strip('chr')
		csvFile.write('\n'.join([format('%s\t%f' % (chrcode, density)) for density in dictDensity[chrname]]))		
	csvFile.close()

def write_kde_density(cgipdf, noncgipdf, valleypdf, cgimax, noncgimax, valleymax, HMthreshold, MLthreshold, filename):
	try:
		outfile = open(filename, 'w')
	except IOError:
		log.info('error: write to file "' + filename + '" failed!')
		sys.exit(-1)
	outfile.write(format('# cgipeak = %f cgipeakvalue = %f \n# noncgipeak = %f noncgipeakvalue = %f\n# valleypeak = %f valleypeakvalue = %f\n# HMthreshold = %f  MLthreshold = %f\n') % 
		(cgimax, cgipdf(cgimax), noncgimax, noncgipdf(noncgimax), valleymax, valleypdf(valleymax), HMthreshold, MLthreshold))
	xs = np.linspace(-0.99, 0.99, 2e4)
	cgiys = cgipdf(xs)
	noncgiys = noncgipdf(xs)
	valleyys = valleypdf(xs)
	outfile.write('\n'.join([str(x) + '\t' +str(y1) + '\t' + str(y2) + '\t' + str(y3) for x, y1, y2, y3 in zip(xs, cgiys, noncgiys, valleyys)]))
	outfile.close()

def write_regions_bed(dictRegion, filename):
	try:
		bedfile = open(filename, 'w')
	except IOError:
		log.info('error: write to file "' + filename + '" failed!')
		sys.exit(-1)
	dictColors = {"L":"255,0,0", "M":"0,255,0", "H":"0,0,255"}
	for chrname in dictRegion:
		for region in dictRegion[chrname]:
			bedfile.write(format('%s\t%d\t%d\t%s\t%d\t%s\t%d\t%d\t%s\n') % 
				(chrname, region[1], region[2], region[6], 0, '+' if region[5] == "PEAK" else '-', region[1], region[2], dictColors[region[6]]))
	bedfile.close()

def write_cgposition_wig(dictRefSeq, filename):
	try:
		wigFile = open(filename, 'w')
	except IOError:
		log.info('error: write to wig file "' + filename + '" failed!')
		sys.exit(-1)
	for chrname in dictRefSeq:
		vpos = get_cg_pos(dictRefSeq[chrname])
		wigFile.write('fixedStep chrom=' + chrname + ' start=1 step=1' + '\n')
		wigFile.write('\n'.join([format(x) for x in vpos]) + '\n')		
	wigFile.close()

def write_density_wig(dictDensity, filename):
	try:
		wigFile = open(filename, 'w')
	except IOError:
		log.info('error: write to wig file "' + filename + '" failed!')
		sys.exit(-1)
	
	for chrname in dictDensity:
		wigFile.write('fixedStep chrom=' + chrname + ' start=1 step=1' + '\n')
		wigFile.write('\n'.join([format(x) for x in dictDensity[chrname]]) + '\n')		
	wigFile.close()

def main():

	# parse command line options

	parser = argparse.ArgumentParser(description = '')
	parser.add_argument('infafile', metavar = 'FaFile', 
		type = str, 
		help='Fasta file of the reference genome')
	parser.add_argument('cgifile', metavar = 'cgifile', 
		type = str, 
		help='CpG island database file in csv format')
	parser.add_argument('-F', '--convfunc', dest = 'convfunc',
		type = str, default = 'guassian', 
		help = 'convolution function')
	parser.add_argument('-W', '--winsize', dest = 'winsize',
		type = float, 
		help = 'convolution window size')
	parser.add_argument('-D', '--delta', dest = 'delta',
		type = float,
		help = 'delta value for peak finder')

	args = parser.parse_args()

	# set up logging system

	baseFileName = os.path.splitext(os.path.basename(args.infafile))[0]
	global log
	log = init_log(baseFileName + '.log')

	# check commandline varabile

	if(not os.path.exists(args.infafile)):
		log.info('error: Reference sequence file "', args.infafile, '"', ' doest not exist.')
		sys.exit(-1)
	if(not os.path.exists(args.cgifile)):
		log.info('error: CpG island database file "', args.cgifile, '"', ' doest not exist.')
		sys.exit(-1)
	
	isWinSizeSet = (args.winsize is not None)
	isDeltaSet = (args.delta is not None)
	
	# load reference sequence

	log.info('[*] loading reference sequences')
	dictRefSeq = load_reference_fa(args.infafile)

	# load CpG Island & calculate convolution window size

	dictCGI = load_cg_island(args.cgifile)
	if isWinSizeSet:
		winsize = args.winsize
	else:
		winsize = cgi_avg_len(dictCGI)

	# get CpG densities

	dictDensity = {}
	log.info('[*] calculating CpG density ...')
	for chrname in dictRefSeq:
		log.info('    calculating CpG density for chromsome ' + chrname)
		log.info('    [window size = ' + str(winsize) + ']')
		cgdensity = get_cg_density(dictRefSeq[chrname], winsize, args.convfunc)
		dictDensity[chrname] = cgdensity

	# get CpG density peaks & valleys

	dictRegion = {}
	log.info('[*] spliting peaks and valleys ...')
	for chrname in dictDensity:
		log.info('    calculating Region for chromsome ' + chrname)
		density = dictDensity[chrname]
		if isDeltaSet:
			delta = args.delta
		else:
			delta = np.max(density) * 0.05
		log.info('    [peak detect delta = ' + str(delta) + ']')
		peaks = get_peaks(density, winsize = winsize, delta = float(delta))
		valleys = get_valley(density, peaks)
		dictRegion[chrname] = peaks + valleys

	# get overlaps with CpG island

	log.info('[*] getting CpG density threshold ...')
	HMthreshold, MLthreshold, cgipdf, cgimax, noncgipdf, noncgimax, valleypdf, valleymax = get_density_threshold(dictRegion, dictCGI)

	# annotate regions 

	log.info('[*] classifying regions ...')
	dictRegion = classfy_regions(dictRegion, HMthreshold, MLthreshold)

	# write output files

	log.info('[*] writting output files ...')

	log.info('    writting cg position wig file')
	write_cgposition_wig(dictRefSeq, baseFileName + '.cgpos.wig')

	log.info('    writting regions csv file')
	write_regions_csv(dictRegion, baseFileName + '.regions.csv')

	log.info('    writting density csv file')
	write_density_csv(dictDensity, baseFileName + '.density.csv')

	log.info('    writting CpG Island and CpG Density for Kernel Density Estimation')
	write_kde_density(cgipdf, noncgipdf, valleypdf, cgimax, noncgimax, valleymax, HMthreshold, MLthreshold, baseFileName + '.kde')

	log.info('    writting regions bed file')
	write_regions_bed(dictRegion, baseFileName + '.bed')

	log.info('    writting wig file')
	write_density_wig(dictDensity, baseFileName + '.cgden.wig')

	log.info('[*] done')

if __name__ == '__main__':
	main()
