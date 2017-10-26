#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
from __future__ import division

import csv

def load_csv(filename):
	dictGeneSwitch = {}
	with open(filename, 'r') as csvfile:
		lines = csv.reader(csvfile, delimiter = '\t')
		next(lines, None)
		for line in lines:
			tissue = line[2].strip()
			if tissue in dictGeneSwitch:
				tdata = dictGeneSwitch[tissue]
			else:
				tdata = {}

			gene = line[0]
			if not gene in tdata:
				tdata[gene] = []
			tdata[gene] += [(int(line[1]), int(line[3]))]

			dictGeneSwitch[tissue] = tdata
	csvfile.close()

	return(dictGeneSwitch)

def _has_open_switch(gene):
	for promoter in gene:
		if(promoter[1] == 0):
			return True
	return False

def _build_link(preGene, preGeneData, gene, promoter):
	link = '{"name":"flare.' + gene + '.'+ gene +'P' + str(promoter[0]) + '",';
	imports = []
	for prePromoter in preGeneData:
		if prePromoter[1] == 0:
			imports += ['"flare.'+ preGene + '.' + preGene + 'P'+ str(prePromoter[0]) + '"']
			break;
	link += '"imports":[' + ','.join(imports) + ']}'
	return link

def generate_d3js_json(dictGeneSwitch):
	links = []
	for tissue in dictGeneSwitch:
		tdata = dictGeneSwitch[tissue]
		genes = tdata.keys()
		genes.sort()
		for i, gene in enumerate(genes):
			
			# where this gene has open switches

			if not _has_open_switch(tdata[gene]):
				continue

			# fine the previous gene

			preGene = None
			j = i - 1
			if(j < 0):
				j += len(genes)
			while j != i:

				if _has_open_switch(tdata[genes[j]]):
					preGene = genes[j]
					break;

				j -= 1
				if(j < 0):
					j += len(genes)

			# build links

			if(preGene):
				for promoter in tdata[gene]:
					if(promoter[1] == 0):
						links += [_build_link(preGene, tdata[preGene], gene, promoter)]

	return ('[' + ','.join(links) + ']')						
														

dictGeneSwitch = load_csv("/home/fsch/Project/CpGDenLowess/test/bonemarrow/multipromoter/geneSwitch.csv")
print(dictGeneSwitch)
links = generate_d3js_json(dictGeneSwitch)

jsonFile = open("/home/fsch/Project/CpGDenLowess/test/bonemarrow/multipromoter/geneSwitch.d3js.json", 'w')
jsonFile.write(links)
jsonFile.close()