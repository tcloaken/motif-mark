#!/usr/bin/env python3
"""
splicing marker


"""


####################
# Import libraries #
####################

import cairocffi as cairo
import argparse as ap
import re
import gzip
import time
import math
import random
import matplotlib as plt
from Bio import Seq
from itertools import product

#####################
# Argument Parser   #
#####################
def get_args():
	parser = ap.ArgumentParser()
	parser = ap.ArgumentParser(description="This is a program that takes two arguments: a FASTA file and a simple text file with a list of splicing motifs.  Returns an SVG documents (one per gene) illustrating where the splicing motifs are in each gene")
	parser.add_argument("-F","--fasta", help="put the path the FASTA file after the flag followed by a space for the GENE sequence, assumes genes are in the orientation"
	,type=str)
	parser.add_argument("-M","--motif", help="put the path the motif file after the flag followed by a space for the motif"
	,type=str)
	
	return(parser.parse_args())

####################
# GLOBAL VARIABLES #
####################

args = get_args()
fasta = args.fasta #fasta file
motifs = args.motif #motif

#set up color palate
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
rgb_colors = []
for c in colors:
	c = c.lstrip('#')
	rgb_colors.append(tuple(int(c[i:i+2], 16) for i in (0, 2 ,4)))

fancy = {"r":["a","g"],"y":["c","t"],"s":["g","c"],"w":["a","t"],"k":["g","t"],
	"m":["a","c"],"b":["c","g","t"],"d":["a","g","t"],"h":["a","c","t"],
	"v":["a","c","g"],"n":["a","c","g","t"], "u":["t"]}
#############
# FUNCTIONS #
#############

def getNextColor(index):
		#get color palate
	
	if index < len(rgb_colors):	
		return rgb_colors[index]
	else:
		index_mod = index%len(rgb_colors)
		r,g,b = rgb_colors[index_mod]
		r = max(min(random.normal(r,25),255),0)
		g = max(min(random.normal(g,25),255),0)
		b = max(min(random.normal(b,25),255),0)
		return (r,g,b)

def find_all(motif_in, sequence_in):
	#find all motifs by switching ambiguous nucleotides
	#calls motif_find on all variations and returns one single list in order
	""" taken from http://www.bioinformatics.org/sms2/iupac.html
	R.................A or G
	Y.................C or T
	S.................G or C
	W.................A or T
	K.................G or T
	M.................A or C
	B.................C or G or T
	D.................A or G or T
	H.................A or C or T
	V.................A or C or G
	N.................any base
	"""
	l = []	#list to hold the start positions
	motif_in = motif_in.upper()
	motif_in = motif_in.replace("U","T") #deals with RNA motifs
	#get list of motifs to find in the sequence
	motif_list = get_ambiguous(motif_in)

	for motif in motif_list:
		#find motif in the sequence
		l.append([m.start() for m in re.finditer(motif,sequence_in.upper())])
	
	flatten = [item for sublist in l for item in sublist]
	
	return flatten


def get_ambiguous(seq):
	# taken from "https://stackoverflow.com/questions/27551921/how-to-extend-ambiguous-dna-sequence"
	# and "https://github.com/biopython/biopython"
	# and "http://biopython.org/DIST/docs/api/Bio.Alphabet.IUPAC.IUPACAmbiguousDNA-class.html"
    d = Seq.IUPAC.IUPACData.ambiguous_dna_values
    ra = []
    for i in product(*[d[j] for j in seq]):
        ra.append("".join(i))
    return ra

def RevComp(nucleotides):
	"""
	Given a DNA sequence (nucleotides) return the reverse compliment
	"""
	swap = {"A":"T","T":"A","G":"C","C":"G","N":"N"}#dictionary
	return ''.join([ swap[x] for x in nucleotides[::-1]])

def DrawIt(exons_pos,ex_lengths,motifs_list,gene,fileName):
	
	width = len(gene)+150
	height = 100+ (len(motifs_list)+2)*50
	
	surface = cairo.SVGSurface(fileName+".svg", width, height)
	context = cairo.Context(surface)
	#draw Intron
	context.set_line_width(5)
	context.move_to(50,height-50)
	context.line_to(len(gene)+50, height-50)
	context.stroke()
	#exon color
	pat = cairo.LinearGradient(0.0,0.0,0.0,1.0)
	pat.add_color_stop_rgba(0,0.3,1,0.2,0.5)
	
	#draw legend
	context.rectangle(5,5,12,12)
	context.set_source(pat)
	context.fill()
	context.select_font_face("pyrus", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
	context.set_font_size(20)
	context.move_to(20, 18)
	context.set_source_rgb(0, 0, 0)
	context.show_text("Exon")
	context.move_to(85, 12)
	context.line_to(97,12)
	context.stroke()
	context.move_to(100,18)
	context.show_text("Intron")
	context.move_to(30,height-50)
	context.show_text("5'")
	context.move_to(len(gene)+60,height-50)
	context.show_text("3'")
	
	#draw exons
	for start, lengths in zip(exons_pos,ex_lengths):
		context.rectangle(int(start[0])+50,height-110,int(lengths),110)
		context.set_source(pat)
		context.fill()
		context.stroke()
	
	
	#draw motifs
	skipped = 0
	for i,motif in enumerate(motifs_list):
		j = i+1-skipped
		positions = find_all(motif,gene)
		if positions == []:
			#just skip the rest if there are no motifs in that gene
			skipped += 1
			continue
		#set color / legend for motif
		mot = cairo.LinearGradient(0.0,0.0,0.0,0.7)
		r,g,b = getNextColor(i)
		mot.add_color_stop_rgba(0,r/255,g/255,b/255,0.75)
		context.rectangle(5,30*j,12,12)
		context.set_source(mot)
		context.fill()
		context.move_to(20, 30*j + 12)
		context.set_source_rgb(0, 0, 0)
		context.show_text(motif)
		for pos in positions:
			context.rectangle(int(pos)+50,height-100,len(motif),100)
			context.set_source(mot)
			context.fill()
			context.stroke()

	
	
	

			
	surface.finish()
		
	

def find_exon_positions(genes):
	#exons need to be capitalized
	exon_pat = "[A-Z]+"
	
	#get list for exons and introns that have 
	exon_list = [re.findall(exon_pat,x) for x in genes]
	
	exon_pos = []
	
	for i,x in enumerate(exon_list):
		for exon in x:
			exon_pos.append([(m.start(), m.end()) for m in re.finditer(exon,genes[i])])
	
	exon_lengths = []
	for x in exon_list:
		exon_lengths.append([ int(len(ex)) for ex in x])
	
	return(exon_pos,exon_lengths)
	
	
def processFiles():
	with open(fasta,"r") as fh, open (motifs,"r") as mot:
		
		gene_titles = []
		genes = []
		string = ""
		for i,line in enumerate(fh):
			if (i == 0) and (">" in line):
				#start
				name = line.split(">")[1].rstrip()
				name = name.replace(" ","_")
				gene_titles.append(name)
			elif ">" in line:
				gene_titles.append(line.split(">")[1].rstrip())
				genes.append(string)
				string = ""
			else:
				string += line.rstrip()
		#append last exon
		genes.append(string)
		#exons are capitalized and introns are not
		
		
		motif_list = []
		for i,line in enumerate(mot):
			#loop through the motifs file and get the list 
			#of motifs
			line = line.rstrip()
			motif_list.append(line)
	
	return (motif_list,genes,gene_titles)
			
def main():
	"""
	Calls functions to execute program
	"""
	motif_list,genes,gene_names = processFiles()
	exon_pos,exon_lengths = find_exon_positions(genes)
	genes = [gene.lower() for gene in genes]
	for i,gene in enumerate(genes):
		DrawIt(exon_pos[i],exon_lengths[i],motif_list,gene,gene_names[i])
		
			
########
# MAIN #
########
main()
