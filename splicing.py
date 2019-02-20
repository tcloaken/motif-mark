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

#####################
# Argument Parser   #
#####################
def get_args():
	parser = ap.ArgumentParser()
	parser = ap.ArgumentParser(description="This is a program that takes two arguments: a FASTA file and a simple text file with a list of splicing motifs.  Returns an SVG documents (one per gene) illustrating where the splicing motifs are in each gene")
	parser.add_argument("-F","--fasta", help="put the path the FASTA file after the flag followed by a space for the GENE sequence"
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
	
	#get list of motifs to find in the sequence
	motif_list = getMotifList([motif_in])

	for motif in motif_list:
		#find motif in the sequence
		l.append([m.start() for m in re.finditer(motif,sequence_in.lower())])
	
	return [item for sublist in l for item in sublist]


def getMotifList(List):
	#given a list of motifs, return a list without any
	#fancy "iupac" characters, by getting ther permutations of each motif
	fancy = {"r":["a","g"],"y":["c","t"],"s":["g","c"],"w":["a","t"],"k":["g","t"],
	"m":["a","c"],"b":["c","g","t"],"d":["a","g","t"],"h":["a","c","t"],
	"v":["a","c","g"],"n":["a","c","g","t"]}

	for pos,motif in enumerate(List):
		for position,letter in enumerate(motif):
			if letter in fancy:
				#keep working: remove the culprite motif from the List
				#replace it with as many instances of that iupac letter
				#as is in it's fancy dictionary
				motif_m = List.pop(pos)
				motif_m = list(motif_m) #make it mutable
				for newLetter in fancy[letter]:
					motif_m[position] = newLetter
					List.append("".join(motif_m))
				return getMotifList(List)
				
	return List


def RevComp(nucleotides):
	"""
	Given a DNA sequence (nucleotides) return the reverse compliment
	"""
	swap = {"A":"T","T":"A","G":"C","C":"G","N":"N"}#dictionary
	return ''.join([ swap[x] for x in nucleotides[::-1]])

def DrawIt(exons_pos,ex_lengths,motifs_list,gene,fileName):
	
	width = len(gene)+100
	height = 100+ (len(motifs_list)+2)*50
	
	surface = cairo.SVGSurface(fileName+".svg", width, height)
	context = cairo.Context(surface)
	context.set_line_width(5)
	context.move_to(0,height-50)
	context.line_to(len(gene), height-50)
	
	context.stroke()
	
	pat = cairo.LinearGradient(0.0,0.0,0.0,1.0)
	#pat.add_color_stop_rgba(0.2,1,0,0.5)
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
	
	#draw exons
	for start, lengths in zip(exons_pos,ex_lengths):
		context.rectangle(int(start[0]),height-100,int(lengths),100)
		context.set_source(pat)
		context.fill()
		context.stroke()
	
	
	#draw motifs

	for i,motif in enumerate(motifs_list):
		j = i+1
		positions = find_all(motif,gene)
		mot = cairo.LinearGradient(0.0,0.0,0.0,0.7)
		r,g,b = getNextColor(i)
		mot.add_color_stop_rgba(0,r/255,g/255,b/255,0.9)
		context.rectangle(5,30*j,12,12)
		context.set_source(mot)
		context.fill()
		context.move_to(20, 30*j + 12)
		context.set_source_rgb(0, 0, 0)
		context.show_text(motif)
		for pos in positions:
			context.rectangle(int(pos),height-100,len(motif),100)
			context.set_source(mot)
			context.fill()
			context.stroke()
	
	surface.finish()
		
	

def find_exon_intron_positions(genes):

	exon_pat = "[A-Z]+"
	intron_pat = "[a-z]+"
	
	#get list for exons and introns that have 
	exon_list = [re.findall(exon_pat,x) for x in genes]
	intron_list = [re.findall(intron_pat,x) for x in genes]
	
	exon_pos = []
	intron_pos = []
	for i,x in enumerate(exon_list):
		for exon in x:
			exon_pos.append([(m.start(), m.end()) for m in re.finditer(exon,genes[i])])
	
	for i,x in enumerate(intron_list):
		for introns in x:
			intron_pos.append([(m.start(), m.end()) for m in re.finditer(introns,genes[i])])	
	
	exon_lengths = []
	for x in exon_list:
		exon_lengths.append([ int(len(ex)) for ex in x])
	
	return(exon_pos,exon_lengths,intron_pos)
	
	
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
	exon_pos,exon_lengths,intron_pos = find_exon_intron_positions(genes)
	genes = [gene.lower() for gene in genes]
	for i,gene in enumerate(genes):
		DrawIt(exon_pos[i],exon_lengths[i],motif_list,gene,gene_names[i])
		
			
########
# MAIN #
########
main()
