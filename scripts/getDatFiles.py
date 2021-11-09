#!/usr/bin/python3

"""
A script to convert the liftover chain to a list of corresponding intervals in two species. 
"""

# Import packages
import subprocess
import argparse
import datetime
import time
import os
import operator
import numpy as np 
import gzip
import sys
from collections import Counter


class interval :
	"""
	A class of interval object, representing 2 corresponding intervals in two species
	"""
	def __init__(self, c1, s1, e1, cid, c2, s2, e2, strand, identifier) :
		# Chrom, start and end coordinate of the interval in specie 1 are stored as direct attributes of the object
		self.chrom = c1
		self.start = s1
		self.end = e1
		# The id of the chain from which the interval comes from is also stored as an attribute 
		self.chain_id = cid
		# The id of the interval (set manually), to be used to remove overlaps (see check_overlap function)
		self.interval_id = identifier
		# The coordinates of the corresponding interval in specie 2 (chrom, start, end and strand) are stored as a tuple in the inter attribute 
		self.inter = (c2, s2, e2, strand)
		
	
	def __eq__(self, other) : 
		return self.__dict__ == other.__dict__
	
	def __str__(self) : 
		return str(self.__dict__)

def load_regions(input_file, chains_score) :
	"""
	This function is made to load the liftOver chain file and translate them into a set of corresponding intervals in specie 1 and 2
	"""
	f = gzip.open(input_file, 'rt')
	# Data will be stored in a dictionnary ordered by specie 1 chromosomes 
	regions = {}
	identifier = 0 
	for line in f :
		# For some reason, there is an empty line at the end of each chain, skip it
		if line.strip() != "" :
			# Except for chain lines, chain files are in tab separated format 
			line = line.strip().split("\t")
			# New chain 
			if "chain" in line[0] :
				# Get data from chain line
				chain = line[0].split()
				score = int(chain[1])
				chrom_1 = chain[2]
				size_1 = int(chain[3])
				strand_1 = chain[4]
				start_1 = int(chain[5])
				end_1 = int(chain[6])
				chrom_2 = chain[7]
				size_2 = int(chain[8])
				strand_2 = chain[9]
				start_2 = int(chain[10])
				end_2 = int(chain[11])
				chain_id = int(chain[12])
				chains_score[chain_id] = score
				# Check that chrom 1 and 2 are well sequenced autosome (we didnt calculate NIEBs on the others). If so, put the flag to True to keep working with these data. If not, put the flag to False to discard these data.
				if len(chrom_1) < 6 and len(chrom_2) < 6 and chrom_1 != "chrY" and chrom_2 != "chrY" and chrom_1 != "chrX" and chrom_2 != "chrX" and chrom_1 != "chrM" and chrom_2 != "chrM" :
					flag = True
				else :
					flag = False

				if flag : 
					# For now, - strand in specie 1 is not managed. We'll see later if adding it is necessary (hope not)
					if strand_1 == "-" :
						print("ERROR : Minus strand for specie 1 in following chain : ")
						print(chain)
						print("Exiting")
						sys.exit()

					# Add the chrom from specie 1 to dictionnary if it is not already in it, and init a list to put the region objects
					if chrom_1 not in regions.keys() :
						regions[chrom_1] = []

					# Establish start specie 1 and specie 2 from line data
					# These variables will next be updated with each region
					# If both strand are +, we are in a classical case, both chromosome must be browsed from start to end 
					if strand_1 == "+" and strand_2 == "+" : 
						s1 = start_1
						s2 = start_2
						# s1 and s2 variables will then be updated by adding the size of each region to the variable (and of gaps if necessary). In the end, we should obtain a value corresponding to the end variable for each specie 
						
					# If strand for specie 2 is -, the region of specie 1 is aligned on the reverse complemented chromosome of specie 2
					# We must then browse the chromosome like in following schema :
					# + : 1 -------> 10  Classical case for specie 1
					# - : 10 <------- 1  Reverse complemented case for specie 2
					# Start (corresponding to 1 in previous schema) can be calculated like this : size_chromosome - start_position.
					elif strand_1 == "+" and strand_2 == "-" :
						s1 = start_1
						s2 = size_2 - start_2
						# In this case, s1 variable will be updated as in classical case while s2 variable will be updated by removing the size of each region, not adding it. In the end, we should obtain a value corresponding to size_chromosome - end_position
					else :
						print(chain)
			# New interval 	
			else :
				if flag :
					c1 = chrom_1
					identifier += 1 
					# Last entry seems to have no gap column
					if len(line) == 3 : 
						# Get data from line
						size = int(line[0])
						gap_1 = int(line[1])
						gap_2 = int(line[2])
					else :
						size = int(line[0])
						gap_1 = 0
						gap_2 = 0
					# Calculate end position of the interval in genome 1 and 2 
					# genome 1 
					e1 = s1 + size
					# genome 2
					c2 = chrom_2
					# If regions is + in genome 2, calculate end position as in genome 1
					if strand_2 == "+" : 
						e2 = s2 + size
					else :
						e2 = s2
						s2 = s2-size
					# Init a new inter object and put it in the specie 1 chromosome list 
					inter = interval(c1, s1, e1, chain_id, c2, s2, e2, strand_2, identifier)
					regions[chrom_1].append(inter)
					# Update s1 and s2 variable by adding (or removing if - strand) size of the interval + gap to next interval for the corresponding specie 
					s1 = e1 + gap_1
					if strand_2 == "+" : 
						s2 = e2 + gap_2
					else :
						s2 = s2 - gap_2

	f.close()
	# Sort each list of object by start position, to put interval in ascending order on chromosome 
	for c in sorted(regions.keys()) :
		# #print(c)
		regions[c].sort(key=lambda x:x.start)

	# At the end, we have a specie 1 chromosome-as-keys dictionnary, with, for each chromosome, a list of interval objects sorted by coordinates in ascending order. Each interval object represents a corresponding interval in specie 1 and 2
	return regions

def check_overlap(regions) :
	"""
	A simple function to check for overlap of intervals. The idea is to verify that none of the specie 1 locus are aligned at two loci in specie 2 and vice versa. It is basically a check for "one interval in specie 1 must correspond to one and one only interval in specie 2. To check that, we first verify that the end position of each interval is inferior or equal to the start position of the next interval in specie 1, and then re-sort intervals based on specie 2 positions, and do the same. If overlapping is found, all overlapping intervals are removed from the data and put in another file  
	"""

	## First, check overlaps in specie 1. No overlap is expected because liftover files are from specie 1 to specie 2, so logically, each specie 1 interval is aligned on one specific specie 2 interval, but only once. We want to check that, and add an overlap attribute to each interval object that will be set to True if the object overlap another, and to False if not 
	print("\t\tCheck overlap in specie 1...")
	for c in sorted(regions.keys()) : 
		for i in range(len(regions[c])-1) :
			# Case 1 : There is overlap between two intervals 
			if regions[c][i].end > regions[c][i+1].start :
				# From this, 2 cases are possible :
					# 1. The region i is also overlapping with region i-1. In this case, we want to add the region i+1 to region i overlapping_interval list attribute without removing i-1 region. We also want to create overlapping_interval list attribute for region i+1, with region i in it
					# 2. The region i is not overlapping with region i-1. In this case, we want to create overlap and overlapping_interval attributes for both region i and region i+1 objects. 
				if hasattr(regions[c][i], 'overlap') :
					# 1 : region i is overlapping with region i-1 
					# Set overlap attribute to True for region i+1. For region i, this attribute is already created and set to True so nothing to do here 
					regions[c][i+1].overlap = True 
					# Add region i+1 to region i overlapping_interval list attribute
					regions[c][i].overlapping_interval.append(regions[c][i+1])
					# Create region i+1 overlapping_interval list attribute with region i in it
					regions[c][i+1].overlapping_interval = [region[c][i]]
				else :
					# 2 : region i is not overlapping with region i-1 
					# Create overlap attributes for both regions and set them to True
					regions[c][i].overlap = True
					regions[c][i+1].overlap = True
					# Create overlapping_interval list attributes for both regions, with corresponding region in it (region i in regions i+1 attribute and vice-versa)
					regions[c][i].overlapping_interval = [regions[c][i+1]]
					regions[c][i+1].overlapping_interval = [regions[c][i]]
			# Case 2 : No overlap between region i and region i+1 
			elif regions[c][i].end <= regions[c][i+1].start :
				# Check the existence of attribute overlap in region i . Each region i has already been seen as a i+1 region in previous iteration. If it was concerned by an overlap back then, it should have the overlap attribute set to True. If it was not concerned by an overlap back then, it should not have overlap attribute at all. In this latter case, that means that region i is not overlapping with region i-1, neither does it with region i+1. So no overlap at all, then we can set its overlap attribute to False. For region i+1, we dont do anything, as it will be treated in next iteration, where it will be region i 
				if not hasattr(regions[c][i], 'overlap') :
					regions[c][i].overlap = False
		# As for loop goes only until len(regions[c]-1), if last interval of chromosome is not overlapping with previous interval, it will not have overlap attribute. We then need to create it and set it to False (it is not overlapping with previous interval and there is not next interval so no overlap here) 
		if not hasattr(regions[c][len(regions[c])-1], 'overlap') :
			regions[c][len(regions[c])-1].overlap = False

	# ~ n_inter = 0 
	# ~ for c in regions.keys() :
		# ~ n_inter += len(regions[c])
	# ~ print("Number of intervals in regions : ", n_inter)
	
	## At this point, we have identified the overlapping intervals in specie 1. Each interval has an overlap attribute, that is False if the interval has not overlap with another interval, and True if the interval has an overlapping interval. If there is an overlap, the overlapping interval of an interval is stored as an attribute. But at this point, no overlap is expected.

	## Now, we want to separate overlapping intervals and non-overlapping intervals, to write them in separated files. The idea is to remove, from regions dictionnary, the overlapping intervals and to store them in another dictionnary (overlapping regions) 
	print("\t\tRemove overlaps in specie 1...")
	overlapping_regions_s1 = {}
	regions = separate_overlaps(regions, overlapping_regions_s1) 

	# ~ n_inter = 0 
	# ~ for c in regions.keys() :
		# ~ n_inter += len(regions[c])
	# ~ print("Number of intervals in regions after remove intervals from specie 1: ", n_inter)
	
	## At this point, we have, in the regions dictionnary, a set of non-overlapping intervals, sorted by chromosomes and position on specie 1 genome. In the overlapping_regions dictionnary, we have a set of overlapping intervals, sorted in the same way. As no overlap were expected for specie 1, overlapping regions should be empty.
	

	## Now, we want to check overlaps in specie 2 on non-overlapping regions in specie 1. Overlap are expected here. Indeed, intervals of specie 1 are unique, but several specie 1 unique intervals can align on the same locus in specie 2. If so, it creates overlaps, with 2 regions of specie 1 corresponding to one region in specie 2. It can be because of duplication in specie 1, or deletions of duplicated regions in specie 2. As for now, we can not distinguish the duplicated regions from the original, we want to remove both regions from the data. The idea is to have a clean set of intervals aligned between specie 1 and specie 2, and another set of intervals containing all the overlapping ones. 

	print("\t\tRe-sort according to specie 2 coordinates...")
	# Re-sorting of regions dictionnary according to specie 2 positions instead of specie 1
	regions_resorted = resort_on_s2(regions)

	# ~ n_inter = 0 
	# ~ for c in regions_resorted.keys() :
		# ~ n_inter += len(regions_resorted[c])
	# ~ print("Number of intervals in regions_resorted after resorting of region on specie 2 coordinates : ", n_inter)

	print("\t\tCheck overlap in specie 2...")
	for c in sorted(regions_resorted.keys()) : 
		for i in range(len(regions_resorted[c])-1) :
			# Case 1 : There is overlap between two intervals 
			if regions_resorted[c][i].end > regions_resorted[c][i+1].start :
				# From this, 2 cases are possible :
					# 1. The region i is also overlapping with region i-1. In this case, we want to add the region i+1 to region i overlapping_interval list attribute without removing i-1 region. We also want to create overlapping_interval list attribute for region i+1, with region i in it
					# 2. The region i is not overlapping with region i-1. In this case, we want to create overlap and overlapping_interval attributes for both region i and region i+1 objects. 
				if hasattr(regions_resorted[c][i], 'overlap') :
					# 1 : region i is overlapping with region i-1 
					# Set overlap attribute to True for region i+1. For region i, this attribute is already created and set to True so nothing to do here 
					regions_resorted[c][i+1].overlap = True 
					# Add region i+1 to region i overlapping_interval list attribute
					regions_resorted[c][i].overlapping_interval.append(regions_resorted[c][i+1])
					# Create region i+1 overlapping_interval list attribute with region i in it
					regions_resorted[c][i+1].overlapping_interval = [regions_resorted[c][i]]
				else :
					# 2 : region i is not overlapping with region i-1 
					# Create overlap attributes for both regions and set them to True
					regions_resorted[c][i].overlap = True
					regions_resorted[c][i+1].overlap = True
					# Create overlapping_interval list attributes for both regions, with corresponding region in it (region i in regions i+1 attribute and vice-versa)
					regions_resorted[c][i].overlapping_interval = [regions_resorted[c][i+1]]
					regions_resorted[c][i+1].overlapping_interval = [regions_resorted[c][i]]
			# Case 2 : No overlap between region i and region i+1 
			elif regions_resorted[c][i].end <= regions_resorted[c][i+1].start :
				# Check the existence of attribute overlap in region i . Each region i has already been seen as a i+1 region in previous iteration. If it was concerned by an overlap back then, it should have the overlap attribute set to True. It it was not concerned by an overlap back then, it should not have overlap attribute at all. In this latter case, that means that region i is not overlapping with region i-1, neither does it with region i+1. So no overlap at all, then we can set its overlap attribute to False. For region i+1, we dont do anything, as it will be treated in next iteration, where it will be region i 
				if not hasattr(regions_resorted[c][i], 'overlap') :
					regions_resorted[c][i].overlap = False
		# As for loop goes only until len(regions[c]-1), if last interval of chromosome is not overlapping with previous interval, it will not have overlap attribute. We then need to create it and set it to False (it is not overlapping with previous interval and there is not next interval so no overlap here) 
		if not hasattr(regions_resorted[c][len(regions_resorted[c])-1], 'overlap') :
			regions_resorted[c][len(regions_resorted[c])-1].overlap = False

	print("\t\tRemove overlaps in specie 2...")
	overlapping_regions_s2 = {}
	regions_resorted = separate_overlaps(regions_resorted, overlapping_regions_s2) 

	# ~ n_inter = 0 
	# ~ for c in regions_resorted.keys() :
		# ~ n_inter += len(regions_resorted[c])
	# ~ print("Number of intervals in regions_resorted after remove overlaps from specie 2 : ", n_inter)

	# ~ n_inter = 0
	# ~ for c in overlapping_regions_s2.keys() :
		# ~ n_inter += len(overlapping_regions_s2[c])
	# ~ print("Number of overlapping intervals : ", n_inter)

	# At this point, we have in regions_resorted dictionnary all the non-overlapping intervals sorted according to specie 2 coordinates. We want to re-sort them according to specie 1 coordinates to write them in output file. We can use the resort_on_s2 function.
	print("\t\tRe-sorting of non-overlapping intervals according to specie 1 coordinates")
	regions = resort_on_s2(regions_resorted)

	print("\t\tRe-sorting of overlapping intervals according to specie 1 coordinates")
	overlapping_regions = merge_overlap(overlapping_regions_s1, overlapping_regions_s2)

	return overlapping_regions

	
def separate_overlaps(regions, olr) :
	"""
	Here, we want to separate the overlapping regions and the non-overlapping regions from the regions dictionnary in two distinct dictionnaries
	"""
	# Create non overlapping regions dictionnary 
	non_olr = {}
	# The idea is to browse the regions dictionnary, and put the intervals in the right dictionnary (overlapping or non overlapping) based on the overlap attribute
	for chrom in sorted(regions.keys()) :
		for interv in regions[chrom] :
			# If the interval has overlap
			if interv.overlap :
				# Add it to olr dictionnary 
				if chrom not in olr.keys() :
					olr[chrom] = [interv]
				else :
					olr[chrom].append(interv)
			# If the interval does not have overlap
			else :
				# Add it to non_olr dictionnary
				if chrom not in non_olr.keys() :
					non_olr[chrom] = [interv]
				else :
					non_olr[chrom].append(interv)
	# Return non_olr dictionnary to replace regions dictionnary in check_overlap function 
	return non_olr

def resort_on_s2(regions) :
	"""
	Here, we want to invert specie 1 and specie 2 of intervals in regions dictionnary, to re-sort it based on specie 2.
	For each interval, we then just want to add the interval to a new dictionnary, but based on specie 2 coordinates this time
	"""
	# Create new regions dictionnary
	new_regions = {}
	# Browse all intervals of regions dictionnary 
	for chrom in sorted(regions.keys()) :
		for interv in regions[chrom] :
			# Get infos for new interval 
			new_chrom = interv.inter[0]	# Get chrom of specie 2
			new_start = interv.inter[1]	# Get start of specie 2
			new_end = interv.inter[2]		# Get end of specie 2 
			strand = interv.inter[3]		# Get strand (not sure if really useful though)
			# Establish new interval based of specie 2 coordinates
			if new_chrom not in new_regions.keys() :
				new_regions[new_chrom] = [interval(new_chrom, new_start, new_end, interv.chain_id, chrom, interv.start, interv.end, strand, interv.interval_id)]
			else : 
				new_regions[new_chrom].append(interval(new_chrom, new_start, new_end, interv.chain_id, chrom, interv.start, interv.end, strand, interv.interval_id))

	# Sort the new dictionnary on coordinates like specie 1 based dictionnary 
	# Sort each list of object by start position, to put interval in ascending order on chromosome 
	for c in sorted(new_regions.keys()) :
		# #print(c)
		new_regions[c].sort(key=lambda x:x.start)

	# Return the new dictionnary to check overlaps in it
	return new_regions

def merge_overlap(olr1, olr2) :
	"""
	A function to sort all overlapping regions according to specie 1 coordinates, and merge the regions found in specie 1 and specie 2.
	As no regions are supposed to be found in specie 1, olr1 is supposed to be empty, and the merging step is kinda useless, but as I am not sure that it is inherent to all liftOver files, i choose the safety here and treat the case where overlaps have been found in specie 1 too.
	"""

	# First, re-sort olr2 dictionnary on specie 1 coordinates. The code is similar to the resort_on_s2 function, with conservation of overlapping intervals 
	olr2_rs = {}
	# Browse all intervals of regions dictionnary 
	for chrom in sorted(olr2.keys()) :
		for interv in olr2[chrom] :
			# First, get infos on specie 1 coordinates
			chrom_s1 = interv.inter[0]
			start_s1 = interv.inter[1]
			end_s1 = interv.inter[2]
			strand_s1 = interv.inter[3]
			# Second, get infos on specie 2 coordinates
			chrom_s2 = chrom
			start_s2 = interv.start
			end_s2 = interv.end

			# Establish new interval with infos gathered earlier 
			new_interval = interval(chrom_s1, start_s1, end_s1, interv.chain_id, chrom_s2, start_s2, end_s2, strand_s1, interv.interval_id)

			# Translate the overlapping intervals in specie 1 coordinates and add them to new_interval object
			new_interval.overlapping_interval = []
			for i in interv.overlapping_interval :
				new_i = interval(i.inter[0], i.inter[1], i.inter[2], i.chain_id, i.chrom, i.start, i.end, i.inter[3], i.interval_id)
				new_interval.overlapping_interval.append(new_i)
		
			# Add new_interval to the overlapping intervals dictionnary 
			if chrom_s1 not in olr2_rs.keys() :
				olr2_rs[chrom_s1] = [new_interval]
			else :
				olr2_rs[chrom_s1].append(new_interval)
	
	# At this point, we have resorted overlapping intervals in specie 2 from coordinates in specie 1. We then need to add eventual overlapping intervals in specie 1 to the data. It should not happen though, cause we should not have overlapping intervals in specie 1 
	for chrom in olr1.keys() :
		for interv in olr1[chrom] :
			if chrom not in olr2_rs.keys() : # This case should not happen 
				olr2_rs[chrom] = [interv]
			else : 
				olr2_rs[chrom].append(interv)

	# At this point, we have in olr2_rs dictionnary all the overlapping intervals, in specie 1 coordinates. We want to sort them according to these coordinates. 

	# Sort each list of object by start position, to put interval in ascending order on chromosome 
	for c in sorted(olr2_rs.keys()) :
		olr2_rs[c].sort(key=lambda x:x.start)

	# Return the new dictionnary of overlapping regions 
	return olr2_rs


def write_non_olr_intervals(regions, output) :
	"""
	A function to write the non overlapping intervals for each specie 1 chromosome.
	One file per chromosome, with the following tab separated format :
	Start s1    End s1    Chrom s2    Start s2    End s2    Strand s2    ID of chain
	"""
	for chrom in sorted(regions.keys()) :
		# One file per chromosome, named chromosome.dat. The date, species, etc are in the directory name 
		filename = output + chrom + ".dat"
		f = open(filename, 'w')
		for i in regions[chrom] :
			f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(str(i.start), str(i.end), i.inter[0], str(i.inter[1]), str(i.inter[2]), i.inter[3], str(i.chain_id)))
		f.close()

def write_olr_intervals(regions, output) : 
	"""
	A function to write the overlapping intervals for each specie 1 chromosome.
	One file per chromosome. On each line, coordinates of specie 1 interval, corresponding specie 2 interval, and all the overlapping couple of intervals (max 2).
	If some column have no data (for intervals with only one overlap for example), they are filled with NA
	"""
	for chrom in sorted(regions.keys()) :
		# One file per chromosome, named chromosome.dat. Dates, species and other informations are in directory name
		filename = output + chrom + ".dat"
		f = open(filename, 'w')
		for i in regions[chrom] :
			# One overlapping interval 
			if len(i.overlapping_interval) == 1 :
				olr = i.overlapping_interval[0]
				# Add the overlapping interval to the line, and fill the rest with NA 
				f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n".format(str(i.start), str(i.end), i.inter[0], str(i.inter[1]), str(i.inter[2]), i.inter[3], str(i.chain_id), str(olr.chrom), str(olr.start), str(olr.end), olr.inter[0], str(olr.inter[1]), str(olr.inter[2]), olr.inter[3], str(olr.chain_id)))
			# Two overlapping intervals 
			elif len(i.overlapping_interval) == 2 :
				olr1 = i.overlapping_interval[0]
				olr2 = i.overlapping_interval[1]
				# Add the two overlapping intervals to the line 
				f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(str(i.start), str(i.end), i.inter[0], str(i.inter[1]), str(i.inter[2]), i.inter[3], str(i.chain_id), str(olr1.chrom), str(olr1.start), str(olr1.end), olr1.inter[0], str(olr1.inter[1]), str(olr1.inter[2]), olr1.inter[3], str(olr1.chain_id), str(olr2.chrom), str(olr2.start), str(olr2.end), olr2.inter[0], str(olr2.inter[1]), str(olr2.inter[2]), olr2.inter[3], str(olr2.chain_id)))
			# More than two overlapping intervals should not be possible 
			elif len(i.overlapping_interval) > 2 :
				print("ERROR : More than two overlapping intervals !")
				print("Interval : ", i)
				for over in i.overlapping_interval :
					print("Overlapping interval : ", over)
			# No overlapping interval should not be possible as it should be considered in non overlapping data 
			elif len(i.overlapping_interval) == 0 :
				print("ERROR : No overlapping interval !")
				print(i)
		f.close()

# Main 
def main():
	
	parser = argparse.ArgumentParser()
	
	# Input
	parser.add_argument('-i', '--input_file', type=str, help='Path to input liftOver chain file', default ="/home/jbarbi02/These/DATA/LiftOver_chains/hg38ToPanTro5.over.chain.gz")
	
	# Output
	parser.add_argument('-o', '--output_dir', type=str, help='Path to output directory', default = "/home/jbarbi02/These/PROJECTS/PhylogeNIEB/2_Liftover_conversion/Test_2/")
	
	args = parser.parse_args()

	# Establish path to overlapping intervals and non-overlapping intervals output directories
	# Check that directory ends with '/' character
	if args.output_dir.endswith('/') :
		# If so, complete path 
		olr_path = args.output_dir + "Overlapping_regions/"
		non_olr_path = args.output_dir + "Non_Overlapping_regions/"
	else :
		#If not, add it to the argument before completing the path
		args.output_dir = args.output_dir + "/"
		olr_path = args.output_dir + "Overlapping_regions/"
		non_olr_path = args.output_dir + "Non_Overlapping_regions/"

	# Create output directories for overlapping and non-overlapping regions 
	subprocess.run(["mkdir", "-p", olr_path])
	subprocess.run(["mkdir", "-p", non_olr_path])

	
	# Load intervals
	print("1. Load liftover file and translate to intervals...")
	chains_score = {} # Init a dictionnary to keep chain scores. Format is : {id_1 : score_1, id_2 : score_2 ...}. This will be needed for the remove overlap function 
	regions = load_regions(args.input_file, chains_score)
	# Check interval overlapping
	print("\tLiftOver file loaded and translated\n2. Looking for overlaps in intervals...")
	olr = check_overlap(regions)
	# Write output
	print("\tOverlap removed\n3. Writing overlapping and non-overlapping intervals in output directory...")
	print("\tWriting non overlapping intervals")
	write_non_olr_intervals(regions, non_olr_path)
	print("\tWriting overlapping intervals")
	write_olr_intervals(olr, olr_path)
	print("\tData written\nAll done, exiting.")




if "__main__" == __name__:
	main()
