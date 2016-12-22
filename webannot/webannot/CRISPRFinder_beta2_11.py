#!/usr/bin/python3
import re
import os
import sys
import time
import datetime

################
# File Classes #
################

# fasta file, which contains header and sequence

class FastaSeq:
	def __init__ (self, header, sequence):
		self.header = header		# header of fasta file
		self.sequence = sequence 	# sequence itself
		

# contains list of all sequences in fasta file

class FastaList:
	def __init__ (self, file_path):
		self.fna_list =	[]		# list of fasta_seq objects
		
		file_content = open(file_path, 'r')
		header=""			# header of file
		sequence=""			# DNA sequence
		
		for file_line in file_content:
			if re.match('>', file_line):		# new sequence in multifasta
				if (header != "" and sequence != ""):	# if not first ">" in file
					self.fna_list.append(FastaSeq(header, sequence))	
					
				header = file_line[0:-1]	# last character is newline
				sequence = ""			# reset sequence if new file
			else:
				sequence = sequence + file_line[0:-1]	# append sequence part

		# append last element to list
		if (header!="" and sequence !=""):	# just in case
			self.fna_list.append(FastaSeq(header, sequence))
			

###########################
# Main CRISPRFinder Class #
###########################

# class containing all variables and functions necessary to find CRISPRs
#	-	file_path contains a path to .fasta file
#	-	output_file defines path where result directory will be saved
#	-	pattern is a combination of '#' (character) and '_' (space) used as seeds 
#		 during second pass. First pass used continuous seed of length k_mer_size (numbner of '#')
#	-	k_mer_size_filter is used to compare segments of DRs and spacers
#	-	window_size defines size of window in which k_mers will be searched for	
#	-	allowed_mismatch defines a number of mismatch allowed in repeat
#	-	spacer_dr_match_limit is a maximum number of matches of length 
#		 k_mer_size_filter per spacer between DR and spacer
# 	-	min_DR defines minimum length in bp of repeat part
#	-	max_DR defines maximum length in bp of repeat part
#	-	min_spacer_DR_ratio is minimum quotient allowed length of spacer versus DR
#	-	max_spacer_DR_ratio is maximum quotient allowed length of spacer versus DR
#	-	first_pass_limit is maximum allowed distance between two regions with repeats
#	-	window_container holds all k-mers currently detected in window, whose 
#		 position is saved in lookup table
#	-	k_mer_counts are counts for each k-mer for sequence
#	-	repeat_list_first_pass contains "interesting zones" with possible repeats 
#		 (thus possible CRISPRs)
#	-	identity_limit specifies limit necessary to extend a repeat during alignment


class FindCRISPRs:
	def __init__ (self, file_path, output_path, k_mer_size_filter, pattern,\
	 window_size, allowed_mismatch, spacer_dr_match_limit, min_DR, max_DR,\
	 min_spacer_DR_ratio, max_spacer_DR_ratio, first_pass_limit, search_tracrRNA):

		# file path
		self.fastas = FastaList(file_path)
	
		# pattern
		self.pattern = Pattern(pattern)

		# CRISPRFinder parameters
		self.k_mer_size = len(self.pattern)
		if self.k_mer_size < 5:
			print("Error: Pattern has too low number of accepted symbols '#'; please enter the pattern with at least 5.")
			sys.exit()
		self.k_mer_size_filter = k_mer_size_filter
		if self.k_mer_size_filter > self.k_mer_size:
			print("Warning!! Size of filtering k-mer bigger than number of accepted symbols '#' in pattern, setting value to %d"%(self.k_mer_size - 1))
		self.window_size = window_size
		self.allowed_mismatch = allowed_mismatch
		self.spacer_dr_match_limit = spacer_dr_match_limit
		self.min_DR = min_DR 
		self.max_DR = max_DR
		self.min_spacer_DR_ratio = min_spacer_DR_ratio
		self.max_spacer_DR_ratio = max_spacer_DR_ratio
		self.first_pass_limit = first_pass_limit
		self.threshold = 0.5
		self.identity_limit = 1
		self.identity_truncated = 0.5 
		
		self.search_tracrRNA = search_tracrRNA

		# variables used to stock results during first pass

		self.window_container = {}		# stocks association k-mer - count
		self.lookup_table = []			# stocks k-mers in order
		self.k_mer_counts = []			# stocks k-mer counts in order
		self.repeat_list_first_pass = []	# stocks possible repeats (objects of Repeat type)

		# variables used to stock results during second pass of analysis		

		self.local_container = {}		# stocks association k-mer - count for 'interesting zone'
		self.local_lookup_table = [] 		# k-mers in order
		self.local_counts = []			# stocks k-mer counts in order
		self.clusters = []			# list of clusters

		self.CRISPRs = []			# list of CRISPR results for actually analyzed sequence
		self.CRISPRs_all = [[] for i in range(len(self.fastas.fna_list))]
		 # list of CRISPRs for all sequences in file
		
		# this method will initialize window container with k-mer occurences and positions
		# once i == window_length/2, results are exported from dictionnary to self.k_mer_counts

		# path to save outputted files
		self.output_path = output_path


	def initiate_window(self, fasta_seq):
		# fasta_seq is FastaSeq object
		
		if len(fasta_seq.sequence) < self.window_size:
			self.window_size = len(fasta_seq.sequence)
			print("Warning, sequence is shorter than desired window length. Window length will be set to %d" %(len(fasta_seq.sequence)))

		pos = 0
		window_middle = int(self.window_size/2)
		
		while pos < self.window_size:
			key_sequence = fasta_seq.sequence[pos : (pos+self.k_mer_size)]
			self.lookup_table.append(key_sequence)
	
			if key_sequence in self.window_container:
				self.window_container[key_sequence][0] = self.window_container[key_sequence][0]+1	# update count
				self.window_container[key_sequence][1].append(pos)					# update list of positions
			else:
				self.window_container[key_sequence] = [1, [pos]]	
			
			# once i == window_length/2, results are exported from dictionnary to self.k_mer_counts
			if pos >= window_middle :
				window_middle_position = pos - window_middle
				self.k_mer_counts.append(self.window_container[self.lookup_table[window_middle_position]][0])	

			pos += 1

		return pos


	# moves window by jump_size
	# inserts new key if absent, updates if present, and removes sorting entry
	# position : where window ends
	def move_window(self, fasta_seq, position):

		key_sequence = fasta_seq.sequence[position:(position+self.k_mer_size)]
		self.lookup_table.append(key_sequence)	

		if key_sequence in self.window_container:
			self.window_container[key_sequence][0] = self.window_container[key_sequence][0]+1	# update count
			self.window_container[key_sequence][1].append(position)					# update list of positions
		else:
			self.window_container[key_sequence] = [1, [position]]
	
		# window start position = window end position - window length		
		window_starting_position = position - self.window_size
		window_middle_position = position - int(self.window_size/2)

		self.k_mer_counts.append(self.window_container[self.lookup_table[window_middle_position]][0])

		# if only one entry, delete this key
		if self.window_container[self.lookup_table[window_starting_position]][0]==1:
			del self.window_container[self.lookup_table[window_starting_position]]
		else:
			self.window_container[self.lookup_table[window_starting_position]][0] =\
			 self.window_container[self.lookup_table[window_starting_position]][0] - 1
			del self.window_container[self.lookup_table[window_starting_position]][1][0]				# remove first element from list

		return position + 1

		
	# transforms dictionnary to list of frequencies when end of sequence is reached
	def finish(self, fasta_seq):
		for pos in range(len(fasta_seq.sequence) - int(self.window_size/2), len(fasta_seq.sequence)):
			self.k_mer_counts.append(self.window_container[self.lookup_table[pos]][0])


	# this function uses previous three (in that order) to obtain k-mer counts in range of window for entire sequence
	def get_sequence_counts(self, fasta_seq):
		# reset stocking tables
		self.window_container = {}
		self.lookup_table = []
		self.k_mer_counts = []
		self.repeat_list_first_pass = []

		# initiate window
		pos = self.initiate_window(fasta_seq)
		# get counts for entire sequence
		while pos < len(fasta_seq.sequence):
			pos = self.move_window(fasta_seq, pos)
			
		# transcribe window dictionary in counts
		self.finish(fasta_seq)


	# finds and returns expanded zones containing repeated segments in fasta_seq
	def get_repeats(self, fasta_seq):
		pos = 0
		repeats = []
		precursor_zone = 0			# indicates if inside of possible CRISPR precursor
		precursor_start = 0			# start of CRISPR precursor 
		precursor_end = 0			# end of CRISPR precursor
		after_repeat_length = 0	 	# length of mismatched space after last repeat
			
		while pos < len(self.k_mer_counts):

			# anytime there is a repeated sequence found, a search is made until either desired length is found (accepted) or number of mismatches is met (dismissed)
			# added case when 		
			if self.k_mer_counts[pos] > 1 and pos < (len(self.k_mer_counts) - 1):
				if precursor_zone == 0:				# only if we're not in repeat zone yet (ie. first repeat of local group)
					precursor_start = pos			# where possible interesting zone containing repeat(s)
					
				consecutive_mismatch = 0			# number of mismatches found in row (not important if 0 or 1 mismatch allowed only, in other cases however...)
				total_repeat_length = self.k_mer_size		# length of repeat (remember)
				
				while consecutive_mismatch < (self.allowed_mismatch + self.k_mer_size)\
				 and pos < (len(self.k_mer_counts) - 1):
					if self.k_mer_counts[pos] <= 1:
						consecutive_mismatch += 1
						
					else:
						total_repeat_length = total_repeat_length + 1 + consecutive_mismatch
						consecutive_mismatch = 0
						
					pos += 1
				
				if total_repeat_length > self.min_DR/2 and total_repeat_length < self.max_DR:
					# if repeat is long enough (but not too long), we mark start of potential precursor to track next repeats
					precursor_zone = 1
					precursor_end = pos - consecutive_mismatch		# marks actual end of zone if repeat is sufficiently long
					after_repeat_length = consecutive_mismatch		# new space between repats is defined by number of mismatches at the extremity
					
				else:
					if precursor_zone == 1:
						after_repeat_length += total_repeat_length + consecutive_mismatch
			
			else:
				if precursor_zone == 1:
					# if space without repeats is longer than allowed or the sequence is ending, export repeat
					if after_repeat_length >= self.first_pass_limit or pos == len(self.k_mer_counts) - 1:
						self.repeat_list_first_pass.append(\
						 Repeat(fasta_seq.sequence[precursor_start : (precursor_end+self.k_mer_size-1)],\
						 precursor_start, precursor_end + self.k_mer_size-1))

						precursor_zone = 0				# once saved we search for next extended zone with repeats
						after_repeat_length = 0			# we leave zone in proximity of repeats (defined by self.first_pass_limit)

					else:
						after_repeat_length += 1
					
				pos += 1


	# uses sliding window and then searches for regions with repeats
	def first_pass (self, fasta_seq):
		self.get_sequence_counts(fasta_seq)
		self.get_repeats(fasta_seq)
		

	# counts repeats locally using other pattern more friendly to mismatches. Precises CRISPRs
	def get_pattern_counts (self, repeat):
		# reset stocking tables
		self.local_container = {}
		self.local_lookup_table = []
		self.local_counts = []
		self.repeat_list_second_pass = []
		
		pos = 0
		while pos<len(repeat):
			key_sequence = self.pattern.extract_pattern(repeat.sequence, pos)
			self.local_lookup_table.append(key_sequence)
	
			if key_sequence in self.local_container:
				self.local_container[key_sequence][0] = self.local_container[key_sequence][0]+1		# update count
				self.local_container[key_sequence][1].append(pos)					# update list of positions
			else:
				self.local_container[key_sequence] = [1, [pos]]

			pos += 1
	
		for patt in self.local_lookup_table:
			self.local_counts.append(self.local_container[patt][0])
		

	# extracts direct repeats for sequences and stocks them into clusters
	def extract_clusters (self, repeat):

		pos = 0					# position within repeat
		DR_pos = 0				# position within DR
		max_kmer_found = 0 			# indicates if maximum kmer was found (k_mer with highest number of repeats)
		max_kmer_exported = 0			# indicates if maximum kmer was exported (in case if repeat turns out to be false - not short or long enough)
		max_kmer_info = []	
		
		repeat_start = repeat.begin		# true starting position of repeat (needed to offset local position pos)	
		max_value = max(self.local_counts)	# maximum number of repeats in table
		threshold_value = round(max_value*self.threshold)	# limit value to distinguish repeats and unique zones (though unique should mean 1)	
		if threshold_value < 2:
			 threshold_value = 2

		spacer_length = 0	# distance between two repeats
		total_repeat_length = 0
		in_cluster = 0		# indicates whether we are in cluster or not
			
		while pos < len(self.local_counts):
			
			repeat_starting_pos = pos
			if self.local_counts[pos] >= 1:

				DR_pos = 0				
				mismatch_count = 0				# total number of mismatches
				consecutive_mismatch = 0			# number of mismatches in a row

				total_repeat_length = self.pattern.n_symbols	# starting length is that of k-mer				
				
				# while less than allowed number of mismatches and sequence is not ending			
				while (consecutive_mismatch < (self.allowed_mismatch+self.pattern.n_symbols)) and (pos < len(self.local_counts)):
					if self.local_counts[pos] < threshold_value:
						consecutive_mismatch += 1
					else:
						total_repeat_length = total_repeat_length + 1 + consecutive_mismatch
						consecutive_mismatch = 0
		
						# first time maximum value is found and not exported yet (maximum repeat)
						if self.local_counts[pos] == max_value and max_kmer_found == 0 and max_kmer_exported == 0:

							max_kmer_found = 1
							actual_repeat = 0
							
							if in_cluster == 0:
								actual_repeat = 0
								
							else:
								actual_repeat = len(self.clusters[-1])
								
							max_kmer_info = [self.pattern.extract_pattern(repeat.sequence, pos), actual_repeat, DR_pos]
							
					pos += 1
					DR_pos += 1						
	
				if total_repeat_length >= self.min_DR and total_repeat_length <= self.max_DR:

					if in_cluster == 0:
						self.clusters.append(Cluster())	
						in_cluster = 1					
					
					spacer_length = consecutive_mismatch
					repeat_ending_pos = repeat_starting_pos + total_repeat_length
					self.clusters[-1].add_repeat(Repeat(repeat.sequence[repeat_starting_pos:repeat_ending_pos], repeat_start + repeat_starting_pos + 1, repeat_start + repeat_ending_pos))	
					if max_kmer_found == 1 and max_kmer_exported == 0:	# kmer was found but not exported yet
			
						max_kmer_exported = 1
						self.clusters[-1].most_repeated = max_kmer_info
	
				else:
					spacer_length += total_repeat_length + consecutive_mismatch
					max_kmer_found = 0

			elif self.local_counts[pos] < threshold_value:
				if in_cluster == 1:
					if (spacer_length <= self.max_DR*self.max_spacer_DR_ratio) and self.clusters[-1]:
						spacer_length += 1
						
					elif (spacer_length > self.max_DR*self.max_spacer_DR_ratio):	
						# if space too long, go to new cluster
						spacer_length = 0
						# if last one consist of only one repeat, delete it
						if len(self.clusters[-1]) == 1 or None in self.clusters[-1].most_repeated:
							del self.clusters[-1]
							
						in_cluster = 0
		
						max_kmer_found = 0		# resetting values for k-mers
						max_kmer_exported = 0
				
				pos += 1

		if self.clusters:
			if len(self.clusters[-1]) == 1  or None in self.clusters[-1].most_repeated:
				del self.clusters[-1]
				

	# aligns segments of repeats of the same cluster
	# uses k-mers appearing most frequently as seeds
	def align_seeds (self, cluster):

		seed_seq = cluster.most_repeated[0]			# most frequently appearing sequence 
		ref_number = cluster.most_repeated[1]			# number of reference repeat
		ref_seed_pos = cluster.most_repeated[2]			# position of most frequent sequence

		seed_found = 0						# indicates if seed was found for given repeat

		seed_positions = [ref_seed_pos]				# position at which the seeds start
		
		comp_number = ref_number + 1				# sequence to compare
		while comp_number < len(cluster):
			
			seed_found = 0
			comp_seed_pos = 0
			while comp_seed_pos < (len(cluster.repeats[comp_number].sequence) - self.pattern.n_symbols) and seed_found == 0:
				if seed_seq == self.pattern.extract_pattern(cluster.repeats[comp_number].sequence, comp_seed_pos):
					seed_found = 1	
					seed_positions.append(comp_seed_pos)
					break
				
				comp_seed_pos += 1

			if seed_found == 0:
				seed_positions.append(-1)

			comp_number += 1

		for i in range(len(seed_positions) - 1, 0, -1):
			if seed_positions[i] == -1:
				del seed_positions[i]
				del cluster.repeats[i + ref_number]				

		if ref_number > 0:	
			for rep_num in range(0, ref_number):	
				del cluster.repeats[0]	

		return seed_positions


	# extends using seeds, extends repeats one by one
	# needs the result from previous function
	# works only if cluster has at least two repeats
	def extend (self, cluster, seed_positions, fasta_seq):
		# initializing mismatch number parameter
		column_mismatch_limit = 0
		total_mismatch_limit = 0
		
		# setting limit to number of mismatches depending on length of cluster
		if len(cluster) < 5:
			column_mismatch_limit = self.identity_limit
		else:
			column_mismatch_limit = self.identity_limit*0.15*len(cluster)

		if len(cluster) == 1:
			self.clusters.remove(cluster)
		else:			
			# first part : filling eventual gaps in spaces of pattern
			total_mismatch_count_left = 0
			total_mismatch_count_right = 0

			for i in self.pattern.return_blanks():
				blank_nucls = [cluster.repeats[x].sequence[i+seed_positions[x]] for x in range(0, len(cluster))]
				most_frequent = max(set(blank_nucls), key = blank_nucls.count)
				if (len(blank_nucls) - blank_nucls.count(most_frequent)) <= self.identity_limit:
					total_mismatch_count_left += len(blank_nucls) - blank_nucls.count(most_frequent)
					
				else:
					# major disagreement between nucleotides in the middle of the pattern - big probability solution is wrong, so it gets erased instead
					self.clusters.remove(cluster)
					return 
					
			total_mismatch_count_right = total_mismatch_count_left		
					
			# second part : extending results themselves
			
			before_seed = 0
			consecutive_mismatch_before = 0 	# used in shorter CRISPRs since these can have mismatches (case extension in front of seed)
			
			# setting up a limit to number of mismatches occuring in CRISPR
			upper_mismatch_limit = (self.allowed_mismatch) * (len(cluster) - 1)

			# check for extension before
			while (total_mismatch_count_left <= upper_mismatch_limit):
					
				if min([(cluster.repeats[x].begin + seed_positions[x] + before_seed - 1)\
				 for x in range(0, len(cluster))]) >= 0:
				# ^ -1 is because cluster.repeats[x].begin is defined as position in sequence, not in list
					before_seed -= 1
					# determine the column of supposedly aligned nucleotides, then calculate mismatches
					extend_nucls_before = [fasta_seq.sequence[cluster.repeats[x].begin +\
					 seed_positions[x] + before_seed - 1] for x in range(0, len(cluster))]
					[mismatches, max_counts] = self.count_mismatches(extend_nucls_before)
					# if either percentage or number of mismatches is lower than limit, extend
					most_frequent = max(set(extend_nucls_before), key = extend_nucls_before.count)
					if (len(extend_nucls_before) - max_counts) <= self.identity_limit or mismatches < column_mismatch_limit:
						total_mismatch_count_left += mismatches
						consecutive_mismatch_before = 0
							
					else:
						if consecutive_mismatch_before == 1:
							before_seed += 2
							consecutive_mismatch_before = 0
						else:
							before_seed += 1
						break # end while
							
				else:
					before_seed += 1
					break # end while		
					
			after_seed = 0
			consecutive_mismatch_after = 0 		# same (case extension after seed)
									
			# alternating between checking for extensions before seed and after seed
			while (total_mismatch_count_right <= upper_mismatch_limit):						

				# checking nucleotides after seed
				if max([(cluster.repeats[x].begin + seed_positions[x] + self.pattern.n_symbols + after_seed - 1)\
				 for x in range(0, len(cluster))]) < len(fasta_seq.sequence):
				# ^ -2 because of string -> list conversion (-1) and since we want to add self.pattern.n_symbols - 1		
					after_seed += 1
					# determine the column of supposedly aligned nucleotides, then calculate mismatches 			
					extend_nucls_after = [fasta_seq.sequence[cluster.repeats[x].begin +\
					 seed_positions[x] + self.pattern.n_symbols + after_seed - 2]\
					 for x in range(0, len(cluster))]	
					[mismatches, max_counts] = self.count_mismatches(extend_nucls_after)
					# if either percentage or number of mismatches is lower than limit, extend
					most_frequent = max(set(extend_nucls_after), key = extend_nucls_after.count)
					if (len(extend_nucls_after) - max_counts) <= self.identity_limit or mismatches < column_mismatch_limit:		
						total_mismatch_count_right += mismatches
						consecutive_mismatch_after = 0
												
					else:
						if consecutive_mismatch_after == 1:
							after_seed -= 2
							consecutive_mismatch_after = 0
						else:
							after_seed -= 1
						break # end while
				else:
					after_seed -= 1
					break # end while
						

			if consecutive_mismatch_before == 1:
				before_seed += 1
			
			if consecutive_mismatch_after == 1:
				after_seed -= 1

			for i in range(0, len(cluster)):
				new_begin = cluster.repeats[i].begin + seed_positions[i] + before_seed
				new_end = cluster.repeats[i].begin + seed_positions[i] + self.pattern.n_symbols + after_seed - 1
				new_sequence = fasta_seq.sequence[(new_begin - 1) : (new_end)]
				new_repeat = Repeat(new_sequence, new_begin, new_end)
				cluster.modify_repeat(new_repeat, i)
				
				
	# compares nucleotides while remembering last actual state, i.e. number of mismatches is
	# number of state changes (ie. A A C C C has single mismatch - A->C transition
	# nucl column is a column of supposedly aligned nucleotides
	def count_mismatches(self, nucl_column):
		last_state = ''
		state_counts = {'A':0, 'C':0, 'G':0, 'T':0}
					
		num_mismatch = 0
		for nucl in nucl_column:
			if nucl.upper() in state_counts:
				state_counts[nucl.upper()] += 1
				
			else:
				num_mismatch += 1
				continue
		
			if nucl != last_state and last_state != '':
				# mismatch is only counted if change between states is introduced
				num_mismatch += 1
				
			last_state = nucl
		
		count = [state_counts['A'], state_counts['C'],state_counts['G'], state_counts['T']]
		
		#percentage = max(count)/len(nucl_column) # max(count)/sum(count)
		max_count = max(count) # counts for most frequent nucleotide
		
		return (num_mismatch, max_count)
				
				

	# uses specific pattern to precise precursor repeats from first pass, regroup them into clusters, then align them	
	def second_pass (self, fasta):
		# get clusters
		for repeat in self.repeat_list_first_pass:
			self.get_pattern_counts(repeat)
			self.extract_clusters(repeat)
					
		# align and extend DRs
		for cluster in reversed(self.clusters):
			seed_pos = self.align_seeds (cluster)
			self.extend (cluster, seed_pos, fasta)
			

	# verifies parameter for each cluster (length of DR, length of spacer)
	# this will also rewrite corect clusters as CRISPRs
	def validate (self, cluster, fasta):
		DR_length = len(cluster.repeats[1])

		# if repeat isn't long enough or is way too long
		if len(cluster.repeats[0]) < self.min_DR or len(cluster.repeats[0]) > self.max_DR :
			self.clusters.remove(cluster)
		else:

			# we the look on the length of spacers 
			if len(cluster) == 2:
				spacer_length = cluster.repeats[1].begin - cluster.repeats[0].end - 1
				if spacer_length < DR_length * self.min_spacer_DR_ratio or spacer_length > DR_length * self.max_spacer_DR_ratio :
					self.clusters.remove(cluster)
				else: 
					new_crispr = CRISPR(cluster.repeats, fasta)
					self.CRISPRs.append(new_crispr)	
			else:
				starting_DR = 0
				for i in range(0, len(cluster) - 1):
					spacer_length = cluster.repeats[i + 1].begin - cluster.repeats[i].end
					if spacer_length < (DR_length * self.min_spacer_DR_ratio) or spacer_length > (DR_length * self.max_spacer_DR_ratio) :
						# exporting repeat as CRISPR if long enough
						if i - starting_DR > 0:
							new_crispr = CRISPR(cluster.repeats[starting_DR : (i+1)], fasta)
							self.CRISPRs.append(new_crispr)		
						starting_DR = i + 1
				
				if len(cluster) - 1 - starting_DR > 0:
					new_crispr = CRISPR(cluster.repeats[starting_DR : ], fasta)
					self.CRISPRs.append(new_crispr)
					
	
	# searches for repeats before CRISPR which could be missed due to 
	# identicity of spacers
	def look_before_CRISPR(self, fasta, CRISPR1, CRISPR2 = None):
		startpos = CRISPR1.begin - len(CRISPR1.DR_consensus)*\
		 self.min_spacer_DR_ratio
		endpos = CRISPR1.begin - len(CRISPR1.DR_consensus)*\
		 self.max_spacer_DR_ratio
		self.search_DR_before(fasta, CRISPR1, CRISPR2, startpos, endpos)
					
				
	# search recursively for consensus after CRISPR
	def search_DR_before(self, fasta, CRISPR1, CRISPR2, startpos, endpos):
		startpos = int(startpos)
		endpos = int(endpos)
		interspace_filled = False
		last_DR_begin = CRISPR1.begin
		
		# if at start of sequence, quit
		if endpos < 0:
			return
				
		# okay if space between two consecutive CRISPRs is smaller
		if CRISPR2 is not None:
			if (last_DR_begin - CRISPR2.end) <= len(CRISPR1.DR_consensus) *\
			 self.max_spacer_DR_ratio:
				if(self.compare_DRs(CRISPR1.DR_consensus, CRISPR2.DR_consensus, 2)):
					for DR in reversed(CRISPR2.DR):
					# fuse this and next spacer					
						CRISPR1.insert_DR(DR, fasta, 0)

					self.CRISPRs.remove(CRISPR2)
				
				return
		
		# check for next DR
		i = startpos
		while i>endpos:
			if self.compare_DRs(fasta.sequence[(i-len(CRISPR1.DR_consensus)):i],\
			 CRISPR1.DR_consensus, 2):
				coords = [(i-len(CRISPR1.DR_consensus)), i]
				newDR = Repeat(fasta.sequence[coords[0]:coords[1]],\
				 (coords[0]+1), (coords[1]+1))
				startpos = i - len(CRISPR1.DR_consensus)*(1+self.min_spacer_DR_ratio)
				endpos = i - len(CRISPR1.DR_consensus)*(1+self.max_spacer_DR_ratio)
				CRISPR1.insert_DR(newDR, fasta, 0)
				self.search_DR_before(fasta, CRISPR1, CRISPR2, startpos, endpos)
				break
			
			i -= 1
			
	
	# searches for repeats after CRISPR which could be missed due to 
	# identicity of spacers
	def look_after_CRISPR(self, fasta, CRISPR1, CRISPR2 = None):
		startpos = CRISPR1.end + len(CRISPR1.DR_consensus)*\
		 self.min_spacer_DR_ratio
		endpos = CRISPR1.end + len(CRISPR1.DR_consensus)*\
		 self.max_spacer_DR_ratio
		 
		self.search_DR_after(fasta, CRISPR1,\
		 CRISPR2, startpos, endpos)
					
				
	# search recursively for consensus after CRISPR
	def search_DR_after(self, fasta, CRISPR1, CRISPR2, startpos, endpos):
		startpos = int(startpos)
		endpos = int(endpos)
		
		# if at end of sequence, quit
		if endpos > len(fasta.sequence):
			return
		
		interspace_filled = False
		last_DR_end = CRISPR1.end
				
		# okay if space between two consecutive CRISPRs is smaller
		if CRISPR2 is not None:
			if (CRISPR2.begin - last_DR_end) <= len(CRISPR1.DR_consensus) *\
			 self.max_spacer_DR_ratio:
				if(self.compare_DRs(CRISPR1.DR_consensus, CRISPR2.DR_consensus, 2)):
					for DR in CRISPR2.DR:
					# fuse this and next spacer					
						CRISPR1.insert_DR(DR, fasta)
						
					self.CRISPRs.remove(CRISPR2)
				
				return
		
		# check for next DR
		i = startpos
		while i<endpos:
			if self.compare_DRs(fasta.sequence[i:(i+len(CRISPR1.DR_consensus))],\
			 CRISPR1.DR_consensus, 2):
				coords = [i,(i+len(CRISPR1.DR_consensus))]
				newDR = Repeat(fasta.sequence[coords[0]:coords[1]],\
				 (coords[0]+1), (coords[1]+1))
				startpos = i + len(CRISPR1.DR_consensus)*(1+self.min_spacer_DR_ratio)
				endpos = i + len(CRISPR1.DR_consensus)*(1+self.max_spacer_DR_ratio)
				CRISPR1.insert_DR(newDR, fasta)
				self.search_DR_after(fasta, CRISPR1, CRISPR2, startpos, endpos)
				break
			
			i += 1
			
			
	# compares two sequences with certain tolerance	
	def	compare_DRs(self, sequence1, sequence2, tolerance):
		errors = 0
		if len(sequence1) == len(sequence2):
			for i in range(len(sequence1)):
				
				if sequence1[i] != sequence2[i]:
					errors += 1
			
				if errors > tolerance:
					return False	
					
		else:
			return False
	
		return True
					

	# uses parts of first or last DR of crispr to look at specified distance before or after CRISPR
	# if one parts align then comparison site by site is realized and if similarity is greater then 50%, first/last DR is added
	def check_truncated (self, crispr, fasta, direction):
		DR_length = len(crispr.DR[0])
		fragment_len = int((DR_length)/4)
		switch_limit = 0.33			# number of times truncated DR can switch between match and mismatch
		DR_found = False
		DR = ""
		DR_start = None
		DR_end = None	
		similarity = None				# percentage of similarity between first - last 
		last_position_status = None		# 0 = mismatch, 1 = identity
		identity_switch = None			# nnumber of alternations between identities and mismatches

		# algorithm description: the DR is cut into three or four separate pieces. We search each piece in sequence 
		# until a match is found. When a match is detected, we compare the rest of sequence base after base then we 
		# determine fraction of part identic to DR (and if that's >limit, we extend sequence by this DR).

		if direction == -1:

			search_range_start = crispr.DR[0].begin 
			search_range_end = crispr.DR[0].begin - DR_length * (self.max_spacer_DR_ratio + 1)
			if search_range_end < 0:
				search_range_end = 0
				
			fragment = 0
			while (fragment + fragment_len) < len(crispr.DR[0]) and DR_found == False:
				pos = crispr.DR[0].begin - 1
				while pos > search_range_end and DR_found == False:
					DR_start = pos - fragment_len - fragment + 1
					DR_end = pos - fragment_len - fragment + len(crispr.DR[0])
					lower_limit = crispr.DR[0].begin - DR_length * self.max_spacer_DR_ratio
					upper_limit = crispr.DR[0].begin - DR_length * self.min_spacer_DR_ratio
					if (crispr.DR[0].sequence[fragment:(fragment+fragment_len)] ==\
					 fasta.sequence[(pos - fragment_len) : pos]) and DR_end > lower_limit and DR_end < upper_limit:
						# checking if supposed start of DR is starting in sequence
						if (pos - fragment) > 0 :
							DR_found = True
							DR = ""
							similarity = 0.0
							identity_switch = 0.0 
							for i in range(0, len(crispr.DR[0])):
								# calculating similarity and number of changes between mismatch - identity and vice versa status
								if crispr.DR[0].sequence[i] == fasta.sequence[pos - fragment_len - fragment + i]:
									similarity += 1
									if last_position_status == 0:
										identity_switch += 1
										
									last_position_status = 1
								else:
									if last_position_status == 1:
										identity_switch += 1
										
									last_position_status = 0		
									
								DR+=fasta.sequence[pos - fragment_len - fragment + i]
								
							similarity /= len(crispr.DR[0])	
							identity_switch /= (len(crispr.DR[0]) - 1)
							if similarity <= self.identity_truncated and identity_switch >= switch_limit: 	# continue searching if conditions aren't respected
								DR_found = False
											
					pos -= 1	
				fragment += fragment_len			 

		elif direction == 1:

			search_range_start = crispr.DR[-1].end
			search_range_end = crispr.DR[-1].end + DR_length * (self.max_spacer_DR_ratio + 1)
			if search_range_end > len(fasta.sequence):
				search_range_end = 0
			fragment = 0
			while (fragment + fragment_len) < len(crispr.DR[0]) and DR_found == False:
				pos = crispr.DR[-1].end
				while pos < search_range_end and DR_found == False:
					DR_start = pos - fragment + 1
					DR_end = pos - fragment + len(crispr.DR[-1])
					lower_limit = crispr.DR[-1].end + DR_length * self.min_spacer_DR_ratio
					upper_limit = crispr.DR[-1].end + DR_length * self.max_spacer_DR_ratio
					if (crispr.DR[-1].sequence[fragment:(fragment+fragment_len)] ==\
					 fasta.sequence[pos : (pos + fragment_len)]) and DR_start > lower_limit and DR_start < upper_limit:
						if (pos + len(crispr.DR[-1]) - fragment) < len(fasta.sequence):
							DR_found = True
							DR = ""
							similarity = 0.0
							identity_switch = 0.0
							for i in range(0, len(crispr.DR[-1])):
								if crispr.DR[-1].sequence[i] == fasta.sequence[pos - fragment + i]:
									similarity += 1
									if last_position_status == 0:
										identity_switch += 1
									last_position_status = 1
								else:
									if last_position_status == 1:
										identity_switch += 1
									last_position_status = 0
								DR+=fasta.sequence[pos - fragment + i]
								
							similarity /= len(crispr.DR[-1])
							identity_switch /= len(crispr.DR[-1]) - 1
							if similarity <= self.identity_truncated and identity_switch >= switch_limit:		# continue searching if conditions aren't respected
								DR_found = False
					pos += 1
						
				fragment += fragment_len
				
		if DR_found and similarity > self.identity_truncated and identity_switch < switch_limit:
			if direction == 1:
				crispr.insert_DR(Repeat(DR, DR_start, DR_end), fasta)
				
			elif direction == -1:
				crispr.insert_DR(Repeat(DR, DR_start, DR_end), fasta, 0)

	
	# checks crispr for nucleotide composition and erases them if they present strong nucleotide bias (low complexity region)
	# returns 1 if crispr was not deleted and 0 otherwise, as to not check in next filter if it doesn't exist anymore
	def filter_low_complexity(self, crispr):
		if len(crispr) < 6:
			nucleotides = ['A', 'C', 'G', 'T']
			nucleotide_counts = [0.0, 0.0, 0.0, 0.0]	# counts for A, C, G and T
			
			# looking on composition by DR + spacer
			for i in range(0, 4):
				for j in range(0, len(crispr.DR) - 1):
					nucleotide_counts[i] += crispr.DR[j].sequence.count(nucleotides[i])
					nucleotide_counts[i] += crispr.spacers[j].sequence.count(nucleotides[i])
				
				nucleotide_counts[i] += crispr.DR[-1].sequence.count(nucleotides[i])

			nucleotide_total = crispr.end - crispr.begin
				
			for i in range(0,4):
				nucleotide_counts[i] /= nucleotide_total

			# comparison of composition A + T vs. C + G
			nucleotide_counts_AT = nucleotide_counts[0] + nucleotide_counts[3]

			if nucleotide_counts_AT < 0.2 or nucleotide_counts_AT > 0.8:			# strong bias to AT or GC in sequence
				self.CRISPRs.remove(crispr)	
				return 0
				
			nucleotide_counts.sort()

			# looking if least prevalent value is very rare or others are very prevalent
			if sum(nucleotide_counts[1:4]) > 0.95:
				self.CRISPRs.remove(crispr)		
				return 0
	
			return 1


	# filters out CRISPRs having too much common small sequences between DR and spacer
	def filter_tandem(self, crispr):
		match_counter = 0.0
		fractional_match = 0.0
		
		# we compare segments of spacer with both DRs if CRISPR has only 2 DRs (due to truncanted DRs)
		if len(crispr) == 2:
			for i in range(0, 2):
				for j in range(0, len(crispr.DR[i].sequence) - self.k_mer_size_filter):
					for k in range(0, len(crispr.spacers[0].sequence) - self.k_mer_size_filter):
						if crispr.DR[i].sequence[j : (j+self.k_mer_size_filter)] == crispr.spacers[0].sequence[k : (k + self.k_mer_size_filter)]:
							match_counter += 1
			
			fractional_match = match_counter/2	
		
		# if multiple spacers are present
		else:		
			reference = crispr.DR[2].sequence	# we take 2nd repeat as reference (since first is often truncated)
			for i in range(0, len(reference) - self.k_mer_size_filter):
				small_seed = reference[i : (i + self.k_mer_size_filter)]		
				for j in range(0, len(crispr.spacers)):
					for k in range(0, len(crispr.spacers[j]) - self.k_mer_size_filter):
						if small_seed == crispr.spacers[j].sequence[k : (k + self.k_mer_size_filter)]:
							match_counter += 1

			fractional_match = match_counter/len(crispr.spacers)

		if fractional_match > self.spacer_dr_match_limit:
			self.CRISPRs.remove(crispr)
			
			
	# once truncated DRs were added, checks if CRISPRs with them added still respect imposed conditions
	# looks two last nucleotides of each repeat using 'lTrim' and 'rTrim' recursively while alterning 
	# between them
	def border_trimmer (self, crispr):
		border_limits = self.identity_limit # number of mismatches allowed per column on borders
		if len(crispr) <= 2:
			border_limits = 0
		repeat_len = len(crispr.DR[0])
		checker_depth = 2
		if(repeat_len > self.min_DR):
			# start by right side 
			self.lTrim(crispr, True, checker_depth, border_limits)
	
					
	# trims crispr with too many mismatches on right side. Used by 'border_trimmer'.				
	def lTrim (self, crispr, rTrim_continue, checker_depth, border_limits):
		pos = 0
		lTrim_continue = False
		while (pos < checker_depth):
			border_nucls = [x.sequence[pos] for x in crispr.DR]
			# get most frequent nucl
			max_freq = max(set(border_nucls), key = border_nucls.count)
			if (len(border_nucls) - border_nucls.count(max_freq)) > border_limits:
				crispr.DR_consensus = crispr.DR_consensus[(pos+1):]
				for i in range(len(crispr.DR)):
					rej_nucls = crispr.DR[i].sequence[:(pos+1)] # nucls deplaced to spacer
					crispr.DR[i].sequence = crispr.DR[i].sequence[(pos+1):]
					crispr.DR[i].begin +=  pos + 1
					if(i>0):
						crispr.spacers[i-1].sequence += rej_nucls
						crispr.spacers[i-1].end += pos + 1
						
				lTrim_continue = True
								
			pos += 1
		
		# if previously there was something eliminated on right side of DR
		if rTrim_continue and len(crispr.DR[0]) > self.min_DR:
			self.rTrim(crispr, lTrim_continue, checker_depth, border_limits)
		
		# else try the same thing from the other side
		elif lTrim_continue and len(crispr.DR[0]) > self.min_DR: 
			self.lTrim(crispr, rTrim_continue, checker_depth, border_limits)	
					
	
	# trims crispr with too many mismatches on left side. Used by 'border_trimmer'.
	def rTrim (self, crispr, lTrim_continue, checker_depth, border_limits):		
		pos = 0
		rTrim_continue = False
		while (pos < checker_depth):
			border_nucls = [x.sequence[-1-pos] for x in crispr.DR]
			# get most frequent nucl
			max_freq = max(set(border_nucls), key = border_nucls.count)
			if (len(border_nucls) - border_nucls.count(max_freq)) > border_limits:
				crispr.DR_consensus = crispr.DR_consensus[:(-1-pos)]
				for i in range(len(crispr.DR)):
					rej_nucls = crispr.DR[i].sequence[(-1-pos):] # nucls deplaced to spacer
					crispr.DR[i].sequence = crispr.DR[i].sequence[:(-1-pos)]
					crispr.DR[i].end -=  (pos + 1)
					if(i < (len(crispr.DR) - 1)):
						crispr.spacers[i].sequence = rej_nucls +\
						 crispr.spacers[i].sequence
						crispr.spacers[i].begin -= (pos + 1)
						
				rTrim_continue = True
								
			pos += 1
		
		# if previously there was something eliminated on left side of DR
		if lTrim_continue and len(crispr.DR[0]) > self.min_DR: 
			self.lTrim(crispr, rTrim_continue, checker_depth, border_limits)
		
		# else try the same thing from the other side	
		elif rTrim_continue and len(crispr.DR[0]) > self.min_DR:
			self.rTrim(crispr, lTrim_continue, checker_depth, border_limits)
			
	
		# function detecting tracrRNA
	def detect_tracrRNA(self, fasta):
		# number of fragments to split DR to. NOTE: less fragments - lower 
		# precision, but (way) faster. Should be 2 or 3.
		DR_splits = 2 
		
		
		DR_consensuses = self.list_consensuses()	# list of DRs
		copyseq = fasta.sequence	# copy of sequence, since masks will be applied to it
		
		# apply masks to positions of CRISPRs
		for CRISPR in self.CRISPRs:
			CRISPR_len = CRISPR.end - CRISPR.begin + 1
			copyseq = copyseq[0:(CRISPR.begin-2)]+'-'*CRISPR_len+copyseq[CRISPR.end-1:]
	
		for consensus in DR_consensuses:
			len_segment = int(len(consensus)/DR_splits)
			for i in range(0, DR_splits):
				found_at_pos = 0
				while found_at_pos != -1:
					search_obj = re.search(consensus[i*len_segment:(i+1)*len_segment],\
					 copyseq[found_at_pos:])
					if search_obj is not None:
						found_at_pos += search_obj.span(0)[1]
						# compute start and end of supposed match
						match_start = found_at_pos - (len_segment*(i+1))
						match_end = match_start + len(consensus)
						if(self.compare_DRs(consensus,\
						 copyseq[match_start:match_end], 2)):
							print('tracrRNA detected at %s - %s for '\
							 %(match_start + 1, match_end)+consensus+'.')
							# mask found sequence
							copyseq = copyseq[0:match_start-1]+'-'*(len(consensus)+1)+\
							 copyseq[match_end:]
							
					else:
						found_at_pos = -1		
			

	# constructs a a list of unique DR consensuses
	# index indicates the number of fasta file (order of treatment) 
	def list_consensuses(self):
		# stores DR for each consensus
		DR_consensuses = []
		for CRISPR in self.CRISPRs:
			DR_consensuses.append(CRISPR.DR_consensus)
			DR_consensuses.append(self.complementary(CRISPR.DR_consensus))
		
		# check for consensus similarity and remove duplicates
		DR_consensuses = list(set(DR_consensuses))
		return DR_consensuses
	
	# creates complementary traduction of given sequence
	def complementary(self, sequence):
		sequence = sequence.upper()
		comp_bases = {'A':'T', 'C':'G', 'G':'C', 'T':'A'} 
		complement = ''
		for base in sequence:
			if base in comp_bases:
				complement += comp_bases[base]
			
			else:
				complement += base
				
		return complement
	

	# creates a driectory with output files in it	
	def output (self, fasta):
		header_split = ''
		version = ''
		try:
			header_split = fasta.header.split("|")		# splits header by "|"
			version = header_split[3].split(".")		# splits number containing CRISPR on RefSeq and version
			
		except IndexError:
			header_split = ['N/A']*5
			version = 'Unknown'
		RefSeq = version[0]
		seq_length = len(fasta.sequence)
		if self.output_path[-1] != "/":
			self.output_path += "/"	
			
		self.output_path += RefSeq + "/"	
		if not os.path.exists(self.output_path):	# create folder if not present already	
			os.makedirs(self.output_path)
		
		CRISPR_number = 0
		hyp_count = 0					# < 3 spacers (or < 4 DRs)
		one_spacer_count = 0				# = 1 spacer (or exactly 2 DRs)
		for CRISPR in self.CRISPRs:
			CRISPR_number += 1	
			if CRISPR.hypothetic == "yes":
				hyp_count += 1
				
			if len(CRISPR.spacers) == 1:
				one_spacer_count += 1

			self.generate_CRISPR_info(CRISPR, RefSeq, CRISPR_number, header_split, seq_length)
			self.generate_spacer_file(CRISPR, CRISPR_number)

		self.generate_main_file(fasta, RefSeq, CRISPR_number, one_spacer_count)	
		

	# generates an info file proper to each crispr:
	#	-	crispr is a concerned crispr
	#	-	RefSes is refseq number of corresponding sequence
	#	-	CRISPR_number is a number of CRISPR
	#	-	header_split is fasta header split by "|"
	#	-	seq_length is a length of sequence
	
	def generate_CRISPR_info(self, crispr, RefSeq, CRISPR_number, header_split, seq_length):
		file_path = self.output_path
		# name is different for hypothetic and "true" CRISPRs
		if crispr.hypothetic == "yes":
			file_path += RefSeq + "_PossibleCrispr_%d"%CRISPR_number
		else:
			file_path += RefSeq + "_Crispr_%d"%CRISPR_number

		actual_time = datetime.datetime.now()			

		info_output = "########################################\n"
		info_output += "# Program: Crispr Finder Program\n"
		info_output += "# Author: Juraj MICHALIK (original version by Ibtissem GRISSA)\n"
		info_output += "# Rundate: %02d/%02d/%04d %02d:%02d:%02d\n"%(actual_time.day, actual_time.month, actual_time.year, actual_time.hour, actual_time.minute, actual_time.second)
		info_output += "# Report_file: " + file_path + "\n"
		info_output += "########################################\n"
		info_output += "#=======================================\n"
		info_output += "#\n"
		info_output += "# Sequence: " + RefSeq + "\n"
		info_output += "# Description: " + header_split[4] + "\n"
		info_output += "# Length: %d\n"%(seq_length)
		info_output += "# Id: " + header_split[0][1:] + "|" + header_split[1] + "|" + header_split[2] + "|" + header_split[3] + "|\n"
		info_output += "#\n"
		info_output += "#=========================================================================\n"
		info_output += "# Crispr Rank in the sequence: %d\n"%CRISPR_number
		info_output += "# Crispr_begin_position: %d\t Crispr_end_position: %d\n"%(crispr.DR[0].begin, crispr.DR[-1].end)
		info_output += "# DR: " + crispr.DR_consensus +"\t DR_length: %d\t Number_of_spacers: %d\n"%(len(crispr.DR_consensus), len(crispr.spacers))
		info_output += "#=========================================================================\n"
		info_output += "Spacer_begin_position\t Spacer_length\t Spacer_sequence\n"
		for spacer in crispr.spacers:
			position_offset = " " * (21 - len(str(spacer.begin)) ) # this and next line are just for formatting purpose
			length_offset = " " * (13 - len(str(len(spacer))) )
			info_output += position_offset + "%d\t "%spacer.begin + length_offset + "%d\t "%len(spacer) + spacer.sequence + "\n" 
			info_output += "\n" 
		info_output += "#=========================================================================\n"
		info_output += "########################################"
		
		new_file = open(file_path, 'w')
		new_file.write(info_output)
		new_file.close()		
	
	# generates a fasta of spacers for each CRISPR
	#	-	crispr is corresponding crispr
	#	-	CRISPR_number is number of crispr
	
	def generate_spacer_file(self, crispr, CRISPR_number):
		file_path = self.output_path + "Spacers_%d"%CRISPR_number

		actual_time = datetime.datetime.now()

		spacer_output = ""		
		spacer_num = 0	
	
		for spacer in crispr.spacers:
			spacer_num += 1
			spacer_output += ">spacer%d\n"%CRISPR_number
			spacer_output += spacer.sequence + "\n"

		new_file = open(file_path, 'w')
		new_file.write(spacer_output)
		new_file.close()
		

	# generates a crispr summary for given sequence
	#	-	fasta is .fasta file of analyzed sequence
	#	-	RefSeq is its RefSeq
	#	-	crispr_total is a totla number of CRISPRs
	#	-	one_spacer_count is a number of CRISPRs with single spacer	

	def generate_main_file(self, fasta, RefSeq, crispr_total, one_spacer_count):
		file_path = self.output_path + RefSeq + "_CRSIPRs"
		actual_time = datetime.datetime.now()	
			
		main_output = "########################################\n"
		main_output += "# Program: Crispr Finder Program\n"
		main_output += "# Author: Juraj MICHALIK (original version by Ibtissem GRISSA)\n"
		main_output += "# Rundate: %02d/%02d/%04d %02d:%02d:%02d\n"%(actual_time.day,\
		 actual_time.month, actual_time.year, actual_time.hour, actual_time.minute,\
		  actual_time.second)
		main_output += "# Report_file: " + file_path + "\n"
		main_output += "########################################\n"
		main_output += "#=======================================\n"
		main_output += "#\n"
		main_output += "# Sequence: " + RefSeq + "\n"
		main_output += "# Nbr_of_good_CRISPRs : %d\n"%crispr_total
		main_output += "# Nbr_of_one_Spacer_CRISPR : %d\n"%one_spacer_count
		main_output += "#\n"		
		main_output += "#=======================================\n"

		new_file = open(file_path, 'w')
		new_file.write(main_output)
		new_file.close()
		

	# makes analysis of every sequence
	def analyze (self):
		fasta_n = 0
		for fasta in self.fastas.fna_list:
			self.first_pass(fasta)
			self.clusters = []			# resets clusters
			self.second_pass(fasta)
			
			for cluster in reversed(self.clusters):
				self.validate(cluster, fasta)
				
			self.CRISPRs = self.CRISPRs[::-1]	# invert result table
			
			for CRISPR in reversed(self.CRISPRs):
				# if 1 then continue with next filter, otherwise source
				check_ok = self.filter_low_complexity(CRISPR)	
				if check_ok == 1:
					self.filter_tandem(CRISPR)
							
			for i in range((len(self.CRISPRs) - 1), -1, -1):
				previousCRISPR = None
				thisCRISPR = self.CRISPRs[i]
				nextCRISPR = None
				if i > 0:
					previousCRISPR = self.CRISPRs[i - 1]
				
				if i < (len(self.CRISPRs) - 1):
					nextCRISPR = self.CRISPRs[i + 1]
						
				self.look_after_CRISPR(fasta, thisCRISPR, nextCRISPR)
				self.look_before_CRISPR(fasta, thisCRISPR, previousCRISPR)
				
			for CRISPR in self.CRISPRs:
				self.check_truncated (CRISPR, fasta, -1)
				self.check_truncated (CRISPR, fasta, 1)
				
			# deleting entries, so reversing passage by list
			for CRISPR in self.CRISPRs:
				self.border_trimmer(CRISPR)	
					
			for CRISPR in self.CRISPRs:
				self.CRISPRs_all[fasta_n].append(CRISPR)
					
			if self.search_tracrRNA:		
				self.detect_tracrRNA(fasta)

			# outputting files
			self.output(fasta)
			
			for CRISPR in reversed(self.CRISPRs):
				self.CRISPRs.remove(CRISPR)
			
			fasta_n += 1
			
		for i in range(len(self.CRISPRs_all)):
			print('-------------')
			print(self.fastas.fna_list[i].header)
			print('-------------')
			for CRISPR in self.CRISPRs_all[i]:
				print(CRISPR)
			
			

# uses classes FindRepeats to analyse every sequence of FastaList
# FindRepeats results are then selected by their proximity and other properties before obtaining CRISPRs

##################
# Object Classes #
##################

# contains repeat :
#	-	sequence
# 	-	starting position
#	-	ending position	
class Repeat:
	def __init__ (self, sequence, begin, end):
		self.sequence = sequence
		self.begin = begin
		self.end = end
		
	def __str__(self):
		return "Sequence: " + self.sequence + " Begin: %d End: %d"%(self.begin, self.end)

	def __len__(self):
		return len(self.sequence)


# contains spacer, functionally identical to repeat :
#	-	sequence
# 	-	starting position
#	-	ending position	
class Spacer:
	def __init__ (self, sequence, begin, end):
		self.sequence = sequence
		self.begin = begin
		self.end = end
		
	def __str__(self):
		return "Sequence: " + self.sequence + " Begin: %d End: %d"%(self.begin, self.end)

	def __len__(self):
		return len(self.sequence)


# a class containing a repeats which are in same cluster:
# 	-	list of repeats in same cluster
#	-	cluster_id is each cluster unique identifier and counts number of created instances of class Cluster
#	-	most_repeated stocks pattern sequence with maximum number of occurences (possible seed), first repeat where it can be found, and its position within repeat 
class Cluster:
	counter = 1

	def __init__ (self):
		self.cluster_id = Cluster.counter
		Cluster.counter += 1
		self.repeats = []
		self.most_repeated = [None]*3

		# adds repeat to a list	
	def add_repeat(self, repeat, index = None):
		if index is None:
			index = len(self.repeats)
		self.repeats.insert(index, repeat)	

	def modify_repeat(self, repeat, repeat_id):
		self.repeats[repeat_id] = repeat

		# len function returns size of cluster
	def __len__ (self):
		return len(self.repeats)

	def __str__ (self):
		string = "Start of Cluster %d \n"%self.cluster_id
		for repeat in self.repeats:
			string = string + str(repeat) + "\n"
		string = string + "End of Cluster %d \n"%self.cluster_id 
		return string

	def __del__ (self):
		del self.cluster_id


# a class containing integrity of information about CRISPR
# 	-	self.DR contains list of DR, which are instances of Repeat class
#	-	self.spacers contains list of spacers, which are instances of Spacer class
class CRISPR:
	#CRISPR_number = 1

	def __init__ (self, repeats, fasta):
		#self.crispr_id = CRISPR_number
		#CRISPR.CRISPR_number += 1
		self.begin = 0
		self.end = 0
		self.sequence_name = ""
		self.DR = repeats
		self.DR_consensus = ""		# contains consensus of DR
		self.spacers = []
		self.hypothetic = ""		# yes if len(self.spacers) < 3, no otherwise
		
		for i in range(0, len(self.DR) - 1):
			spacer_begin = self.DR[i].end + 1
			spacer_end = self.DR[i + 1].begin - 1
			new_spacer = Spacer(fasta.sequence[(spacer_begin - 1):spacer_end], spacer_begin, spacer_end)
			self.spacers.append(new_spacer)
		
		# calculate DR consensus			
		self.calculate_consensus()

		# saving start of CRISPR and its end separately
		self.begin = self.DR[0].begin
		self.end = self.DR[-1].end

		if len(self.spacers) < 3:
			self.hypothetic = "yes"
		else:
			self.hypothetic = "no"
			
	# determines consensus		
	def calculate_consensus(self):
		self.DR_consensus = ''
		for i in range(0, len(self.DR[0])):
			alphabet = ['A', 'C', 'G', 'T', 'a', 'c', 'g', 't']	# alpabet used in DNA
			nucl_counts = [0, 0, 0, 0] 				# counts of A, C, G, T per position
			for j in range(0, len(self.DR)):
				k = 0
				letter_found = 0
				while letter_found == 0 and k<4:
					# if character equals capital or small letter of DNA alphabet	
					if(self.DR[j].sequence[i] == alphabet[k] or self.DR[j].sequence[i] == alphabet[k+4]):	
						nucl_counts[k] += 1 
						letter_found == 1
						
					k += 1
					
			self.DR_consensus += alphabet[nucl_counts.index(max(nucl_counts))]

					
	# with adding DR, we have to add also corresponding Spacer as well as to update
	def insert_DR(self, DR, fasta, index = None):
		if index is None:
			index = len(self.DR)
			
		self.DR.insert(index, DR)
		# adding DR to start
		if index == 0:
			new_spacer = Spacer(fasta.sequence[self.DR[0].end : (self.DR[1].begin - 1)],\
			 self.DR[0].end + 1, self.DR[1].begin - 1)
			self.spacers.insert(0, new_spacer)
			self.begin = self.DR[0].begin		# updating start limit
			
		# adding DR to end
		elif index == (len(self.DR) - 1) or index == None : 
			new_spacer = Spacer(fasta.sequence[self.DR[len(self.DR) - 2].end\
			 : (self.DR[len(self.DR) - 1].begin - 1)],\
			 self.DR[len(self.DR) - 2].end + 1, self.DR[len(self.DR) - 1].begin - 1)
			self.spacers.insert(len(self.spacers), new_spacer)
			self.end = self.DR[-1].end		# updating end limit
		
		# update DR consensus
		self.calculate_consensus()
			

	def __str__ (self):
		string = ""	
		for i in range(0, len(self.DR) - 1):
			string += "DR [%d - %d] : "%(self.DR[i].begin,  self.DR[i].end) + self.DR[i].sequence + " | Spacer [%d - %d] : "%(self.spacers[i].begin,  self.spacers[i].end) + self.spacers[i].sequence + "\n"
		string += "DR [%d - %d] : "%(self.DR[-1].begin,  self.DR[-1].end) + self.DR[-1].sequence + "\n"
		return string

	def __len__ (self):
		return len(self.DR)
		

# a class containing pattern and functions of conversion pattern->string and string->pattern
#	-	self.pattern contains a string like ###_###, with # indicating a character to take into account and _ one to ignore		

class Pattern:
	def __init__ (self, pattern):
		self.pattern = pattern
		self.length = 0
		self.n_symbols = len(self.pattern)
		for symbol in self.pattern:
			if symbol == '#':
				self.length += 1
				
	def extract_pattern (self, input_string, position):
		output_string = ""
		actual_pos = 0
		# in case of end of string
		upper_limit = len(self.pattern)
		if (position + len(self.pattern)) > len(input_string):
			upper_limit = len(input_string) - position
		while actual_pos < upper_limit:
			if self.pattern[actual_pos] == '#':
				output_string = output_string + input_string[actual_pos+position]
			actual_pos = actual_pos + 1
		return output_string
	
	# returns positions of "_" within pattern
	def return_blanks (self):
		list_of_blanks = []
		i = 0
		for symbol in self.pattern:
			if symbol == '_':
				list_of_blanks.append(i)	
			i += 1
		return list_of_blanks		

	def __len__ (self):
		return self.length

if __name__ == '__main__':
	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument('--fasta', type=str, required=True)
	parser.add_argument('--output_dir', type=str, required=True)
	parser.add_argument('--kmer_size_filter', type=int, help='is used to compare segments of DRs and spacers', default=3)
	parser.add_argument('--pattern', type=str, help="is a combination of '#' (character) and '_' (space) used as seeds during second pass. First pass used continuous seed of length k_mer_size (number of '#')", default="####_####")
	parser.add_argument('--window_size', type=int, help="defines size of window in which k_mers will be searched for", default=200)
	parser.add_argument('--allowed_mismatch', type=int, help="defines a number of mismatch allowed in repeat", default=1)
	parser.add_argument('--spacer_dr_match_limit', type=int, help="is a maximum number of matches of length k_mer_size_filter per spacer between DR and spacer", default=20)
	parser.add_argument('--min_dr', type=int, help="defines minimum length in bp of repeat part", default=23)
	parser.add_argument('--max_dr', type=int, help="defines maximum length in bp of repeat part", default=55)
	parser.add_argument('--min_spacer_dr_ratio', type=float, help="is minimum quotient allowed length of spacer versus DR", default=0.6)
	parser.add_argument('--max_spacer_dr_ratio', type=float, help="is maximum quotient allowed length of spacer versus DR", default=2.5)
	parser.add_argument('--first_pass_limit', type=int, help="is maximum allowed distance between two regions with repeats", default=200)
	parser.add_argument('--search_tracrrna', action='store_true', default=False)

	args = vars(parser.parse_args())
	args = [args['fasta'],
			args['output_dir'],
			args['kmer_size_filter'],
			args['pattern'],
			args['window_size'],
			args['allowed_mismatch'],
			args['spacer_dr_match_limit'],
			args['min_dr'],
			args['max_dr'],
			args['min_spacer_dr_ratio'],
			args['max_spacer_dr_ratio'],
			args['first_pass_limit'],
			args['search_tracrrna']]

	findCRISPRs = FindCRISPRs(*args)
	findCRISPRs.analyze()


