#!/usr/bin/env python2.7

import os
import argparse
import gzip
from Bio import SeqIO

BASE_PAIRING = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}

###
# parse the input parameters
###
parser = argparse.ArgumentParser(description="SeqTailor DNA sequence extraction for genomic ranges in BED format")

parser.add_argument("-g", "--genome", help="genome sequence filename (FASTA / FASTA.GZ)")
parser.add_argument("-c", "--coordinate", choices=['0', '1'], default='1', help="coordinate indexing")
parser.add_argument("-s", "--strand", choices=['BOTH', 'FORWARD', 'REVERSE'], default='BOTH', help="strand")
parser.add_argument("-i", "--input", help="input filename (BED)")
parser.add_argument("-o", "--output", help="output filename. default: suffix (.DNA.fasta)")
parser.add_argument("-r", "--report", help="report filename. default: suffix (.report.txt)")

args = parser.parse_args()
GENOME = args.genome
COORDINATE = args.coordinate
STRAND = args.strand
INPUT = args.input
OUTPUT = args.output
REPORT = args.report

# generate the filenames for output and report, when these are not provided by the user
if not OUTPUT:
	OUTPUT = INPUT + '.DNA.fa'
if not REPORT:
	REPORT = INPUT + '.report.txt'

# check if the given output or report filenames already exist in the current folder
OUTPUT_exist = os.path.isfile('./'+OUTPUT)
REPORT_exist = os.path.isfile('./'+REPORT)
if OUTPUT_exist and REPORT_exist:
	print 'OUTPUT and REPORT filenames exist. Please try again using different file names.\n'
elif OUTPUT_exist:
	print 'OUTPUT filename exists. Please try again using different file names.\n'
elif REPORT_exist:
	print 'REPORT filename exists. Please try again using different file names.\n'
else:
	# proceeed when the given output or report filenames do not exist in the current folder

	# provide all defined parameters and filenames in the output file
	SeqTailor_PARAMETER = '# DNA sequence extraction for genomic ranges in BED format\n' \
						+ '# GENOME: ' + GENOME + '\n'\
						+ '# COORDINATE: ' + COORDINATE + '\n' \
						+ '# STRAND: ' + STRAND + '\n' \
						+ '# INPUT FILE: ' + INPUT + '\n' \
						+ '# OUTPUT FILE: ' + OUTPUT + '\n' \
						+ '# REPORT FILE: ' + REPORT + '\n\n'

	try:
		###
		# read input files, generate output files, identify chromosomes, and load sequences
		###
		if GENOME.endswith('fasta') or GENOME.endswith('fa'):
			FILE_GENOME = SeqIO.parse(open(GENOME), 'fasta')
		elif GENOME.endswith('gz'):
			GENOME_gz = gzip.open(GENOME, 'r')
			FILE_GENOME = SeqIO.parse(GENOME_gz, 'fasta')

		if INPUT.endswith('bed'):
			FILE_INPUT = open(INPUT, 'U')
		elif INPUT.endswith('bed.gz'):
			FILE_INPUT = gzip.open(INPUT, 'r')

		FILE_OUTPUT = open(OUTPUT, 'w')
		FILE_OUTPUT.write(SeqTailor_PARAMETER)
		FILE_REPORT = open(REPORT, 'w')
		FILE_REPORT.write(SeqTailor_PARAMETER)

		###
		# check through the input file, and generate a report for each input data
		# DNA sequence extraction will only be applied to the input data that passed this check
		###
		CHROM_set = set()
		data_all_list = []
		data_all_report_comment = dict()
		data_all_report_action = dict()
		data_all_report_submit = dict()
		data_range_set = set()
		data_range_final_list = []
		data_range_final_dict = dict()
		for eachline in FILE_INPUT:
			if not eachline.startswith('#') and eachline.strip():
				eachline = eachline.strip()
				data_all_list.append(eachline)
				field = eachline.split()
				try:
					# read the 3 mantatory fileds of data for each genomic ranges
					chrom = field[0]
					start = int(field[1])
					end = int(field[2])
					each_range = chrom + '_' + str(start) + '_' + str(end)

					# check each line of input, and recored the comment, the action, and the submitted data
					if each_range not in data_range_set:
						data_range_set.add(each_range)

						if start < 0 and end < 0:
							data_all_report_comment[eachline] = 'negative START and negative END'
							data_all_report_action[eachline] = 'skipped'
							data_all_report_submit[eachline] = 'none'

						elif end < 0:
							data_all_report_comment[eachline] = 'negative END'
							data_all_report_action[eachline] = 'skipped'
							data_all_report_submit[eachline] = 'none'

						elif start < 0:
							start = 0
							if end - start > 10000:
								each_range_new = chrom + '_' + str(start) + '_' + str(start+10000)
								data_all_report_comment[eachline] = 'negative START, and too long extraction range'
								data_all_report_action[eachline] = 'set START to 0, and cut sequence at 10,000bp from START'
								data_all_report_submit[eachline] = each_range_new
								CHROM_set.add(chrom)
								data_range_final_list.append(each_range)
								data_range_final_dict[each_range] = each_range_new
							else:
								each_range_new = chrom + '_' + str(start) + '_' + str(end)
								data_all_report_comment[eachline] = 'negative START'
								data_all_report_action[eachline] = 'set START to 0'
								data_all_report_submit[eachline] = each_range_new
								CHROM_set.add(chrom)
								data_range_final_list.append(each_range)
								data_range_final_dict[each_range] = each_range_new

						elif start > end:
							temp = start
							start = end
							end = temp
							if end - start > 10000:
								each_range_new = chrom + '_' + str(start) + '_' + str(start+10000)
								data_all_report_comment[eachline] = 'larger START than END, and too long extraction range'
								data_all_report_action[eachline] = 'swap START with END, and trim at 10,000bp from START'
								data_all_report_submit[eachline] = each_range_new
								CHROM_set.add(chrom)
								data_range_final_list.append(each_range)
								data_range_final_dict[each_range] = each_range_new
							else:
								each_range_new = chrom + '_' + str(start) + '_' + str(end)
								data_all_report_comment[eachline] = 'larger START than END'
								data_all_report_action[eachline] = 'swap START with END'
								data_all_report_submit[eachline] = each_range_new
								CHROM_set.add(chrom)
								data_range_final_list.append(each_range)
								data_range_final_dict[each_range] = each_range_new

						elif end - start > 10000:
							each_range_new = chrom + '_' + str(start) + '_' + str(start+10000)
							data_all_report_comment[eachline] = 'too long extraction range'
							data_all_report_action[eachline] = 'trim at 10,000bp from START'
							data_all_report_submit[eachline] = each_range_new
							CHROM_set.add(chrom)
							data_range_final_list.append(each_range)
							data_range_final_dict[each_range] = each_range_new

						else:
							data_all_report_comment[eachline] = 'normal'
							data_all_report_action[eachline] = 'normal'
							data_all_report_submit[eachline] = each_range
							CHROM_set.add(chrom)
							data_range_final_list.append(each_range)
							data_range_final_dict[each_range] = 'normal'

					else:
						data_all_report_comment[eachline] = 'duplicated'
						data_all_report_action[eachline] = 'skipped'
						data_all_report_submit[eachline] = 'none'

				except Exception:
					data_all_report_comment[eachline] = 'field error'
					data_all_report_action[eachline] = 'skipped'
					data_all_report_submit[eachline] = 'none'

		###
		# write the recorded comments, actions taken, submitted data, for all input data, into a report
		###
		for eachdata in data_all_list:
			comment = '\t> COMMENT: ' + data_all_report_comment[eachdata]
			action = '\t> ACTION: ' + data_all_report_action[eachdata]
			submit = '\t> SUBMITTED: ' + data_all_report_submit[eachdata]
			FILE_REPORT.write(eachdata+'\n'+comment+'\n'+action+'\n'+submit+'\n\n')

		###
		# read the genome sequence, and load the sequences of the given chromosomes
		###
		CHROM_SEQ_dict = dict()
		for each_genome in FILE_GENOME:
			chromosome, sequence = each_genome.id, each_genome.seq
			if chromosome in CHROM_set:
				if COORDINATE == '0':
					CHROM_SEQ_dict[chromosome] = sequence
				elif COORDINATE == '1':
					CHROM_SEQ_dict[chromosome] = 'N' + sequence

		###
		# extract DNA sequences for genomic ranges that have passed the checking step
		###
		for each_range in data_range_final_list:
			# when the submitted variant is exactly the input
			if data_range_final_dict[each_range] == 'normal':
				output_header_forward = '>' + each_range + '|+'
				output_header_reverse = '>' + each_range + '|-'
			# when the submitted variant is modified from the input
			else:
				each_range_new = data_range_final_dict[each_range]
				output_header_forward = '>' + each_range + ':' + each_range_new + '|+'
				output_header_reverse = '>' + each_range + ':' + each_range_new + '|-'
				each_range = each_range_new

			# obtain the mandatory information for genomic ranges
			range_field = each_range.split('_')
			chrom = range_field[0]
			start = int(range_field[1])
			end = int(range_field[2])
			chrom_sequence = CHROM_SEQ_dict[chrom]

			# extract the sequence within the genomic ranges from the forward strand
			# convert it to the reverse strand
			output_seq_forward = chrom_sequence[start: end + 1]
			output_seq_reverse = ''
			for i in range(0, len(output_seq_forward)):
				output_seq_reverse += BASE_PAIRING[output_seq_forward[len(output_seq_forward) - i - 1]]

			# wrap output sequences to 80nt per line, and
			output_seq_forward_wrap = ''
			output_seq_reverse_wrap = ''
			for i in range(1, len(output_seq_forward)+1):
				if i % 80 == 0:
					output_seq_forward_wrap += output_seq_forward[i-1] + '\n'
					output_seq_reverse_wrap += output_seq_reverse[i-1] + '\n'
				else:
					output_seq_forward_wrap += output_seq_forward[i-1]
					output_seq_reverse_wrap += output_seq_reverse[i-1]

			# write the output sequence to fasta file, according to the user's option on strand
			if STRAND == 'BOTH' or STRAND == 'FORWARD':
				FILE_OUTPUT.write(output_header_forward + '\n' + output_seq_forward_wrap + '\n')
			if STRAND == 'BOTH' or STRAND == 'REVERSE':
				FILE_OUTPUT.write(output_header_reverse + '\n' + output_seq_reverse_wrap + '\n')
			FILE_OUTPUT.write('\n')

	# if any IOerror, give an error message
	except IOError:
		FILE_OUTPUT = open(OUTPUT, 'w')
		FILE_REPORT = open(REPORT, 'w')
		FILE_OUTPUT.write('FILE OPEN ERROR\n')
		FILE_REPORT.write('FILE OPEN ERROR\n')
