#!/usr/bin/env python2.7
        
import os
import argparse
import gzip
from Bio import SeqIO

BASE_PAIRING = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}

###
# parse the input parameters
###
parser = argparse.ArgumentParser(description="SeqTailor DNA sequence extraction for genomic variants (neighbourhood) in VCF format")

parser.add_argument("-g", "--genome", help="genome sequence filename (FASTA / FASTA.GZ)")
parser.add_argument("-c", "--coordinate", choices=['0', '1'], default='1', help="coordinate indexing")
parser.add_argument("-s", "--strand", choices=['BOTH', 'FORWARD', 'REVERSE'], default='BOTH', help="strand")
parser.add_argument("-wd", "--window_down", type=int, default=25, help="window size downstream")
parser.add_argument("-wu", "--window_up", type=int, default=25, help="window size upstream")
parser.add_argument("-q", "--seq_type", choices=['BOTH', 'WT', 'MT'], default='BOTH', help="output sequence type")
parser.add_argument("-i", "--input", help="input filename (VCF)")
parser.add_argument("-o", "--output", help="output filename. default: suffix (.DNA.fasta)")
parser.add_argument("-r", "--report", help="report filename. default: suffix (.report.txt)")
parser.add_argument('--is_debug', default=False, type=lambda x: (str(x).lower() in ['true', '1', 'yes']),
        help="Debug (true) or production (false, default to production)")

args = parser.parse_args()
GENOME = args.genome
COORDINATE = args.coordinate
STRAND = args.strand
WINDOW_DOWN = args.window_down
WINDOW_UP = args.window_up
SEQUENCE_TYPE = args.seq_type
INPUT = args.input
OUTPUT = args.output
REPORT = args.report

# generate the filenames for output and report, when these are not provided by the user
if not OUTPUT:
    OUTPUT = INPUT + '.DNA.fa'
if not REPORT:
    REPORT = INPUT + '.report.txt'


# concatenate a list of nucleotides into sequence
def seqlist_to_seqstr(sequence_list):
    sequence_string = ''
    for each_seq in sequence_list:
        sequence_string += each_seq
    return sequence_string


# reverse the input sequence into complementary strand in reverse order
def seqreverse(sequence_input):
    sequence_length = len(sequence_input)
    sequence_reversed = ''
    for i in range(0, sequence_length):
        sequence_reversed += BASE_PAIRING[sequence_input[sequence_length - i - 1].upper()]
    return sequence_reversed

# multiplicative effect function:
# input: a queue of variants, current_index
# output: cartian products of all variants with the current sliding window
def output_seq_forward_cartesian_product_sets(current_data_variant_neigbors_list, current_data_variant_neigbors_dict):
    products = [ [] ]

    for each_variant in current_data_variant_neigbors_list:
        temp_products = []
        for product in products:
            for alt in current_data_variant_neigbors_dict[each_variant]: 
                temp_products.append(product + [alt])

        products = temp_products

    return products

## check if the given output or report filenames already exist in the current folder
#OUTPUT_exist = os.path.isfile('./'+OUTPUT)
#REPORT_exist = os.path.isfile('./'+REPORT)
#if OUTPUT_exist and REPORT_exist:
#    print ('OUTPUT and REPORT filenames exist. Please try again using different file names.\n')
#elif OUTPUT_exist:
#    print ('OUTPUT filename exists. Please try again using different file names.\n')
#elif REPORT_exist:
#    print ('REPORT filename exists. Please try again using different file names.\n')
#else:
#    # proceeed when the given output or report filenames do not exist in the current folder

def main():
    global WINDOW_DOWN, WINDOW_UP

    # when the window size (downstream or upstream) is negative, skip it
    if WINDOW_DOWN < 0 or WINDOW_UP < 0:
        FILE_OUTPUT = open(OUTPUT, 'w')
        FILE_REPORT = open(REPORT, 'w')
        FILE_OUTPUT.write('NEGATIVE WINDOW SIZE\n')
        FILE_REPORT.write('NEGATIVE WINDOW SIZE\n')

    else:
        # when any provied window size is greater than 5,000bp, trim it to 5,000bp
        if WINDOW_DOWN > 5000:
            WINDOW_DOWN = 5000
        if WINDOW_UP > 5000:
            WINDOW_UP = 5000

        # provide all defined parameters and filenames in the output file
        SeqTailor_PARAMETER = '# DNA sequence extraction for genomic variants (neighbourhood) in VCF format\n' \
                            + '# GENOME: ' + GENOME + '\n' \
                            + '# COORDINATE: ' + COORDINATE + '\n' \
                            + '# STRAND: ' + STRAND + '\n' \
                            + '# WINDOW_DOWN: ' + str(WINDOW_DOWN) + ' bp\n' \
                            + '# WINDOW_UP: ' + str(WINDOW_UP) + ' bp\n' \
                            + '# SEQUENCE_TYPE: ' + SEQUENCE_TYPE + '\n' \
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

            if INPUT.endswith('vcf'):
                FILE_INPUT = open(INPUT, 'r')
            elif INPUT.endswith('vcf.gz'):
                FILE_INPUT = gzip.open(INPUT, 'r')

            FILE_OUTPUT = open(OUTPUT, 'w')
            FILE_OUTPUT.write(SeqTailor_PARAMETER)
            FILE_REPORT = open(REPORT, 'w')
            FILE_REPORT.write(SeqTailor_PARAMETER)

            # initial scan of the chromosomes of the given variants
            CHROM_set = set()
            for eachline in FILE_INPUT:
                if not eachline.startswith('#') and eachline.strip():
                    field = eachline.strip().split()
                    chrom = field[0]
                    CHROM_set.add(chrom)

            # read the genome sequence, and load the sequences of the given chromosomes
            CHROM_SEQ_dict = dict()
            for each_genome in FILE_GENOME:
                chromosome, sequence = each_genome.id, each_genome.seq
                if chromosome in CHROM_set:
                    if COORDINATE == '0':
                        CHROM_SEQ_dict[chromosome] = sequence
                    elif COORDINATE == '1':
                        CHROM_SEQ_dict[chromosome] = 'N' + sequence

            ###
            # check through the input file, and generate a report for each input data
            # DNA sequence extraction will only be applied to the input data that passed this check
            ###
            data_all_list = []
            data_all_report_comment = dict()
            data_all_report_action = dict()
            data_all_report_submit = dict()
            data_variant_set = set()
            data_variant_final_list = []
            data_variant_neighbor_dict = dict()
            data_variant_final_dict = dict()
            data_variant_dict = dict()
            variant_to_line_dict = dict()

            # the first base pair of ref is the mid-point
            # the window size is fixed and the window is sliding forward.
            # Every time a variant's all combinatorial/multiplicative variants are computed,
            # the window move forward to focus on the next variant.
            WINDOW_SIZE = WINDOW_UP + 1 + WINDOW_DOWN
            
            # all neighors of the current variant within WINDOW_SIZE base pairs
            # key: pos, value: list of all alterative alleles 
            current_data_variant_pos = -1
            current_data_variant_index = -1
            window_pos_left = -1
            window_pos_right = -1
            current_data_variant_neigbors_list = []
            current_data_variant_neigbors_dict = dict()
            output_seq_position_debug_dict = dict()

            # variant combinations for the current variant
            output_seq_forward_sets_dict = dict()
            
            # enter each_variant into the queue, for each alt allele, create one variant.
            def enqueue_variant(each_variant):
                variant_field = each_variant.split('_')
                ref = variant_field[2]
                alt = variant_field[3]

                # treat ref as one of the variants and put it in the first place always
                alt =  ref + "," + alt
                if ',' in alt:
                    each_variant_new_list = []
                    alt_list = alt.split(',')

                    for each_alt in alt_list:
                        each_variant_temp = chrom + '_' + str(pos) + '_' + ref + '_' + each_alt
                        each_variant_new_list.append(each_variant_temp)

                    current_data_variant_neigbors_dict[each_variant] = each_variant_new_list
                                

            # the holder for the chr number for previous line
            previous_chrom = 'chr_invalid'
            
            # when cross_chrom is true, the current line is the last line in the current
            # chromosome. It is the time to empty out the queue.
            cross_chrom = False
            current_line_num = 0

            def get_last_valid_line_number(FILE_INPUT):
                line_num = 0
                try:
                    if FILE_INPUT != None:
                        FILE_INPUT.seek(0)
                        for eachline in FILE_INPUT:
                            if not eachline.startswith('#') and eachline.strip():
                                eachline = eachline.strip()
                                data_all_list.append(eachline)
                                field = eachline.split()

                                chrom = field[0]
                                pos = int(field[1])
                                ref = field[3]
                                alt = field[4]
                                
                                each_variant = chrom + '_' + str(pos) + '_' + ref + '_' + alt
                                variant_to_line_dict[each_variant] = eachline
                                
                                # check each line of input, and recored the comment, the action, and the submitted data
                                # dict may be faster than set: O(1) vs O(1) average or O(n) worst case.
                                # reference: https://www.geeksforgeeks.org/internal-working-of-set-in-python/
                                #if each_variant not in data_variant_set:
                                #    data_variant_set.add(each_variant)
                                if each_variant not in data_variant_dict:
                                    data_variant_dict[each_variant] = 1
        
                                    if ref == '.' or ref == '-':
                                        data_all_report_comment[eachline] = 'empty REF allele'
                                        data_all_report_action[eachline] = 'skipped'
                                        data_all_report_submit[eachline] = 'none'
        
                                    elif ref != CHROM_SEQ_dict[chrom][pos: pos+len(ref)].upper():
                                        data_all_report_comment[eachline] = 'unmatched REF allele to the reference genome'
                                        data_all_report_action[eachline] = 'skipped'
                                        data_all_report_submit[eachline] = 'none'
        
                                    elif ref == alt:
                                        data_all_report_comment[eachline] = 'identical REF and ALT allele'
                                        data_all_report_action[eachline] = 'skipped'
                                        data_all_report_submit[eachline] = 'none'
        
                                    elif ',' in alt:
                                        line_num += 1
                                        
                                        alt_list = alt.split(',')
                                        each_variant_new = ''
                                        for each_alt in alt_list:
                                            each_variant_temp = chrom + '_' + str(pos) + '_' + ref + '_' + each_alt
                                            each_variant_new += each_variant_temp + ','
                                            # data_variant_final_list.append(each_variant_temp)
                                            #data_variant_neighbor_dict[each_variant_temp] = set()
                                        each_variant_new = each_variant_new[0:-1]
                                        data_all_report_comment[eachline] = 'normal'
                                        data_all_report_action[eachline] = 'normal'
                                        data_all_report_submit[eachline] = each_variant_new
                                        data_variant_final_dict[each_variant] = 1
        
                                    else:
                                        line_num += 1
                                        
                                        data_all_report_comment[eachline] = 'normal'
                                        data_all_report_action[eachline] = 'normal'
                                        data_all_report_submit[eachline] = each_variant
                                        # data_variant_final_list.append(each_variant)
                                        #data_variant_neighbor_dict[each_variant] = set()
                                        data_variant_final_dict[each_variant] = 1
        
                                else:
                                    data_all_report_comment[eachline] = 'duplicated'
                                    data_all_report_action[eachline] = 'skipped'
                                    data_all_report_submit[eachline] = 'none'

                except IOError:
                    assert False, "Error happens in file I/O."
                except Exception:
                    data_all_report_comment[eachline] = 'field error'
                    data_all_report_action[eachline] = 'skipped'
                    data_all_report_submit[eachline] = 'none'

                return line_num
            
            total_line_num = get_last_valid_line_number(FILE_INPUT)

            # output to the file
            # generate variant combinations for one line (one position) each time
            def output_seq_to_file(current_variant):
                # indicate whether the combination is the first one to dump.
                first_combination_forward = True
                first_combination_reverse = True
                
                ###
                # extract DNA sequences for genomic variants that have passed the checking step
                ###
                # extract the mandatory information for each variant
                each_variant = current_variant
                variant_field = each_variant.split('_')
                chrom = variant_field[0]
                pos = int(variant_field[1])
                ref = variant_field[2]
                alt = variant_field[3]
                chrom_sequence = CHROM_SEQ_dict[chrom]

                # generate the header for each type of output sequence
                output_header_forward_ref = '>' + each_variant + '|+|ref'
                output_header_forward_alt = '>' + each_variant + '|+|alt'
                output_header_reverse_ref = '>' + each_variant + '|-|ref'
                output_header_reverse_alt = '>' + each_variant + '|-|alt'

                # positions on the forward strand for forward strand sequence extraction
                forward_end_up = pos - WINDOW_UP
                forward_end_down = pos + WINDOW_DOWN + 1
                pos_new_forward = WINDOW_UP
                pos_new_reverse = WINDOW_DOWN

                ###
                # extract the reference DNA sequence in forward strand, and convert it to reverse strand
                # write each nucleotide into a list of reference nucleotides
                ###
                output_seq_forward_ref = chrom_sequence[forward_end_up: forward_end_down]
                output_seq_forward_ref_list = list(output_seq_forward_ref)
                output_seq_reverse_ref = seqreverse(output_seq_forward_ref)
                output_seq_reverse_ref_list = list(output_seq_reverse_ref)

                # modify the list of reference nucleotides in forward strand and reverse strand
                # replicate list of ref nucleotides to list of alt nucleotides and empty the corresponding cell(s)
                output_seq_forward_alt_list = output_seq_forward_ref_list
                for temp_pos in range(pos_new_forward, pos_new_forward + len(ref)):
                    output_seq_forward_alt_list[temp_pos] = ''

                output_seq_reverse_alt_list = output_seq_reverse_ref_list
                for temp_pos in range(pos_new_reverse - len(ref) + 1, pos_new_reverse + 1):
                    output_seq_reverse_alt_list[temp_pos] = ''

                ## generatt alternative nucleotides in forward strand and reverse strand
                #if alt == '.' or alt == '-':
                #    # already emptyed
                #    pass
                #else:
                #    output_seq_forward_alt_list[pos_new_forward] = alt
                #    alt_reverse = seqreverse(alt)
                #    output_seq_reverse_alt_list[pos_new_reverse] = alt_reverse
                output_header_forward_ref_origin = output_header_forward_ref
                output_header_forward_alt_origin = output_header_forward_alt
                output_header_reverse_ref_origin = output_header_reverse_ref
                output_header_reverse_alt_origin = output_header_reverse_alt
                
                for variant_neighbor_list in output_seq_forward_sets_dict[current_variant]:
                    ###
                    # if there are neighbor variants detected within the window size of sequence extraction
                    ###
                    # extend the header for each type of output sequence
                    output_header_forward_ref = output_header_forward_ref_origin + '|Neighbor:'
                    output_header_forward_alt = output_header_forward_alt_origin + '|Neighbor:'
                    output_header_reverse_ref = output_header_reverse_ref_origin + '|Neighbor:'
                    output_header_reverse_alt = output_header_reverse_alt_origin + '|Neighbor:'

                    # position information for current variant combinations
                    if args.is_debug:
                        seq_details = dict()
                        seq_details['num_neighbors'] = len(variant_neighbor_list)
                        seq_details['center_index'] = -1
                        
                        # coordinates in chromesome: (start, end, length) for each neighbor
                        seq_details['neighbor_coordinates'] = list()
                        
                        # indices of each neighbor in the final sequence string: (start, end, length)
                        seq_details['seq_list_indices'] = list()
                    
                    # iteratively process the neighbor variants that fall inside the window
                    current_index = 0
                    for each_neighbor in variant_neighbor_list:
                        # append each neighbor variant to the header
                        # if each_neighbor != current_variant:
                        output_header_forward_ref += each_neighbor + ';'
                        output_header_forward_alt += each_neighbor + ';'
                        output_header_reverse_ref += each_neighbor + ';'
                        output_header_reverse_alt += each_neighbor + ';'

                        neighbor_field = each_neighbor.split('_')
                        neighbor_pos = int(neighbor_field[1])
                        neighbor_ref = neighbor_field[2]
                        neighbor_alt = neighbor_field[3]

                        # positioning the neighbor variant on the extracted forward strand
                        neighbor_pos_new_forward = WINDOW_UP + neighbor_pos - pos
                        neighbor_pos_new_reverse = WINDOW_DOWN - neighbor_pos + pos

                        # replicate list of ref nucleotides to list of alt nucleotides and empty the corresponding cell(s)
                        for neighbor_temp_pos in range(neighbor_pos_new_forward, neighbor_pos_new_forward + len(neighbor_ref)):
                            output_seq_forward_alt_list[neighbor_temp_pos] = ''

                        for neighbor_temp_pos in range(neighbor_pos_new_reverse - len(neighbor_ref) + 1, neighbor_pos_new_reverse + 1):
                            output_seq_reverse_alt_list[neighbor_temp_pos] = ''

                        # generate alternative nucleotides in forward strand and reverse strand
                        if neighbor_alt == '.' or neighbor_alt == '-':
                            # already emptyed
                            pass
                        else:
                            output_seq_forward_alt_list[neighbor_pos_new_forward] = neighbor_alt
                            neighbor_alt_reverse = seqreverse(neighbor_alt)
                            output_seq_reverse_alt_list[neighbor_pos_new_reverse] = neighbor_alt_reverse
                        
                        if args.is_debug:
                            if pos == neighbor_pos:
                                seq_details['center_index'] = current_index
                            
                            current_index += 1
                            
                            # position information
                            alt_length = len(neighbor_alt)
                            seq_details['neighbor_coordinates'].append((each_neighbor, neighbor_alt, alt_length, neighbor_pos, 
                                                                        neighbor_pos + alt_length -1))
                            seq_details['seq_list_indices'].append((neighbor_pos_new_forward, neighbor_pos_new_forward + alt_length-1))
                    
                    # remove the last comma in the appended header
                    output_header_forward_ref = output_header_forward_ref[0:-1]
                    output_header_forward_alt = output_header_forward_alt[0:-1]
                    output_header_reverse_ref = output_header_reverse_ref[0:-1]
                    output_header_reverse_alt = output_header_reverse_alt[0:-1]
    
                    ###
                    # concantenate the list of modified nucleotides, to generate the alternative sequences
                    ###
                    output_seq_forward_alt = seqlist_to_seqstr(output_seq_forward_alt_list)
                    output_seq_reverse_alt = seqlist_to_seqstr(output_seq_reverse_alt_list)
                    
                    if args.is_debug:
                        seq_details['overall'] = (forward_end_up, forward_end_down-1, len(output_seq_forward_alt_list))
    
                    ## debug information for foward sequence only
                    def dump_seq_forward_debug_info():
                        if args.is_debug:
                            center_index = seq_details['center_index']
                            FILE_OUTPUT.write(f"#overall: {seq_details['overall']}\n")
                            FILE_OUTPUT.write(f"#num of neighbors: {seq_details['num_neighbors']}\n")
                            FILE_OUTPUT.write(f"#center_index: {center_index}, list index: {seq_details['seq_list_indices'][center_index]}\n")
                            for i in range(0, seq_details['num_neighbors']): 
                                FILE_OUTPUT.write(f"#{i}-th neighbor coordinates: {seq_details['neighbor_coordinates'][i]}, ")
                                FILE_OUTPUT.write(f"#list indices: {seq_details['seq_list_indices'][i]}\n")
                                
                    ###
                    # process output sequences
                    # wrap the sequences to 80nt per line
                    ###
                    output_seq_forward_ref_wrap = ''
                    for i in range(1, len(output_seq_forward_ref)+1):
                        if i % 80 == 0:
                            output_seq_forward_ref_wrap += output_seq_forward_ref[i-1] + '\n'
                        else:
                            output_seq_forward_ref_wrap += output_seq_forward_ref[i-1]
    
                    output_seq_forward_alt_wrap = ''
                    for i in range(1, len(output_seq_forward_alt)+1):
                        if i % 80 == 0:
                            output_seq_forward_alt_wrap += output_seq_forward_alt[i-1] + '\n'
                        else:
                            output_seq_forward_alt_wrap += output_seq_forward_alt[i-1]
    
                    output_seq_reverse_ref_wrap = ''
                    for i in range(1, len(output_seq_reverse_ref)+1):
                        if i % 80 == 0:
                            output_seq_reverse_ref_wrap += output_seq_reverse_ref[i-1] + '\n'
                        else:
                            output_seq_reverse_ref_wrap += output_seq_reverse_ref[i-1]
    
                    output_seq_reverse_alt_wrap = ''
                    for i in range(1, len(output_seq_reverse_alt)+1):
                        if i % 80 == 0:
                            output_seq_reverse_alt_wrap += output_seq_reverse_alt[i-1] + '\n'
                        else:
                            output_seq_reverse_alt_wrap += output_seq_reverse_alt[i-1]
    
                    # write the output sequence to fasta file, according to the user's option on strand and sequence type
                    if STRAND == 'BOTH' or STRAND == 'FORWARD':
                        # only reference
                        if SEQUENCE_TYPE == 'WT':
                            FILE_OUTPUT.write(output_header_forward_ref + '\n')
                            dump_seq_forward_debug_info()
                            FILE_OUTPUT.write(output_seq_forward_ref_wrap + '\n')
                            FILE_OUTPUT.write('\n')
                            break
                        # the first variant combination is always 'WT'. 
                        # For 'MT' only, we need to skip the first combination.
                        elif SEQUENCE_TYPE == 'MT':
                            if first_combination_forward:
                                first_combination_forward = False
                            else:
                                FILE_OUTPUT.write(output_header_forward_alt + '\n')
                                dump_seq_forward_debug_info()
                                FILE_OUTPUT.write(output_seq_forward_alt_wrap + '\n\n')
                        elif SEQUENCE_TYPE == 'BOTH':
                            if first_combination_forward:
                                FILE_OUTPUT.write(output_header_forward_ref + '\n')
                                dump_seq_forward_debug_info()
                                FILE_OUTPUT.write(output_seq_forward_ref_wrap + '\n\n')
                                first_combination_forward = False
                            else:
                                FILE_OUTPUT.write(output_header_forward_alt + '\n')
                                dump_seq_forward_debug_info()
                                FILE_OUTPUT.write(output_seq_forward_alt_wrap + '\n\n')
    
                    if STRAND == 'BOTH' or STRAND == 'REVERSE':
                        if SEQUENCE_TYPE == 'WT':
                            FILE_OUTPUT.write(output_header_reverse_ref + '\n' + output_seq_reverse_ref_wrap + '\n\n')
                            break
                        # the first variant combination is always 'WT'. 
                        # For 'MT' only, we need to skip the first combination.
                        elif SEQUENCE_TYPE == 'MT':
                            if first_combination_reverse:
                                first_combination_reverse = False
                            else:
                                FILE_OUTPUT.write(output_header_reverse_alt + '\n' + output_seq_reverse_alt_wrap + '\n\n')
                        elif SEQUENCE_TYPE == 'BOTH':
                            if first_combination_reverse:
                                FILE_OUTPUT.write(output_header_reverse_ref + '\n' + output_seq_reverse_ref_wrap + '\n\n')
                                first_combination_reverse = False
                            else:
                                FILE_OUTPUT.write(output_header_reverse_alt + '\n' + output_seq_reverse_alt_wrap + '\n\n')
    
            

            ###
            # write the recorded comments, actions taken, submitted data, into a report
            ###
            def write_report_for_unprocessed_line(eachdata):
                comment = '\t> COMMENT: ' + data_all_report_comment[eachdata]
                action = '\t> ACTION: ' + data_all_report_action[eachdata]
                submit = '\t> SUBMITTED: ' + data_all_report_submit[eachdata]
                FILE_REPORT.write(eachdata+'\n'+comment+'\n'+action+'\n'+submit+'\n\n')
            
            def write_report_for_processed_line(current_variant, SEQUENCE_TYPE):
                eachdata = variant_to_line_dict[current_variant]
                comment = '\t> COMMENT: ' + data_all_report_comment[eachdata]
                action = '\t> ACTION: ' + data_all_report_action[eachdata]
                submit = '\t> SUBMITTED: ' + data_all_report_submit[eachdata]
                
                # Do not distinguish strand status
                if SEQUENCE_TYPE == 'WT':
                    neighbors = output_seq_forward_sets_dict[current_variant][0]
                    neighborhood = '\t> Neighbor: ' + ','.join(neighbors)
                    FILE_REPORT.write(current_variant+'\n'+comment+'\n'+action+'\n'+submit+'\n'+neighborhood+'\n\n')
                elif SEQUENCE_TYPE == 'MT':
                    for neighbors in output_seq_forward_sets_dict[current_variant][1:]:
                        neighborhood = '\t> Neighbor: ' + ','.join(neighbors)
                        FILE_REPORT.write(current_variant+'\n'+comment+'\n'+action+'\n'+submit+'\n'+neighborhood+'\n\n')
                        
                    FILE_REPORT.write('\n\n')
                elif SEQUENCE_TYPE == 'BOTH':
                    for neighbors in output_seq_forward_sets_dict[current_variant]:
                        neighborhood = '\t> Neighbor: ' + ','.join(neighbors)
                        FILE_REPORT.write(current_variant+'\n'+comment+'\n'+action+'\n'+submit+'\n'+neighborhood+'\n\n')
                        
                    FILE_REPORT.write('\n')
                   

            FILE_INPUT.seek(0)
            for eachline in FILE_INPUT:
                if not eachline.startswith('#') and eachline.strip():
                    eachline = eachline.strip()
                    data_all_list.append(eachline)
                    field = eachline.split()
                    try:
                        # read the 4 mantatory fields of data for each genomic variant
                        chrom = field[0]
                        pos = int(field[1])
                        ref = field[3]
                        alt = field[4]
                        each_variant = chrom + '_' + str(pos) + '_' + ref + '_' + alt
                        
                        ## if centered is true, use the center of the current variant as
                        # the center of the final output. This will only affect seq position
                        # w/ len(alt) > 1, specifically,
                        # 1. if centered and len(alt) is even, then the window_left = len(alt)/2

                        if each_variant in data_variant_final_dict:                        
                            # current line being processed in the current file
                            current_line_num += 1

                            if previous_chrom != chrom:
                                # the first line
                                if previous_chrom != 'chr_invalid':
                                    # time to move to a new chromosome
                                    # we need to empty out the queue, reset the window, and neigbor dict, etc.
                                    cross_chrom = True
                                else:
                                    cross_chrom = False
                                
                            # the last line
                            if current_line_num == total_line_num:
                                cross_chrom = True
                                
                            previous_chrom = chrom
                                
                            # if this variant is within the sliding window of the variant being processed,
                            # add it into the queue;
                            # otherwise, dequeue the first variant in the queue because it has finished processing.
                            if len(current_data_variant_neigbors_list) == 0:
                                # this is the first in the queue and is the current variant being processed,
                                # the invariants have to be updated.
                                current_data_variant_neigbors_list.append(each_variant)
                                current_data_variant_pos = pos
                                current_data_variant_index = 0
                                
                                # for negative position values, patch "N" as base pairs
                                window_pos_left = pos - WINDOW_UP

                                # position values great than current chrom size, patch "N" as base pairs
                                window_pos_right = pos + WINDOW_DOWN
                                
                                # put the variant into the dict window
                                enqueue_variant(each_variant)

                            # if the variant is within the current sliding window
                            elif pos >= window_pos_left and pos <= window_pos_right:
                                # the current line is not the current variant being processed,
                                # the invariants do not need to be updated here.
                                current_data_variant_neigbors_list.append(each_variant)
                                enqueue_variant(each_variant)
                            
                            # if the variant is to the right side of the current sliding window,
                            # it is time to process the variants in the queue.
                            # after processing, the sliding window and queue have to be updated.
                            else:
                                while pos > window_pos_right:
                                    # time to process the current variant: process and update the queue/window
                                    # iteratively until the current line can enter the queue.
                                    # cartesian products of all variant combinations for the current variant
                                    current_variant = current_data_variant_neigbors_list[current_data_variant_index]

                                    # output to the file
                                    output_seq_forward_sets_dict[current_variant] = output_seq_forward_cartesian_product_sets(current_data_variant_neigbors_list, 
                                                                                    current_data_variant_neigbors_dict)
                                    
                                    print(f"output the variant combinations for {current_variant}")
                                    output_seq_to_file(current_variant)
                                    write_report_for_processed_line(current_variant, SEQUENCE_TYPE)

                                    ## after processing, maintain the sliding window
                                    # there is at least one more unprocessed variant in the queue
                                    if current_data_variant_index+1 < len(current_data_variant_neigbors_list):
                                        # this variant is the current variant being processed, the invariants
                                        # need to be updated.
                                        # update the current variant index
                                        current_data_variant_index += 1

                                        # extract the mandatory information for current variant
                                        current_variant = current_data_variant_neigbors_list[current_data_variant_index]
                                        variant_field = current_variant.split('_')
                                        chrom = variant_field[0]

                                        # update sliding window
                                        current_data_variant_pos = int(variant_field[1])
                                        window_pos_left = current_data_variant_pos - WINDOW_UP
                                        window_pos_right = current_data_variant_pos + WINDOW_DOWN
                                        
                                        # update the queue: step 1
                                        # update the left-hand side of the variant queue: 
                                        # remove any variant moving out of the left of the new sliding window
                                        for i in range(0, len(current_data_variant_neigbors_list)):
                                            # extract the mandatory information for each variant
                                            variant = current_data_variant_neigbors_list[0]
                                            variant_field = variant.split('_')
                                            chrom = variant_field[0]
                                            temp_pos = int(variant_field[1])
                                            
                                            if temp_pos < window_pos_left:
                                                # fall to the left of the sliding window,
                                                # none of the variants in the queue would fall to the right,
                                                # because the updated window moves to the right
                                                current_data_variant_neigbors_list.pop(0)
                                                current_data_variant_index -= 1

                                                # update neigbor dict
                                                current_data_variant_neigbors_dict.pop(variant, None)
                                            else:
                                                break

                                        # update the queue: step 2
                                        # update the right-hand side of the variant queue: 
                                        # try to enqueue the current line of variant in I/O
                                        if pos >= window_pos_left and pos <= window_pos_right:
                                            # enqueue the current line
                                            current_data_variant_neigbors_list.append(each_variant)
                                            
                                            # update its neigbor variants
                                            # Note: since it is not the current variant being processed,
                                            # the invariants do not need to update. In fact, the invariants
                                            # have been updated above in this while loop.
                                            enqueue_variant(each_variant)

                                            # once the current line enters the queue, we are ready to proceed
                                            # to the next line since the invariants are maintained.
                                            # break out of the loop: while pos > window_pos_right
                                            break

                                    # the last variant in the queue has been processed, empty out the queue
                                    else:
                                        # update queue and neighbor dict
                                        # pop all variants in the queue
                                        for i in range(0, len(current_data_variant_neigbors_list)):
                                            variant = current_data_variant_neigbors_list[0]
                                            current_data_variant_neigbors_list.pop(0)
                                            current_data_variant_neigbors_dict.pop(variant, None)

                                        # enqueue the current line, since it is the current variant being processed
                                        # in the queue, we need to maintain the invariants.
                                        current_data_variant_neigbors_list.append(each_variant)
                                        
                                        # update the index of current variant in the queue
                                        # the current line is the first in the queue
                                        current_data_variant_index = 0
                                        
                                        # update its neigbor variants
                                        enqueue_variant(each_variant)

                                        # update sliding window
                                        current_data_variant_pos = pos
                                        window_pos_left = current_data_variant_pos - WINDOW_UP
                                        window_pos_right = current_data_variant_pos + WINDOW_DOWN
                                        
                                        # this break statement is merely for clarity: after this if branch is
                                        # hit, the while loop will exit immediately anyway.
                                        break
                            
                            # At this point, the current line has entered the queue; 
                            # the current variant has been processed;
                            # the invariants have been updated.
                            # If the current line is the last position for the current chromosome,
                            # it is time to empty out the queue: there is no more incoming variants for the current chromesome.
                            if cross_chrom:
                                # all variants before current_data_variant_index in the queue have been processed 
                                # (variant combinations have been dumped into the output file), here only the left-over
                                # are to be processed.
                                    
                                for variant in current_data_variant_neigbors_list[current_data_variant_index:]:
                                    # current_data_variant_index do not need to be maintained here, since we are looping
                                    # over each of the variant in the queue one by one
                                    
                                    # extract the mandatory information for current variant
                                    current_variant = variant
                                    variant_field = current_variant.split('_')
                                    chrom = variant_field[0]

                                    # update sliding window
                                    current_data_variant_pos = int(variant_field[1])
                                    window_pos_left = current_data_variant_pos - WINDOW_UP
                                    window_pos_right = current_data_variant_pos + WINDOW_DOWN
                                    
                                    # update the queue: left-hand side only because there is no incoming variant to extend
                                    # the right-hand side.
                                    # update the left-hand side of the variant queue: 
                                    # remove any variant moving out of the left of the new sliding window
                                    for i in range(0, len(current_data_variant_neigbors_list)):
                                        # extract the mandatory information for each variant
                                        variant = current_data_variant_neigbors_list[0]
                                        variant_field = variant.split('_')
                                        chrom = variant_field[0]
                                        temp_pos = int(variant_field[1])
                                        
                                        if temp_pos < window_pos_left:
                                            # fall to the left of the sliding window,
                                            # none of the variants in the queue would fall to the right,
                                            # because the updated window moves to the right
                                            current_data_variant_neigbors_list.pop(0)
                                            current_data_variant_index -= 1

                                            # update neigbor dict
                                            current_data_variant_neigbors_dict.pop(variant, None)
                                        else:
                                            break                        
                                    
                                    output_seq_forward_sets_dict[current_variant] = output_seq_forward_cartesian_product_sets(current_data_variant_neigbors_list, 
                                                                                    current_data_variant_neigbors_dict)
                                    
                                    print(f"output the variant combinations for {current_variant}")
                                    output_seq_to_file(current_variant)
                                    write_report_for_processed_line(current_variant, SEQUENCE_TYPE)

                                # empty the queue: everything in the queue has been processed.
                                for i in range(0, len(current_data_variant_neigbors_list)):
                                    variant = current_data_variant_neigbors_list[0]
                                    current_data_variant_neigbors_list.pop(0)
                                    current_data_variant_neigbors_dict.pop(variant, None)
                                
                                # we are back to a fresh start:
                                # update the invariants to reflect this status.
                                current_data_variant_index = -1
                                current_data_variant_pos = -1
                                window_pos_left = -1
                                window_pos_right = -1
                                cross_chrom = False

                            # each_variant not in data_variant_final_dict
                            # output the report directly
                        else:
                            write_report_for_unprocessed_line(eachline)
                        
                    except Exception:
                        assert False, "Error happens in preprocessing code"


        # if any IOerror, give an error message
        except IOError:
            FILE_OUTPUT = open(OUTPUT, 'w')
            FILE_REPORT = open(REPORT, 'w')
            FILE_OUTPUT.write('FILE OPEN ERROR\n')
            FILE_REPORT.write('FILE OPEN ERROR\n')

if __name__ == '__main__':
    main()
