#!/usr/bin/env python3.9

import argparse
import gzip
from Bio import SeqIO
from tqdm import tqdm

BASE_PAIRING = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N', '*': '*'}

###
# parse the input parameters
###
parser = argparse.ArgumentParser(description="SeqTailor DNA seqeuence extraction for genomic variants (independent) in VCF format")

parser.add_argument("-g", "--genome", help="genome sequence filename (FASTA / FASTA.GZ)")
parser.add_argument("-c", "--coordinate", choices=['0', '1'], default='1', help="coordinate indexing")
parser.add_argument("-s", "--strand", choices=['BOTH', 'FORWARD', 'REVERSE'], default='BOTH', help="strand")
parser.add_argument("-wd", "--window_down", type=int, default=25, help="window size downstream")
parser.add_argument("-wu", "--window_up", type=int, default=25, help="window size upstream")
parser.add_argument("-q", "--seq_type", choices=['BOTH', 'WT', 'MT'], default='BOTH', help="output sequence type")
parser.add_argument("-i", "--input", help="input filename (VCF)")
parser.add_argument("-o", "--output", help="output filename (FASTA)")
parser.add_argument("-r", "--report", help="report filename (TEXT)")

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
def seqlist_to_seqstr(seqeunce_list):
    sequence_string = ''
    for each_seq in seqeunce_list:
        sequence_string += each_seq
    return sequence_string


# reverse the input sequence into complementary strand in reverse order
def seqreverse(sequence_input):
    sequence_length = len(sequence_input)
    sequence_reversed = ''
    for i in range(0, sequence_length):
        sequence_reversed += BASE_PAIRING[sequence_input[sequence_length - i - 1]]
    return sequence_reversed


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
    SeqTailor_PARAMETER = '# DNA seqeuence extraction for genomic variants (independent) in VCF format\n' \
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

        FILE_INPUT = open(INPUT, 'r')
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
        FILE_INPUT.seek(0)
        for eachline in FILE_INPUT:
            if not eachline.startswith('#') and eachline.strip():
                eachline = eachline.strip()
                data_all_list.append(eachline)
                field = eachline.split()
                try:
                    # read the 4 mantatory fileds of data for each genomic variant
                    chrom = field[0]
                    pos = int(field[1])
                    ref = field[3]
                    alt = field[4]
                    each_variant = chrom + '_' + str(pos) + '_' + ref + '_' + alt

                    # check each line of input, and recored the comment, the action, and the submitted data
                    if each_variant not in data_variant_set:
                        data_variant_set.add(each_variant)

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

                        elif len(ref) > WINDOW_DOWN or len(ref) > WINDOW_UP:
                            data_all_report_comment[eachline] = 'REF longer than winodow size'
                            data_all_report_action[eachline] = 'skipped'
                            data_all_report_submit[eachline] = 'none'

                        elif ',' in alt:
                            alt_list = alt.split(',')
                            each_variant_new = ''
                            for each_alt in alt_list:
                                each_variant_temp = chrom + '_' + str(pos) + '_' + ref + '_' + each_alt
                                each_variant_new += each_variant_temp + ','
                                data_variant_final_list.append(each_variant_temp)
                            each_variant_new = each_variant_new[0:-1]
                            data_all_report_comment[eachline] = 'normal'
                            data_all_report_action[eachline] = 'normal'
                            data_all_report_submit[eachline] = each_variant_new

                        else:
                            data_all_report_comment[eachline] = 'normal'
                            data_all_report_action[eachline] = 'normal'
                            data_all_report_submit[eachline] = each_variant
                            data_variant_final_list.append(each_variant)

                    else:
                        data_all_report_comment[eachline] = 'duplicated'
                        data_all_report_action[eachline] = 'skipped'
                        data_all_report_submit[eachline] = 'none'

                except Exception:
                    data_all_report_comment[eachline] = 'field mistake'
                    data_all_report_action[eachline] = 'skipped'
                    data_all_report_submit[eachline] = 'none'

        ###
        # write the recorded comments, actions taken, submitted data, into a report
        ###
        for eachdata in data_all_list:
            comment = '\t> COMMENT: ' + data_all_report_comment[eachdata]
            action = '\t> ACTION: ' + data_all_report_action[eachdata]
            submit = '\t> SUBMITTED: ' + data_all_report_submit[eachdata]
            FILE_REPORT.write(eachdata+'\n'+comment+'\n'+action+'\n'+submit+'\n\n')

        ###
        # extract DNA sequences for genomic variants that have passed the checking step
        ###
        with tqdm(total=len(data_variant_final_list)) as pbar:
            for each_variant in data_variant_final_list:
                # extract the mandatory information for each variant
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
                # extract the reference DNA seqeunce in forward strand, and convert it to reverse strand
                # write each nucleotide into a list of reference nucleotides
                ###
                output_seq_forward_ref = chrom_sequence[forward_end_up: forward_end_down]
                output_seq_forward_ref_list = list(output_seq_forward_ref)
                output_seq_reverse_ref = seqreverse(output_seq_forward_ref.upper())
                output_seq_reverse_ref_list = list(output_seq_reverse_ref)

                # modify the list of reference nucleotides in forward strand and reverse strand
                # replicate list of ref nucleotides to list of alt nucleotides and empty the corresponding cell(s)
                output_seq_forward_alt_list = output_seq_forward_ref_list
                for temp_pos in range(pos_new_forward, pos_new_forward + len(ref)):
                    output_seq_forward_alt_list[temp_pos] = ''

                output_seq_reverse_alt_list = output_seq_reverse_ref_list
                for temp_pos in range(pos_new_reverse - len(ref) + 1, pos_new_reverse + 1):
                    output_seq_reverse_alt_list[temp_pos] = ''

                # generatt alternative nucleotides in forward strand and reverse strand
                if alt == '.' or alt == '-':
                    # already emptyed
                    pass
                else:
                    output_seq_forward_alt_list[pos_new_forward] = alt
                    alt_reverse = seqreverse(alt)
                    output_seq_reverse_alt_list[pos_new_reverse] = alt_reverse

                ###
                # concantenate the list of modified nucleotides, to generate the alternative sequences
                ###
                output_seq_forward_alt = seqlist_to_seqstr(output_seq_forward_alt_list)
                output_seq_reverse_alt = seqlist_to_seqstr(output_seq_reverse_alt_list)

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

                # write the output seqeuence to fasta file, according to the user's option on strand and sequence type
                if STRAND == 'BOTH' or STRAND == 'FORWARD':
                    if SEQUENCE_TYPE == 'BOTH' or SEQUENCE_TYPE == 'WT':
                        FILE_OUTPUT.write(output_header_forward_ref + '\n' + output_seq_forward_ref_wrap + '\n')
                    if SEQUENCE_TYPE == 'BOTH' or SEQUENCE_TYPE == 'MT':
                        FILE_OUTPUT.write(output_header_forward_alt + '\n' + output_seq_forward_alt_wrap + '\n')

                if STRAND == 'BOTH' or STRAND == 'REVERSE':
                    if SEQUENCE_TYPE == 'BOTH' or SEQUENCE_TYPE == 'WT':
                        FILE_OUTPUT.write(output_header_reverse_ref + '\n' + output_seq_reverse_ref_wrap + '\n')
                    if SEQUENCE_TYPE == 'BOTH' or SEQUENCE_TYPE == 'MT':
                        FILE_OUTPUT.write(output_header_reverse_alt + '\n' + output_seq_reverse_alt_wrap + '\n')

                FILE_OUTPUT.write('\n')

                # progress report
                pbar.update(1)

    # if any IOerror, give an error message
    except IOError:
        FILE_OUTPUT = open(OUTPUT, 'w')
        FILE_REPORT = open(REPORT, 'w')
        FILE_OUTPUT.write('FILE OPEN ERROR\n')
        FILE_REPORT.write('FILE OPEN ERROR\n')
