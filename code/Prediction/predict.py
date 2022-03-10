'''
Model wrapper for basepairmodel

It takes a sequence file (in the format defined by
our universal genome sequence tailor tool) and makes prediction for TF-binding signal profile
for each input sequence.

Written by Tao Wang, Department of Genetics, School of Medicine, Stanford University, CA, USA.

Contributors: 
    Jay Luo, Johns Hopkins University, Baltimore, MD, USA 
    Shannon White, Department of Genetics, School of Medicine, Stanford University, CA, USA.
'''

from concise.preprocessing import encodeDNA
import requests
import argparse
import os
import sys
import logging

from basepairmodels.cli import bigwigutils
from basepairmodels.cli import logger
from basepairmodels.cli import MTBatchGenerator
from basepairmodels.cli.batchgenutils import *
from basepairmodels.cli.bpnetutils import *
from basepairmodels.cli.losses import MultichannelMultinomialNLL, multinomial_nll

from keras.models import load_model
from keras.utils import CustomObjectScope
from mseqgen import quietexception
from tqdm import tqdm
import pandas as pd
from scipy.special import logsumexp
import numpy as np

# argparse: to accept command line options
# input: the input sequence will be read from a file.
def argsparser_input():
    """ Command line arguments for the predict script
    for our wrapper
    """

    parser = argparse.ArgumentParser()
    
    # path to the pre-trained basepair model
    parser.add_argument('--model', '-m', type=str,
                        help="path to the .h5 model file.")
    
    # predict modes
    parser.add_argument('--predict-peaks', action='store_true', 
                        help="generate predictions only on the peaks "
                        "contained in the peaks.bed files")
    
    # reference params
    parser.add_argument('--reference-genome', '-g', type=str, required=True,
                        help="the path to the reference genome fasta file")
    
    parser.add_argument('--chrom-sizes', '-s', type=str, required=True,
                        help="path to chromosome sizes file")

    # input sequence file name
    parser.add_argument('--target-sequence-file', '-f', type=str, required=True,
                        help='the input file containing target sequence. Default to example.tsv')

    # input sequence length
    parser.add_argument('--sequence-length', '-l', type=int, required=False, default=2114,
                        help='the length of the single-line input. Default to 2114.')

    
    # output params
    parser.add_argument('--output-len', type=int, 
                        help="length of output profile", default=1000)
    
    parser.add_argument('--output-dir', '-o', type=str, required=True,
                        help="destination directory to store predictions as a "
                        "bigWig file")
    
    parser.add_argument('--automate-filenames', action='store_true', 
                        help="specify if the predictions output should "
                        "be stored in a timestamped subdirectory within "
                        "--output-dir")
    
    parser.add_argument('--time-zone', type=str,
                        help="time zone to use for timestamping model "
                        "directories", default='US/Pacific')
    
    parser.add_argument('--exponentiate-counts', action='store_true', 
                        help="specify if the predicted counts should be "
                        "exponentiated before writing to the bigWig files")

    parser.add_argument('--output-window-size', type=int,
                        help="size of the central window of the output "
                        "profile predictions that will be written to the "
                        "bigWig files", default=1000)
    
    # misc params
    parser.add_argument('--write-buffer-size', type=int,
                        help="size of the write buffer to store predictions "
                        "before writing to bigWig files", default=10000)

    return parser

def predict(args, input_data, pred_dir):
    # load the model
    model = load_model(args.model)

    # tsv
    TSV_FILE_OUTPUT = open(pred_dir+'/prediction.tsv', 'w')
    TSV_FILE_HEADER = "Sequence_num\tRef_total_count_N\tVariant_total_count_N\tRef_signal_values\tVariant_signal_values\n"
    TSV_FILE_OUTPUT.write(TSV_FILE_HEADER)

    output_flank = int(args.output_len/2)
    batch_size = 1
    control_smoothing = 2
    control_profile = np.zeros((batch_size, output_flank*2, control_smoothing), dtype=np.float32)
    control_profile_counts = np.zeros((1), dtype=np.float32)

    # predict
    input_df = pd.read_csv(input_data, sep='\t')
    for row in range(0, input_df.shape[0]):
        # validate the input sequence 
        seq = input_df.loc[row]['Window_sequence']
        seq_ref = input_df.loc[row]['Window_reference']
        row_num = input_df.loc[row]['Row_num']
        
        if len(seq) != 2114 or len(seq_ref) != 2114:
            raise quietexception.QuietException("Sequence length is not 2114 for line {}".format(row_num))

        # one-hot encoding for sequence:
        # references: 
        # 1. https://www.cmm.in.tum.de/public/docs/concise/preprocessing/sequence/
        # 2. https://colab.research.google.com/drive/1VNsNBfugPJfJ02LBgvPwj-gPK0L_djsD#scrollTo=ovPKEvmyu-Zl
        # seq = 'N'*23 + 'GGAGGAACTGGGTGTGGGGAGGTTGTAGCCCGACCCTGCCCCTCCCCCCAGGGAGGTTGAGAGTTCTGGGCAGACGGCAGATGCATAACAAAGGTGCATGATAGCTCTGCCCTGGGGGCAGAGAAGATGGTTGGGGAGGGGTCCCTCTCGTCCTA' + 'N'*22
        seq_onehot = encodeDNA([seq.upper()]) # one-hot encode
        seq_ref_onehot = encodeDNA([seq_ref.upper()])
       
        # format the input
        ## reference: https://github.com/kundajelab/basepairmodels/blob/8aaff358a4e36041e322829cacc2cddf9bf8d3cd/basepairmodels/cli/MTBatchGenerator.py#L379
        seq_input_dict = {'sequence': seq_onehot, 
                           'control_profile': control_profile,
                           'control_logcount': control_profile_counts}
        
        seq_ref_input_dict = {'sequence': seq_ref_onehot, 
                           'control_profile': control_profile,
                           'control_logcount': control_profile_counts}
        
        # predict
        predictions = model.predict(seq_input_dict)
        predictions_ref = model.predict(seq_ref_input_dict)
        
        # extract profile and counts
        ## Note: 
        ##  1. there is only 1000-bp long prediction for the center of 2114 input window
        ##  2. the profile contains logit values (raw output from the model) after softmax operation gives
        ##      the probabilities of each base pair (see https://github.com/twang15/K562-Analysis/issues/17#issuecomment-967701919)
        ##  3. the counts are the same number log(N) for positions in the 1000 prediction window and needs to be exponentiated
        ##      before final p*N calculation for signal values
        total_count_N = np.exp(predictions[1][0][0])
        logits_profile = predictions[0][0][:, 0]
        ref_total_count_N = np.exp(predictions_ref[1][0][0])
        ref_logits_profile = predictions_ref[0][0][:, 0]

        # softmax operation to convert logits to log(probabilities), then exponentiate to probabilities
        probability_profile = np.exp(logits_profile - logsumexp(logits_profile))
        ref_probability_profile = np.exp(ref_logits_profile - logsumexp(ref_logits_profile))

        # get per-base count (binding strength) profile
        signal_profile = np.multiply(total_count_N, probability_profile) 
        ref_signal_profile = np.multiply(ref_total_count_N, ref_probability_profile) 

        # output I/O
        # Analysis:
        ## 1. bigwig file format is easy to visualize 
        ## 2. pandas dataframe: Sequence_num, Ref_total_count_N, Ref_signal_value, Variant_total_count_N, Variant_signal_value
        ##      Sequence_num: index (0-based) of the sequence in the input dataframe input_df
        ##      Ref_total_count_N: total count N for the reference sequence
        ##      Ref_signal_values: p*Ref_total_count_N
        ##      Variant_total_count_N: total count N for the sequence with variants
        ##      Variant_signal_values: p*Variant_total_count_N
        Sequence_num_tsv = row
        Ref_total_count_N_tsv = ref_total_count_N
        Variant_total_count_N_tsv = total_count_N
        Ref_signal_values_tsv = ",".join(["%.6f" % e for e in ref_signal_profile.tolist()])
        Variant_signal_values_tsv = ",".join(["%.6f" % e for e in signal_profile.tolist()])
        
        statistics = False
        if statistics:
            # As for the metric for filtering most interesting sequence variant, area under curve could be one candidate.
            ## There are also other candidates, such as
            ## 1. Log2 fold-change Predicted Read Counts between Variant and reference sequence, and P-value for the two groups.
            ## 2. similarity metrics: https://pypi.org/project/similaritymeasures/
            ## 3. Goodness of fit
            ## 4. Correlation
            pass
        else:
            TSV_FILE_OUTPUT.write(f'{Sequence_num_tsv}\t{Ref_total_count_N_tsv}\t{Variant_total_count_N_tsv}'
                              f'\t{Ref_signal_values_tsv}\t{Variant_signal_values_tsv}\n')

        #with open(args.target_sequence_file, 'r') as target_file:
            #lines = target_file.readlines()
        ## extract the context sub-sequence from reference genome for target motif
        #response = requests.get(url="http://togows.org/api/ucsc/hg38/chr1:100000-100010")
        #print(response.text)

def predict_main():
    """ The main entry to our wrapper predictor.
    """

    # parse the command-line arguments
    parser = argsparser_input()
    args = parser.parse_args()

    # check if the output directory exists
    if not os.path.exists(args.output_dir):
        logging.error("Directory {} does not exist".format(args.output_dir))
    
    if args.automate_filenames:
        # create a new directory using current date/time to store the
        # predictions and logs
        date_time_str = local_datetime_str(args.time_zone)
        pred_dir = '{}/{}'.format(args.output_dir, date_time_str)
        os.mkdir(pred_dir)
    elif os.path.isdir(args.output_dir):
        pred_dir = args.output_dir        
    else:
        logging.error("Directory does not exist {}.".format(args.output_dir))
        return
    
    # set up the loggers
    logfname = "{}/predict.log".format(pred_dir)
    logger.init_logger(logfname)

    # make sure that the target sequence file exists
    if not os.path.isfile(args.target_sequence_file):
        raise quietexception.QuietException(
                "File not found: {} OR you may have accidentally "
                "specified a directory path.".format(args.target_sequence_file))

    input_data = args.target_sequence_file
    logging.info("Target Sequence DATA -\n{}".format(input_data))
    
    # check sequence length
    if args.sequence_length != 2114:
        raise quietexception.QuietException("Sequence length should be 2114.")

    # predict
    logging.info("Loading {}".format(args.model))
    with CustomObjectScope({'MultichannelMultinomialNLL': 
                            MultichannelMultinomialNLL}):
            
        predict(args, input_data, pred_dir)

if __name__ == '__main__':
    predict_main()
