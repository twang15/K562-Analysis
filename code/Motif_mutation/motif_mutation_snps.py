# Mar 2, 2022
""" Code for generation of VCF file(s) for motif perturbation.
Motif disruption serves as the positive control for evaluation of signal-change detection metrics
VCF files are the output for downstream sequence construction

Methods of motif disruption:
1) complete motif deletion (VCF deletion)
2) motif scrambling (VCF in-place substitution)
3) deletion of important / less important bases of motif (VCF single-bp deletion)
4) mutation of important / less important bases of motif (VCF SNP)

The below code generates a motif disruption VCF for CTCF first.
VCF <version here> format used throughout.

File sources:
- hg38 (reference genome): https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/

- ENCODE 4 CTCF ChIP-seq data from https://www.encodeproject.org/experiments/ENCSR000EGM/

- CTCF motif (Discovered by MEME, bed format ) locations obtained from https://factorbook.org/tf/human/CTCF/motif/ENCSR000EGM/
- TF peaks (bed narrowPeak): https://www.encodeproject.org/files/ENCFF660GHM/
"""

import gzip
import subprocess
import random
import sys
import pandas as pd


def intersect_peaks_motif(bedtools_path, peak_path, motif_path):
    """ Intersect peak and motif BED files using bedtools

    Below code from k562_all_tf_stats.py

    :return: list of coordinates of motifs in peaks
    """
    # TF peaks and Motifs are both using BED format: 0-based, half-open
    # Position format (1-based fully closed) vs. BED format (0-based half open)
    ## https://genome-blog.soe.ucsc.edu/blog/2016/12/12/the-ucsc-genome-browser-coordinate-counting-systems/
    
    # keep motif, -wa
    # only when motif 100% overlap with TF peaks, -f 1
    # only report the motif once, -u
    ## see: https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html 
    overlap = subprocess.run([bedtools_path, 'intersect', '-a', motif_path, '-b', peak_path, '-wa', '-f', '1', '-u'],
                                 stdout=subprocess.PIPE)

    overlap_decoded = overlap.stdout.decode('utf-8').split('\n')    # some redundancy, i.e. multiple motifs in same peak
    motifs_in_peaks_l = [i.split('\t')[:4] for i in list(set(overlap_decoded))]  # all unique motif coordinates
    
    # keep only + strand motifs
    motifs_in_peaks_final = [[i[0], int(i[1]), int(i[2]), i[3]] for i in motifs_in_peaks_l if len(i) > 1 and i[3] == '+']

    return motifs_in_peaks_final


def load_hg38(fa_dir, chrom):
    """ Load fasta file of reference genome (hg38). Adapted from k562_custom_genome.py

    :param fa_dir: directory containing all gzipped chromosome-level fasta files (i.e. directory containing chr1.fa.gz, chr2.fa.gz, etc.)
    :param chrom: chromosome string (e.g. 'chr1')
    :return: header-excluded reference sequence of the given chromosome
    """
    fa_name = '%s.fa.gz' % chrom
    fa = subprocess.run('ls', stdout=subprocess.PIPE, cwd=fa_dir)     # fasta files list
    fa_l = fa.stdout.decode('utf-8').split('\n')
    if fa_name in fa_l:
        print('fasta file for %s found' % chrom)
        with gzip.open(fa_dir + fa_name, 'r') as fa:
            fa_str = fa.read()
    else:
        raise FileNotFoundError('Fasta file for %s cannot be found' % chrom)

    fa_str_decoded = fa_str.decode('utf-8')
    chrom_label = fa_str_decoded.split('\n')[0]
    fa_str_final = fa_str_decoded.replace(chrom_label, '').replace('\n', '')

    return fa_str_final


def make_vcf_motif_deletion(motifs_in_peaks_final, fa_dir):
    """ Generate dataframe of vcf entries from lists
    Parts from from k562_custom_genome.py

    BED format

    :param motifs_in_peaks_final: Output of intersect_peaks_motif()
    :param fa_dict: dictionary for chromosomal hg38 sequences (each key a chromosome, each value the corresponding hg38 sequence string)
    :return: Nested list, each element a vcf entry (list form)
    """
    # Double check indexing (0-based, 1-based, subtraction for motif_len)
    vcf_lines = []
    for i in range(23):
        if (i + 1) != 23:
            chr_i = 'chr%d' % (i + 1)
        else:
            chr_i = 'chrX'
        fa_chr_i = load_hg38(fa_dir, chr_i)   # Obtain chromosomal hg38 sequence
        motifs_in_peaks_chr_i = [i for i in motifs_in_peaks_final if i[0] == chr_i]     # motifs in peaks on chromosome chr_i
        chrom_i, id_all = chr_i, '.'
        for motif in motifs_in_peaks_chr_i:
            motif_sta, motif_end, strand = motif[1], motif[2], motif[3]
            ref = fa_chr_i[motif_sta - 1: motif_end].upper()
            alt = fa_chr_i[motif_sta - 1].upper()     # already '+' strand
            if ref[0] != alt:
                raise ValueError('First base of REF does not match ALT!')
            vcf_lines.append([chrom_i, motif_sta-1, id_all, ref, alt])  # first 5 VCF columns

        print('Motif deletion VCF entries for %s generated' % chr_i)

    return vcf_lines
  
def get_mutation(index, reference, mutation_map):
    '''
    :param index: the 0-based index of the current base
    :param reference: the reference base
    :param mutation_map: pandas data frame, the map to transform each motif base.
    :return: all the possible bases for the current reference base to transform into. 
    :Note: it could be empty if the index of a base does not exist in mutation map.
    '''
    try:
        # the original index was 1-based R index.
        result = mutation_map.loc[index+1]
    except KeyError:
        result = None

    return result

def mutate_motif_bases(motifs_in_peaks_final, fa_dir, mutation_map):
    """
    For every base in motif (crude approach), mutate to all 3 other bases (comma delimited ALT entries)

    :param motifs_in_peaks_final:
    :param pos: 1-based desired motif (not hg38) position (will be ocnverted to 0-based) for mutation
    :param fa_dir:
    :param mutation_map: the path to the mutation map tsv file
    :return:
    """

    # read the mutation map
    df = pd.read_csv(mutation_map, sep='\t', index_col='index')

    vcf_mutation_lines = []
    for i in range(23):
        if (i + 1) != 23:
            chr_i = 'chr%d' % (i + 1)
        else:
            chr_i = 'chrX'
        fa_chr_i = load_hg38(fa_dir, chr_i)   # Obtain chromosomal hg38 sequence
        motifs_in_peaks_chr_i = [i for i in motifs_in_peaks_final if i[0] == chr_i]     # motifs in peaks on chromosome chr_i
        chrom_i, id_all = chr_i, '.'
        for motif in motifs_in_peaks_chr_i:
            motif_sta, motif_end, strand = motif[1], motif[2], motif[3]
            motif_len = motif_end - motif_sta       # absolute length (double check this)
            for motif_pos in range(motif_len):
                ref = fa_chr_i[motif_sta + motif_pos].upper()
                # possible_bases = ['A', 'C', 'G', 'T']
                # possible_bases.remove(ref)
                # vcf_mutation_lines.append([chrom_i, motif_sta+motif_pos, id_all, ref, ','.join(possible_bases)])  # first 5 VCF columns
                
                # 5 + [probability_source, probability_target]
                mutation = get_mutation(index=motif_pos, reference=ref, mutation_map=df)
                if mutation is not None:
                    vcf_mutation_lines.append([chrom_i, motif_sta+motif_pos, id_all, ref, 
                                            mutation['target'], 
                                            motif_pos,
                                            mutation['probability_source'],
                                            mutation['probability_target']])

            # all alternative bases in ALT entry (comma delimited)
            # motif_sta+1 converts motif_sta into 1-based motif start coordinate
            # motif_sta+1+motif_pos is the 1-based hg38 coordinate of the base to mutate
            # motif_pos is "0-based" because the first base would dbe motif_sta + 0, and so on.

        print('Motif mutations VCF entries for %s generated' % chr_i)

    return vcf_mutation_lines


def to_vcf(vcf_lines_l, vcf_output_dir, factor_name, motif_perturb_type):
    """ Generate vcf file

    :param vcf_output_path: full path for VCF file output directory (not file), including final '/'
    :param factor_name: name of TF (e.g. 'CTCF')
    :return:
    """
    with open(vcf_output_dir + '%s_%s.vcf' % (factor_name, motif_perturb_type), 'w') as vcf:
        nrows = len(vcf_lines_l)
        for row_i in range(nrows):
            row = vcf_lines_l[row_i]
            num_cols = len(row)
            col_i = 0
            for entry in row:
                if col_i < num_cols - 1:
                    col_i += 1
                    vcf.write(str(entry) + '\t')
                elif col_i == num_cols - 1:
                    vcf.write(str(entry) + '\n')


# Example for running code

""" File paths """
#bedtools_path = '/scg/apps/legacy/software/bedtools/2.27.1/bedtools2/bin/bedtools'
bedtools_path = '/opt/local/bin/bedtools'

prefix='/oak/stanford/scg/prj_ENCODE/ChIP/basepairmodel/basepairmodels/tutorial/Data_from_Shannon/'
prefix='/Users/taowang9/Documents/2022Spring/ChIP/Motif_mutation/'
peak_path =  prefix + sys.argv[3] + '/' + sys.argv[1]
motif_path = prefix + sys.argv[3] + '/' + sys.argv[2]

fa_directory = '/oak/stanford/scg/prj_ENCODE/Long-read-RNA/hg38/gziped/'
fa_directory = prefix +'hg38/gziped/'
vcf_outdir = prefix + 'Output.mutation/' + sys.argv[3] + '/'

mutation_map = prefix + sys.argv[3] + '/' + sys.argv[4]

overlaps_final = intersect_peaks_motif(bedtools_path, peak_path, motif_path)
#vcf_lines_l = make_vcf_motif_deletion(overlaps_final, fa_directory)
vcf_lines_l = mutate_motif_bases(overlaps_final, fa_directory, mutation_map)
#to_vcf(vcf_lines_l, vcf_outdir, sys.argv[3], 'deletion')
to_vcf(vcf_lines_l, vcf_outdir, sys.argv[3], 'mutation')
