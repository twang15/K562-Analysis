""" Generate useful statistics for all hg19 TF datasets for K562
"""
import time
import numpy as np
import vcf_encode
import subprocess
import gzip
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt



def parse_metadata(path_metadata_file):
    """ Inspect metadata, identify multiple experiments on same factor, keep newest one
    """
    factors_id_dict = {}
    factors = []
    id_keep = []
    with open(path_metadata_file, 'r') as metadata:
        metadata_lines = metadata.readlines()
        metadata_lines_new = [i.rstrip().split('\t') for i in metadata_lines]
        metadata_df = pd.DataFrame(metadata_lines_new[1:])
        metadata_df.columns = metadata_lines_new[0]
        metadata_df = metadata_df.sort_values(by=['Experiment date released'], ascending=False).reset_index(drop=True)
        num_lines_metadata = metadata_df.shape[0]
        for row_ind in range(num_lines_metadata):
            factor_i = metadata_df.iloc[row_ind, 22].split('-')[0]
            id_i = metadata_df.iloc[row_ind, 0]
            if factor_i not in factors:
                factors.append(factor_i)
                id_keep.append(id_i)
                factors_id_dict[id_i] = factor_i
            else:
                continue

    return id_keep, factors_id_dict


def filenames(narrowpeak_dir, keep_id):
    """ Get filenames of narrowPeak files
    """
    keep_l = [i + '.bed.gz' for i in keep_id]
    list_files = subprocess.run('ls', stdout=subprocess.PIPE, cwd=narrowpeak_dir)
    list_files_l = list_files.stdout.decode('utf-8').split('\n')        # files list
    list_filenames_peaks = [i for i in list_files_l if '.bed' in i]
    list_filenames_peaks_filt = [i for i in list_filenames_peaks if i in keep_l]

    return list_filenames_peaks_filt


def regions(filename, narrowpeak_dir):
    """ Return dataframe of regions for a gzipped narrowPeak file
    'filename' includes .bed.gz extension

    Note: since variant coordinates are converted to 0-based in vcf_encode.py and the coordinates
    in narrowPeak files are 0-based (excluding end base), no further changes to the coordinate themselves
    are needed.
    """
    with gzip.open(narrowpeak_dir + filename, 'r') as bed:
        bed_lines = bed.readlines()
        bed_lines = [i.decode('utf-8') for i in bed_lines]
        bed_lines_final = [i for i in bed_lines if 'chr' in i]        # in case any blank lines
        bed_reg = [i.split('\t')[:3] for i in bed_lines_final]

        return pd.DataFrame(bed_reg)    # all chromosomes for dataset 'filename'


def filt_chr(reg_df, chrom):
    """ Filter regions dataframe (containing all chromsomes) for chromosome 'chrom'
    """
    return reg_df[reg_df.iloc[:, 0] == chrom]


def num_peaks_per_dataset(filename, narrowpeak_dir):
    """ Count number of peaks for one ENCODE K562 TF dataset (hg19)
    """
    return regions(filename, narrowpeak_dir).shape[0]


def num_peaks_per_chr_per_dataset(filename, narrowpeak_dir):
    num_peaks_per_chr = []
    df = regions(filename, narrowpeak_dir)
    for i in range(1, 23):
        if i != 23:
            chr_i = 'chr%d' % i
            num_peaks_per_chr.append(df[df.iloc[:, 0] == chr_i].shape[0])
        else:
            chr_i = 'chrX'
            num_peaks_per_chr.append(df[df.iloc[:, 0] == chr_i].shape[0])

    return num_peaks_per_chr


def get_var_all(vcf_paths, col_names_l):
    return vcf_encode.load_vcf(vcf_paths, 'all', col_names_l)


def get_snp_df(vcf_paths, col_names_l):
    var_all = get_var_all(vcf_paths, col_names_l)
    snps_only = vcf_encode.vcf_filt(var_all, 'SNP', 'all')

    return snps_only


def get_coor(var_df_all, gt, chrom):
    """ Generate variant coordinates from VCF files
    """
    # var_df_all = vcf_encode.load_vcf(vcf_paths, chrom, col_names_l)
    snps_only = vcf_encode.vcf_filt(var_df_all, 'SNP', chrom)
    ins_only = vcf_encode.vcf_filt(var_df_all, 'INS', chrom)
    del_only = vcf_encode.vcf_filt(var_df_all, 'DEL', chrom)
    dup_only = vcf_encode.vcf_filt(var_df_all, 'DUP', chrom)
    inv_only = vcf_encode.vcf_filt(var_df_all, 'INV', chrom)  # TODO: check why it's not used (error?) [done]
    bnd_only = vcf_encode.vcf_filt(var_df_all, 'BND', chrom)
    bnd_intrachromosomal = vcf_encode.filt_bnd(bnd_only, chrom)[0]
    bnd_interchromosomal = vcf_encode.filt_bnd(bnd_only, chrom)[1]

    # Add '<INV>' alleles to bnd_intrachromosomal
    bnd_intrachromosomal = pd.concat([bnd_intrachromosomal, inv_only]).reset_index(drop=True)

    # Get maternal and paternal dataframes
    bnd_inter_maternal, bnd_inter_paternal = vcf_encode.bnd_subset_gt(bnd_interchromosomal)
    bnd_intra_maternal, bnd_intra_paternal = vcf_encode.bnd_subset_gt(bnd_intrachromosomal)

    ins_mult_alt_m, ins_mult_alt_p = vcf_encode.mult_alt_df_gt(ins_only)  # still containing some SNPs
    snps_mult_alt_m = vcf_encode.vcf_filt(ins_mult_alt_m, 'SNP', chrom)
    snps_mult_alt_p = vcf_encode.vcf_filt(ins_mult_alt_p, 'SNP', chrom)

    ins_maternal, ins_paternal = vcf_encode.filt_gt(ins_only)
    ins_only_mult_alt_m, ins_only_mult_alt_p = vcf_encode.vcf_filt(ins_mult_alt_m, 'INS', chrom), vcf_encode.vcf_filt(ins_mult_alt_p, 'INS',
                                                                                                chrom)
    ins_merged_m, ins_merged_p = vcf_encode.merge_df(ins_maternal, ins_only_mult_alt_m), vcf_encode.merge_df(ins_paternal,
                                                                                       ins_only_mult_alt_p)

    # Deletions
    del_maternal, del_paternal = vcf_encode.filt_gt(del_only)

    # SNPs
    snps_maternal, snps_paternal = vcf_encode.filt_gt(snps_only)
    snps_merged_m, snps_merged_p = vcf_encode.merge_df(snps_maternal, snps_mult_alt_m), vcf_encode.merge_df(snps_paternal, snps_mult_alt_p)

    # Duplications
    dups_maternal, dups_paternal = vcf_encode.filt_gt(dup_only)

    # Get coordinates for maternal and paternal dataframes
    coord_inter_m, coord_inter_p = vcf_encode.ind_bnd_inter(bnd_inter_maternal), vcf_encode.ind_bnd_inter(bnd_inter_paternal)
    coord_intra_m, coord_intra_p = vcf_encode.ind_bnd_intra(bnd_intra_maternal), vcf_encode.ind_bnd_intra(bnd_intra_paternal)
    coord_ins_m, coord_ins_p = vcf_encode.ind_ins(ins_merged_m), vcf_encode.ind_ins(ins_merged_p)
    coord_del_m, coord_del_p = vcf_encode.ind_del(del_maternal), vcf_encode.ind_del(del_paternal)
    coord_snps_m, coord_snps_p = vcf_encode.ind_snps(snps_merged_m), vcf_encode.ind_snps(snps_merged_p)
    coord_dups_m, coord_dups_p = vcf_encode.ind_dups(dups_maternal), vcf_encode.ind_dups(dups_paternal)

    # # Get filtered coordinates
    # coord_intra_filt_m, rm_ind_intra_m = vcf_encode.coor_rm_intra(coord_inter_m, coord_intra_m)
    # coord_intra_filt_p, rm_ind_intra_p = vcf_encode.coor_rm_intra(coord_inter_p, coord_intra_p)
    #
    # coord_ins_filt_m, rm_ind_ins_m = vcf_encode.coor_rm_ins(coord_intra_filt_m, coord_ins_m)
    # coord_ins_filt_p, rm_ind_ins_p = vcf_encode.coor_rm_ins(coord_intra_filt_p, coord_ins_p)
    #
    # coord_del_filt_m, rm_ind_del_m = vcf_encode.coor_rm_del(coord_inter_m, coord_ins_filt_m, coord_del_m)
    # coord_del_filt_p, rm_ind_del_p = vcf_encode.coor_rm_del(coord_inter_p, coord_ins_filt_p, coord_del_p)

    # Rationale for using 'raw' coordinates: each peak-variant overlap is one separate instance

    if gt == 'm':
        return coord_intra_m, coord_snps_m, coord_ins_m, \
           coord_del_m, coord_dups_m
    elif gt == 'p':
        return coord_intra_p, coord_snps_p, coord_ins_p, \
            coord_del_p, coord_dups_p


# individual functions to compute number of each variant type in peaks

def num_intra_in_peaks(peaks_df, coor_intra):
    """ Overlap reference (*): https://stackoverflow.com/questions/3269434/whats-the-most-efficient-way-to-test-two-integer-ranges-for-overlap
    """
    overlap = 0
    num_peaks = peaks_df.shape[0]
    for i in range(num_peaks):
        row = peaks_df.iloc[i, :]
        peak_chr, peak_sta, peak_end = row[0], row[1], row[2] - 1
        for j in coor_intra:
            var_sta, var_end, var_type = j[0], j[1], j[2]
            if (var_type == 'del') or (var_type == 'inv'):
                if int(var_sta) <= int(peak_end) and int(peak_sta) <= int(var_end):     # * (reference in docstring)
                    overlap += 1
            else:
                continue
        # print('%d of %d peaks processed' % (i+1, num_peaks))

    return overlap


def num_snps_in_peaks(peaks_df, coor_snp):
    overlap = 0
    num_peaks = peaks_df.shape[0]
    for i in range(num_peaks):
        row = peaks_df.iloc[i, :]
        peak_chr, peak_sta, peak_end = row[0], int(row[1]), int(row[2]) - 1
        for j in coor_snp:
            snp_coor = int(j[0])
            if peak_sta <= snp_coor <= peak_end:
                overlap += 1
            else:
                continue

    return overlap


def num_ins_in_peaks(peaks_df, coor_ins):
    overlap = 0
    num_peaks = peaks_df.shape[0]
    for i in range(num_peaks):
        row = peaks_df.iloc[i, :]
        peak_chr, peak_sta, peak_end = row[0], int(row[1]), int(row[2]) - 1
        for j in coor_ins:
            ins_coor = int(j[0])      # 'ins_coor' is base preceding actual insertion
            if peak_sta <= ins_coor <= peak_end:
                overlap += 1
            else:
                continue

    return overlap


def num_del_in_peaks(peaks_df, coor_del):
    """ non-BND represented deletions

    Overlap reference (*): https://stackoverflow.com/questions/3269434/whats-the-most-efficient-way-to-test-two-integer-ranges-for-overlap
    """
    overlap = 0
    num_peaks = peaks_df.shape[0]
    for i in range(num_peaks):
        row = peaks_df.iloc[i, :]
        peak_chr, peak_sta, peak_end = row[0], row[1], row[2] - 1
        for j in coor_del:
            var_sta, var_end = j[0], j[1]
            if int(var_sta) <= int(peak_end) and int(peak_sta) <= int(var_end):  # * (reference in docstring)
                    overlap += 1
            else:
                continue

    return overlap


def num_dups_in_peaks(peaks_df, coor_dups):
    """ For now, for a duplication to be "in a peak", the to-be-duplicated region must have some overlap with a peak

    :param peaks_df:
    :param coor_dups:
    :return:
    """
    overlap = 0
    num_peaks = peaks_df.shape[0]
    for i in range(num_peaks):
        row = peaks_df.iloc[i, :]
        peak_chr, peak_sta, peak_end = row[0], row[1], row[2] - 1
        for j in coor_dups:
            var_sta, var_end = j[0], j[1]
            if int(var_sta) <= int(peak_end) and int(peak_sta) <= int(var_end):  # * (reference in docstring)
                    overlap += 1
            else:
                continue

    return overlap


def deduplicate_coor(coor):
    new = []
    for i in coor:
        if i not in new:
            new.append(i)

    return new


def coor_all_chr(var_df_all, gt):
    intra_coor_all, snps_coor_all, ins_coor_all, del_coor_all, dups_coor_all, inv_coor_all = [], [], [], [], [], []
    for i in range(1, 24):
        if i != 23:
            chrom_i = 'chr%d' % i
        else:
            chrom_i = 'chrX'

        # chromosome-specific variant coordinates
        coor_intra, coor_snps, coor_ins, coor_del, coor_dups = get_coor(var_df_all, gt, chrom_i)

        coor_intra_del = [i[:2] for i in coor_intra if i[2] == 'del']   # i[:3] to exclude 'del' label for BED file coordinates
        coor_del = coor_del + coor_intra_del
        coor_inv = [i[:2] for i in coor_intra if i[2] == 'inv']     # i[:3] to exclude 'inv' label for BED file coordinates
        coor_snps = [[i[0], i[0] + 1] for i in coor_snps]     # for 3-column BED file (i[1] + 1 to follow BED coordinate conventions)
        coor_ins = [[i[0], i[0] + 1] for i in coor_ins]

        # for reg in coor_intra:        # already partitioned coor_intra into deletions or inversions
        #     reg.insert(0, chrom_i)
        for reg in coor_snps:
            reg.insert(0, chrom_i)
        for reg in coor_ins:
            reg.insert(0, chrom_i)
        for reg in coor_del:
            reg.insert(0, chrom_i)
        for reg in coor_dups:
            reg.insert(0, chrom_i)
        for reg in coor_inv:
            reg.insert(0, chrom_i)

        # intra_coor_all.append(coor_intra)
        snps_coor_all.append(coor_snps)
        ins_coor_all.append(coor_ins)
        del_coor_all.append(coor_del)
        dups_coor_all.append(coor_dups)
        inv_coor_all.append(coor_inv)

    return snps_coor_all, ins_coor_all, del_coor_all, dups_coor_all, inv_coor_all


def gen_bed_df_var(snps_coor_all, ins_coor_all, del_coor_all, dups_coor_all, inv_coor_all):
    """ Generate BED dataframes of each variant type

    Input is output of coor_all_chr()
    """
    # df_list_chr_intra = []        # list of dataframes, each item a variant dataframe for one chromosome
    df_list_chr_snps = []
    df_list_chr_ins = []
    df_list_chr_del = []
    df_list_chr_dups = []
    df_list_chr_inv = []

    # for chr_var in intra_coor_all:
    #     df_list_chr_intra.append(pd.DataFrame(pd.DataFrame(chr_var)))
    for chr_var in snps_coor_all:
        df_list_chr_snps.append(pd.DataFrame(pd.DataFrame(chr_var)))    # TODO: pd.DataFrame() twice?
    for chr_var in ins_coor_all:
        df_list_chr_ins.append(pd.DataFrame(pd.DataFrame(chr_var)))
    for chr_var in del_coor_all:
        df_list_chr_del.append(pd.DataFrame(pd.DataFrame(chr_var)))
    for chr_var in dups_coor_all:
        df_list_chr_dups.append(pd.DataFrame(pd.DataFrame(chr_var)))
    for chr_var in inv_coor_all:
        df_list_chr_inv.append(pd.DataFrame(pd.DataFrame(chr_var)))

    # intra_df_final = pd.concat(df_list_chr_intra).reset_index(drop=True)
    snps_df_final = pd.concat(df_list_chr_snps).drop_duplicates().reset_index(drop=True)
    ins_df_final = pd.concat(df_list_chr_ins).drop_duplicates().reset_index(drop=True)
    del_df_final = pd.concat(df_list_chr_del).drop_duplicates().reset_index(drop=True)
    dups_df_final = pd.concat(df_list_chr_dups).drop_duplicates().reset_index(drop=True)
    inv_df_final = pd.concat(df_list_chr_inv).drop_duplicates().reset_index(drop=True)

    return snps_df_final, ins_df_final, del_df_final, dups_df_final, inv_df_final


def gen_bed_files(df_l, output_dir, gt):
    """ Convert 3-column BED dataframe to BED file (3 columns as required by BPNet)

    :param df_l: [snps_df_final, ins_df_final, del_df_final, dups_df_final, inv_df_final] (list of output from gen_bed_df_var())
    :param output_dir:
    :return:
    """
    bed_filenames = ['snps_%s' % gt, 'ins_%s' % gt, 'del_%s' % gt, 'dups_%s' % gt, 'inv_%s' % gt]
    for bed_df, bed_name in zip(df_l, bed_filenames):
        vcf_encode.df_to_bed(bed_df, output_dir, bed_name)


def percent_var_in_peaks(total_num, peaks_df, cols, paths, gt):
    """ Compute percent of all variants in peaks (inclusive start/end coordinates) for all datasets. Also
    return percent of each variant in peaks.
    """
    num_intra, num_snps, num_ins, num_del, num_dups = 0, 0, 0, 0, 0
    var_all = get_var_all(paths, cols)

    # TODO: use completely unfiltered coordinates?
    intra_coor_all, snps_coor_all, ins_coor_all, del_coor_all, dups_coor_all = coor_all_chr(var_all, gt)
    for i in range(1, 24):
        if i != 23:
            chrom_i = 'chr%d' % i
        else:
            chrom_i = 'chrX'

        # chromosome-specific variant coordinates
        # coor_intra, coor_snps, coor_ins, coor_del, coor_dups = get_coor(paths, cols, gt, chrom_i)
        # print('Coordinates generated for %s of %s!' % (chrom_i, gt))
        coor_intra, coor_snps, coor_ins, coor_del, coor_dups = \
            intra_coor_all[i], snps_coor_all[i], ins_coor_all[i], del_coor_all[i], dups_coor_all[i]
        # compute number of variants in peaks per chromosome
        peaks_df = peaks_df[peaks_df.iloc[:, 0] == chrom_i].reset_index(drop=True)
        t_1 = time.time()
        num_intra_chr = num_intra_in_peaks(peaks_df, deduplicate_coor(coor_intra))
        print('num_intra_in_peaks(): %.3f seconds' % (time.time() - t_1))

        t_2 = time.time()
        num_snps_chr = num_snps_in_peaks(peaks_df, deduplicate_coor(coor_snps))
        print('num_snps_in_peaks(): %.3f seconds' % (time.time() - t_2))

        t_3 = time.time()
        num_ins_chr = num_ins_in_peaks(peaks_df, deduplicate_coor(coor_ins))
        print('num_ins_in_peaks(): %.3f seconds' % (time.time() - t_3))

        t_4 = time.time()
        num_del_chr = num_del_in_peaks(peaks_df, deduplicate_coor(coor_del))
        print('num_del_in_peaks(): %.3f seconds' % (time.time() - t_4))

        t_5 = time.time()
        num_dups_chr = num_dups_in_peaks(peaks_df, deduplicate_coor(coor_dups))
        print('num_dups_in_peaks(): %.3f seconds' % (time.time() - t_5))

        num_intra += num_intra_chr
        num_snps += num_snps_chr
        num_ins += num_ins_chr
        num_del += num_del_chr
        num_dups += num_dups_chr

    num_total_in_peaks = num_intra + num_snps + num_ins + num_del + num_dups

    # variant coordinates for all chromosomes (for counting total number of variants per variant type)
    coord_intra_all, coord_snps_all, coord_ins_all, coord_del_all, coord_dups_all = get_coor(peaks_df, gt, chrom_i)

    percent_intra = (num_intra / len(coord_intra_all))*100
    percent_snps = (num_snps / len(coord_snps_all))*100
    percent_ins = (num_ins / len(coord_ins_all))*100
    percent_del = (num_del / len(coord_del_all))*100
    percent_dups = (num_dups / len(coord_dups_all))*100

    return (num_total_in_peaks / total_num) * 100, [percent_intra, percent_snps, percent_ins,
                                                    percent_del, percent_dups]


def to_excel(df, excel_file_path):
    """ Store statistics in excel file
    """
    df.to_excel(excel_file_path)


def count_num_peaks_all(narrowpeak_filenames, dir_narrowpeak, excel_dir_path, id_factor_dict):
    all_lines = []
    # num_total_peaks = []
    for filename in narrowpeak_filenames:
        num_peaks = num_peaks_per_dataset(filename, dir_narrowpeak)
        # num_total_peaks.append(num_peaks)
        all_lines.append([filename, id_factor_dict[filename.split('.')[0]], num_peaks])

    num_total_peaks_df = pd.DataFrame(all_lines)
    to_excel(num_total_peaks_df, excel_dir_path + 'peaks_per_sample.xlsx')
    print('Number of peaks per sample excel file generated!')

    return all_lines


def count_peaks(path_metadata, dir_narrowpeak, excel_dir_path):
    keep_id, id_factor_dict = parse_metadata(path_metadata)
    narrowpeak_filenames = filenames(dir_narrowpeak, keep_id)
    num_total_peaks = count_num_peaks_all(narrowpeak_filenames, dir_narrowpeak, excel_dir_path, id_factor_dict)   # all chromosomes

    return num_total_peaks      # list, each item for one ENCODE dataset


def plot_sample_num_peaks_var_in_peaks(num_total_peaks, percentages_l, fig_dir):
    """ Get scatterplot relating number of peaks and % of variant containing peaks per sample

    :param num_total_peaks: from output of count_peaks()
    :param percentages_l: from output of compute_percent_peaks_per_dataset()
    :return:
    """
    peak_num_l = []      # nested list, [<# of peaks for file 1>, <percent of peaks containing variants for file 1> and so on]
    percent_l = []
    # do sanity checks of peak file correspondence
    for i, j in zip(num_total_peaks, percentages_l):
        filename_i, filename_j = i[0], j[0]
        if filename_i == filename_j:
            num_peaks_i, percentage_j = i[2], j[2]
            peak_num_l.append(num_peaks_i)
            percent_l.append(percentage_j)

    plt.scatter(peak_num_l, percent_l)
    plt.suptitle('Relation between number of peaks and percent of peaks containing variants')
    plt.xlabel('Number of peaks')
    plt.ylabel('Percent variant-containing')
    plt.savefig(fig_dir + 'num_peaks_percent_contain_var.png')


def run_bedtools(narrowpeak_dir, bedtools_path, filenames_filt, variants_bed_dir, num_peaks_l, id_factor_dict, gt):
    """

    :param bedtools_path: path to bedtools
    :param filenames_filt: filtered list of ENCODE narrowPeak filenames (extension excluded)
    :param variants_bed_dir: path to directory containing BED files of variant coordinates
    :param num_peaks_l: list containing total number of peaks per ENCODE dataset
    :param gt: maternal or paternal
    :return:
    """
    num_instances_per_dataset = []  # nested list, each item [<number of variant overlap instances>, <total number of peaks>]
    num_per_variant_overlap_per_dataset = []    # nested list, each item [<number of variants overlapping peaks>, <total number of variants overlapping peaks>]
    var_bed_filenames = []
    filenames_filt = [narrowpeak_dir + i for i in filenames_filt]   # absolute paths to narrowPeak files
    if gt == 'm':
        var_bed_filenames = ['snps_m.bed', 'ins_m.bed', 'del_m.bed', 'dups_m.bed', 'inv_m.bed']
        var_bed_filenames = [variants_bed_dir + i for i in var_bed_filenames]
    elif gt == 'p':
        var_bed_filenames = ['snps_p.bed', 'ins_p.bed', 'del_p.bed', 'dups_p.bed', 'inv_p.bed']
        var_bed_filenames = [variants_bed_dir + i for i in var_bed_filenames]

    snps, ins, dels, dups, inv \
        = var_bed_filenames[0], var_bed_filenames[1], var_bed_filenames[2], var_bed_filenames[3], var_bed_filenames[4]

    for narrowPeak, num_total_peaks in zip(filenames_filt, num_peaks_l):
        # below command uses '-u', since computing number of variant-containing peaks (not number of overlap instances)
        peaks_overlap_all = subprocess.run([bedtools_path, 'intersect', '-a', narrowPeak, '-b', snps, ins, dels, dups, inv, '-u'],
                                               stdout=subprocess.PIPE)
        # below commands use '-wa', since 1) computing number of overlap instances and 2) for ease of understanding (if checking needed)
        snp_overlap = subprocess.run([bedtools_path, 'intersect', '-a', narrowPeak, '-b', snps, '-wa'], stdout=subprocess.PIPE)
        ins_overlap = subprocess.run([bedtools_path, 'intersect', '-a', narrowPeak, '-b', ins, '-wa'], stdout=subprocess.PIPE)
        dels_overlap = subprocess.run([bedtools_path, 'intersect', '-a', narrowPeak, '-b', dels, '-wa'], stdout=subprocess.PIPE)
        dups_overlap = subprocess.run([bedtools_path, 'intersect', '-a', narrowPeak, '-b', dups, '-wa'], stdout=subprocess.PIPE)
        inv_overlap = subprocess.run([bedtools_path, 'intersect', '-a', narrowPeak, '-b', inv, '-wa'], stdout=subprocess.PIPE)

        # process subprocess output
        num_peaks_overlap_decoded = peaks_overlap_all.stdout.decode('utf-8').split('\n')
        num_snp_overlap_decoded = snp_overlap.stdout.decode('utf-8').split('\n')
        num_ins_overlap_decoded = ins_overlap.stdout.decode('utf-8').split('\n')
        num_dels_overlap_decoded = dels_overlap.stdout.decode('utf-8').split('\n')
        num_dups_overlap_decoded = dups_overlap.stdout.decode('utf-8').split('\n')
        num_inv_overlap_decoded = inv_overlap.stdout.decode('utf-8').split('\n')

        # count total number of variant-containing peaks
        num_peaks_overlap_all = len([i for i in num_peaks_overlap_decoded if 'chr' in i])
        num_snp_overlap = len([i for i in num_snp_overlap_decoded if 'chr' in i])
        num_ins_overlap = len([i for i in num_ins_overlap_decoded if 'chr' in i])
        num_dels_overlap = len([i for i in num_dels_overlap_decoded if 'chr' in i])
        num_dups_overlap = len([i for i in num_dups_overlap_decoded if 'chr' in i])
        num_inv_overlap = len([i for i in num_inv_overlap_decoded if 'chr' in i])

        # TODO: if a peak contains multiple variants, does bedtools intersect output one line of said peak in output
        # or multiple lines of said peak?
        # I think each instance is one line, so multiple lines of same peak is produced is multiple
        # overlaps w/ variants within that peak

        # Also, percentages of all instances might be more informative than percent of all variant-containing peaks
        # the 'all instances' metric probably contains duplicate peaks

        # e.g. 3 peaks, peak 1 has 2 SNPs, peak 2 has 1 ins / 1 del, and peak 3 has 5 SNPs. It'd be more helpful
        # to know that 7/9 (approx. 78%) of overlap instances are SNPs than 2/3 peaks contain SNPs

        num_instances_per_dataset.append([narrowPeak, id_factor_dict[narrowPeak.split('/')[-1].split('.')[0]], num_peaks_overlap_all, num_total_peaks[2]])
        num_per_variant_overlap_per_dataset.append([narrowPeak, id_factor_dict[narrowPeak.split('/')[-1].split('.')[0]], num_snp_overlap, num_ins_overlap, num_dels_overlap, num_dups_overlap, num_inv_overlap])
        print('Overlap statistics computed for %s' % narrowPeak)

    return num_instances_per_dataset, num_per_variant_overlap_per_dataset


def compute_percent_peaks_per_dataset(num_instances_per_dataset):
    """ Using 'num_peaks_overlap_per_dataset' from output of run_bedtools(), compute percentages of
    variant-containing peaks for each ENCODE dataset.
    """
    percentages_l = []
    for i in num_instances_per_dataset:
        filename = i[0].split('/')[-1]
        percentages_l.append([filename, i[1], (i[2] / i[3]) * 100])

    return percentages_l


def compute_percent_overlapping_var_per_dataset(num_per_variant_overlap_per_dataset):
    """ Using 'num_peaks_overlap_per_dataset' from output of run_bedtools(), compute percentages of
    each variant of in variant-containing peaks. run_bedtools() identified the number of variant-containing
    peaks, and compute_percent_peaks_per_dataset() computes the percentages of said peaks as a fraction of
    the total number of peaks (for each dataset). This function asks: of the peaks that contain variants,
    what percent of those peaks contain SNPs, insertions, deletions, duplications, and inversions?

    Input is from 'num_per_variant_overlap_per_dataset' output of run_bedtools()
    each item from said output is:
    [num_snp_overlap, num_ins_overlap, num_dels_overlap, num_dups_overlap, num_inv_overlap]
    """
    percent_per_var = []    # nested list, each item [<percent of (var1)-containing peaks>, <percent of (var2)-containing peaks>, etc.] (for one dataset)

    for i in num_per_variant_overlap_per_dataset:
        total_instances = sum(i[2:])    # total number of variant-containing peaks for ENCODE dataset
        num_instances_snps, num_instaces_ins, num_instances_del, \
        num_instances_dups, num_peaks_inv = i[2], i[3], i[4], i[5], i[6]
        file_id, factor_name = i[0], i[1]       # TODO: check that factor_name is factor name

        percent_snps = (num_instances_snps / total_instances) * 100 # percent of variant-containing peaks containing SNPs
        percent_ins = (num_instaces_ins / total_instances) * 100
        percent_del = (num_instances_del / total_instances) * 100
        percent_dups = (num_instances_dups / total_instances) * 100
        percent_inv = (num_peaks_inv / total_instances) * 100

        percent_per_var.append([file_id, factor_name, percent_snps, percent_ins, percent_del, percent_dups,
                                percent_inv])

    return percent_per_var


def compute_num_var_per_peak(filenames_filt, bedtools_path, narrowpeak_dir, variants_bed_dir, gt):
    """ Compute number of variants per peak for one dataset, output boxplot?
    """
    var_bed_filenames = []
    num_var_in_peaks_counts_l = []
    if gt == 'm':
        var_bed_filenames = ['snps_m.bed', 'ins_m.bed', 'del_m.bed', 'dups_m.bed', 'inv_m.bed']
        var_bed_filenames = [variants_bed_dir + i for i in var_bed_filenames]
    elif gt == 'p':
        var_bed_filenames = ['snps_p.bed', 'ins_p.bed', 'del_p.bed', 'dups_p.bed', 'inv_p.bed']
        var_bed_filenames = [variants_bed_dir + i for i in var_bed_filenames]

    snps, ins, dels, dups, inv \
        = var_bed_filenames[0], var_bed_filenames[1], var_bed_filenames[2], var_bed_filenames[3], var_bed_filenames[4]

    for narrowPeak in filenames_filt:
        peaks_overlap_all = subprocess.run([bedtools_path, 'intersect', '-a', narrowpeak_dir + narrowPeak, '-b', snps, ins, dels, dups, inv, '-wa'],
                                               stdout=subprocess.PIPE)
        peaks_overlap_decoded = peaks_overlap_all.stdout.decode('utf-8').split('\n')
        peaks_overlap_l = [i for i in peaks_overlap_decoded if 'chr' in i]
        peaks_overlap_df = pd.DataFrame(peaks_overlap_l)
        num_var_in_peaks_counts = list(peaks_overlap_df.value_counts())

        num_var_in_peaks_counts_l.extend(num_var_in_peaks_counts)
        print('Counted number of peak-overlapping variants in %s' % narrowPeak)

    # TODO: does it make sense to use raw coordinates (for biological interpretation) since
    # some variants overlap?

    return sorted(num_var_in_peaks_counts_l, reverse=True)  # each item is number of overlapping variants for a peak (from any dataset)


def plot_num_var_per_peak_all_datasets(num_var_per_peak_l, output_plot_dir, gt):
    plt.bar([i for i in range(len(num_var_per_peak_l))], num_var_per_peak_l)
    plt.suptitle('Number of peak-variant overlaps (per peak) for all peaks in all datasets')
    plt.xlabel('Peaks')
    plt.ylabel('Number of peak-variant overlaps')
    plt.savefig(output_plot_dir + 'num_var_in_peaks_counts_sorted_%s.png' % gt)


def compute_snp_positions(narrowpeak_dir, bedtools_path, filenames_filt, var_bed_dir, gt):
    """ Compute SNP positions relative to peak center for one dataset. Just for w/ SNPs for now, not sure about
    implementation for other variants.

    1) Use '-wa' and '-wb' options for bedtools intersect, as described in the manual (https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html).
    This tells us which SNP overlaps with which peak in that dataset.
    2) Since many peaks have different widths, the metric to use is a width-normalized offset from peak center.
    e.g. if peak width is 100 bp and a SNP is 10bp downstream of peak center, the value is 10 / 100 = 0.1
    If said SNP is 10bp upstream of peak center, the value is - (10 / 100) = -0.1, the negative indicating
    upstream of peak center to attribute 'directionality' to SNP position.
    3) Based on computations from 2), all SNP-containing peaks for a given dataset should have 'offset scores'
    ranging from ~ -0.5 to ~ 0.5, with 0 being perfectly at peak center, 0.5 being at the downstream boundary
    of the peak, and ~ -0.5 being at the upstream boundary of the peak.
    4) For a given dataset (not sure how to scale for all datasets yet), plot these 'offset scores' on a
    density plot, maybe using seaborn? For visualization, the density plot can be matched with
    a generalized "peak" (maybe a Gaussian distribution; not actual peak from experiment) going from -0.5 to 0.5

    Process one dataset first. For all datasets, maybe use a heatmap, where each row is a dataset, and each the number
    of entries per row = number of bins from ~ -0.5 to ~ 0.5

    * Asymmetric peaks will have values < -0.5 or > 0.5
    """
    offset_scores_l = []     # nested list, each item contains offset scores for a dataset
    snps = 0
    if gt == 'm':
        snps = 'snps_m.bed'
    elif gt == 'p':
        snps = 'snps_p.bed'

    for narrowPeak in filenames_filt:
        offset_scores_dataset = []
        snp_overlap_instances = subprocess.run([bedtools_path, 'intersect', '-a', narrowpeak_dir + narrowPeak, '-b', var_bed_dir + snps, '-wa', '-wb',
                                                '-filenames'], stdout=subprocess.PIPE)
        snp_overlap_instances_decoded = snp_overlap_instances.stdout.decode('utf-8').split('\n')
        snp_overlap_instances_decoded = [i for i in snp_overlap_instances_decoded if 'chr' in i]
        snp_overlap_instances_decoded_l = [i.split('\t') for i in snp_overlap_instances_decoded]
        snp_overlap_instances_decoded_df = pd.DataFrame(snp_overlap_instances_decoded_l)
        num_overlap_instances = snp_overlap_instances_decoded_df.shape[0]
        for i in range(num_overlap_instances):
            row_i = snp_overlap_instances_decoded_df.iloc[i, :]
            chr_i_narrowpeak, chr_i_snp = row_i[0], row_i[10]
            if chr_i_narrowpeak == chr_i_snp:       # just in case
                peak_sta, peak_end = int(row_i[1]), int(row_i[2])   # don't think int(row_i[2]) - 1 is needed
                snp_coor = int(row_i[11])
                peak_width = peak_end - peak_sta
                offset_from_sta = int(row_i[9])
                peak_summit = peak_sta + offset_from_sta
                offset_score = (snp_coor - peak_summit) / peak_width
                offset_scores_dataset.append(offset_score)
            else:
                raise ValueError('NarrowPeak peak chromosome does not match SNP chromosome for overlap entry %d!' % i)

        offset_scores_l.append([narrowPeak, offset_scores_dataset])
        print('Offset scores computed for %s' % narrowPeak)

    return offset_scores_l


def parse_offset_scores(offset_scores_l, counts_cutoff=100):
    datasets_scores_keep_l = []
    for i in offset_scores_l:
        dataset, scores = i[0], i[1]
        scores_center_region = [i for i in scores if -0.01 <= i <= 0.01]
        if len(scores_center_region) > counts_cutoff:
            datasets_scores_keep_l.append([dataset, scores])

    # top_datasets = [i for i in offset_scores_l if len(i[1]) > num_scores]
    return datasets_scores_keep_l


def plot_top_samples_offset_scores(scores_l, output_dir, gt):
    """ Generate distribution plots for select datasets
    """
    if gt == 'm':
        for i in scores_l:
            sample_filename, scores = i[0], i[1]
            plt_scores = sns.displot(np.array(scores), binwidth=0.01, kde=True)
            plt_scores.savefig(output_dir + '%s_offset_scores_m.png' % sample_filename)
    elif gt == 'p':
        for i in scores_l:
            sample_filename, scores = i[0], i[1]
            plt_scores = sns.displot(np.array(scores), binwidth=0.01, kde=True)
            plt_scores.savefig(output_dir + '%s_offset_scores_p.png' % sample_filename)

    print('Offset scores plot(s) for %s saved to %s' % (gt, output_dir))


def load_edit_chromhmm(chromhmm_data_path, chromhmm_subset_bed_file):
    """ Get chromHMM data and subset for the BED and annotation parts (first 4 columns)

    Note: I think the chromHMM BED file(s) use the UCSC conventions (0-based start, 1-based end)
    """
    with gzip.open(chromhmm_data_path, 'r') as chromhmm:
        chromhmm_lines = chromhmm.readlines()
        chromhmm_lines_decoded = [i.decode('utf-8').rstrip().split('\t') for i in chromhmm_lines]
        chromhmm_lines_decoded_df = pd.DataFrame(chromhmm_lines_decoded)
        chromhmm_lines_decoded_df_subset = chromhmm_lines_decoded_df.iloc[:, :4]

    with open(chromhmm_subset_bed_file, 'w') as chromhmm_bed_subset:
        num_rows_chromhmm = chromhmm_lines_decoded_df_subset.shape[0]
        for i in range(num_rows_chromhmm):
            row_i = chromhmm_lines_decoded_df_subset.iloc[i, :]
            chr_i, sta_i, end_i, annotation = row_i[0], int(row_i[1]), int(row_i[2]), row_i[3]
            if (i + 1) != num_rows_chromhmm:
                chromhmm_bed_subset.write('%s\t%d\t%d\t%s\n' % (chr_i, sta_i, end_i, annotation))
            elif (i + 1) == num_rows_chromhmm:
                chromhmm_bed_subset.write('%s\t%d\t%d\t%s' % (chr_i, sta_i, end_i, annotation))


def bed_var_in_peaks(narrowpeak_dir, bedtools_path, filenames_filt, variants_bed_dir, num_peaks_l, outdir, gt):
    var_bed_filenames = []
    snp_in_peak_coor_df_l = []
    ins_in_peak_coor_df_l = []
    dels_in_peak_coor_df_l = []
    dups_in_peak_coor_df_l = []
    inv_in_peak_coor_df_l = []
    filenames_filt = [narrowpeak_dir + i for i in filenames_filt]   # absolute paths to narrowPeak files
    if gt == 'm':
        var_bed_filenames = ['snps_m.bed', 'ins_m.bed', 'del_m.bed', 'dups_m.bed', 'inv_m.bed']
        var_bed_filenames = [variants_bed_dir + i for i in var_bed_filenames]
    elif gt == 'p':
        var_bed_filenames = ['snps_p.bed', 'ins_p.bed', 'del_p.bed', 'dups_p.bed', 'inv_p.bed']
        var_bed_filenames = [variants_bed_dir + i for i in var_bed_filenames]

    snps, ins, dels, dups, inv \
        = var_bed_filenames[0], var_bed_filenames[1], var_bed_filenames[2], var_bed_filenames[3], var_bed_filenames[4]

    for narrowPeak, num_total_peaks in zip(filenames_filt, num_peaks_l):
        snp_overlap = subprocess.run([bedtools_path, 'intersect', '-a', narrowPeak, '-b', snps, '-wa', '-wb'], stdout=subprocess.PIPE)
        ins_overlap = subprocess.run([bedtools_path, 'intersect', '-a', narrowPeak, '-b', ins, '-wa', '-wb'], stdout=subprocess.PIPE)
        dels_overlap = subprocess.run([bedtools_path, 'intersect', '-a', narrowPeak, '-b', dels, '-wa', '-wb'], stdout=subprocess.PIPE)
        dups_overlap = subprocess.run([bedtools_path, 'intersect', '-a', narrowPeak, '-b', dups, '-wa', '-wb'], stdout=subprocess.PIPE)
        inv_overlap = subprocess.run([bedtools_path, 'intersect', '-a', narrowPeak, '-b', inv, '-wa', '-wb'], stdout=subprocess.PIPE)

        snp_overlap_decoded = snp_overlap.stdout.decode('utf-8').split('\n')
        ins_overlap_decoded = ins_overlap.stdout.decode('utf-8').split('\n')
        dels_overlap_decoded = dels_overlap.stdout.decode('utf-8').split('\n')
        dups_overlap_decoded = dups_overlap.stdout.decode('utf-8').split('\n')
        inv_overlap_decoded = inv_overlap.stdout.decode('utf-8').split('\n')

        snp_in_peak_coor = [i.split('\t')[-3:] for i in snp_overlap_decoded if 'chr' in i]
        ins_in_peak_coor = [i.split('\t')[-3:] for i in ins_overlap_decoded if 'chr' in i]
        dels_in_peak_coor = [i.split('\t')[-3:] for i in dels_overlap_decoded if 'chr' in i]
        dups_in_peak_coor = [i.split('\t')[-3:] for i in dups_overlap_decoded if 'chr' in i]
        inv_in_peak_coor = [i.split('\t')[-3:] for i in inv_overlap_decoded if 'chr' in i]

        snp_in_peak_coor = list(snp_in_peak_coor)
        ins_in_peak_coor = list(ins_in_peak_coor)
        dels_in_peak_coor = list(dels_in_peak_coor)
        dups_in_peak_coor = list(dups_in_peak_coor)
        inv_in_peak_coor = list(inv_in_peak_coor)

        snp_in_peak_coor_df_l.append(pd.DataFrame(snp_in_peak_coor))
        ins_in_peak_coor_df_l.append(pd.DataFrame(ins_in_peak_coor))
        dels_in_peak_coor_df_l.append(pd.DataFrame(dels_in_peak_coor))
        dups_in_peak_coor_df_l.append(pd.DataFrame(dups_in_peak_coor))
        inv_in_peak_coor_df_l.append(pd.DataFrame(inv_in_peak_coor))

        print('In-peak variants identified for %s' % narrowPeak)

    # TODO: double check fmt parameter for other np.savetxt() commands in this program (and in other programs if exists)
    np.savetxt(outdir + 'snps_in_peak_%s.bed' % gt, np.array(pd.concat(snp_in_peak_coor_df_l).drop_duplicates().reset_index(drop=True)), fmt='%s', delimiter='\t')
    np.savetxt(outdir + 'ins_in_peak_%s.bed' % gt, np.array(pd.concat(ins_in_peak_coor_df_l).drop_duplicates().reset_index(drop=True)), fmt='%s', delimiter='\t')
    np.savetxt(outdir + 'del_in_peak_%s.bed' % gt, np.array(pd.concat(dels_in_peak_coor_df_l).drop_duplicates().reset_index(drop=True)), fmt='%s', delimiter='\t')
    np.savetxt(outdir + 'dups_in_peak_%s.bed' % gt, np.array(pd.concat(dups_in_peak_coor_df_l).drop_duplicates().reset_index(drop=True)), fmt='%s', delimiter='\t')
    np.savetxt(outdir + 'inv_in_peak_%s.bed' % gt, np.array(pd.concat(inv_in_peak_coor_df_l).drop_duplicates().reset_index(drop=True)), fmt='%s', delimiter='\t')


def classify_var_chromhmm(bedtools_path, chromhmm_bed_subset_path, chromhmm_segway_bed_subset_path, variants_bed_dir, gt, in_peaks=True, chromhmm_segway_combined=True):
    """ Overlap variant coordinates with chromHMM annotations

    Note: Large variants might span multiple chromHMM states, count all overlapping states (I don't see
    a reason to choose just one, and I think that's what the bedtools command in this function does)

    States listed in chromhmm_annot are from Figure 1b of https://www.nature.com/articles/nprot.2017.124

    variants_bed_dir could be in-peak variants or not
    """
    chromhmm_annot = ['TssA', 'TssFlnk', 'TssFlnkU', 'TssFlnkD', 'Tx', 'TxWk',
                      'EnhG1', 'EnhG2', 'EnhA1', 'EnhA2', 'EnhWk', 'ZNF/Rpts',
                      'Het', 'TssBiv', 'EnhBiv', 'ReprPC', 'ReprPCWk', 'Quies']
    chromhmm_segway_annot = ['TSS', 'PF', 'E', 'WE', 'CTCF', 'T', 'R']

    counts_states = []
    var_bed_filenames = []
    if gt == 'm':
        if in_peaks is True:
            var_bed_filenames = ['snps_in_peak_m.bed', 'ins_in_peak_m.bed', 'del_in_peak_m.bed', 'dups_in_peak_m.bed', 'inv_in_peak_m.bed']
            var_bed_filenames = [variants_bed_dir + i for i in var_bed_filenames]
        else:
            var_bed_filenames = ['snps_m.bed', 'ins_m.bed', 'del_m.bed', 'dups_m.bed', 'inv_m.bed']
            var_bed_filenames = [variants_bed_dir + i for i in var_bed_filenames]
    elif gt == 'p':
        if in_peaks is True:
            var_bed_filenames = ['snps_in_peak_p.bed', 'ins_in_peak_p.bed', 'del_in_peak_p.bed', 'dups_in_peak_p.bed', 'inv_in_peak_p.bed']
            var_bed_filenames = [variants_bed_dir + i for i in var_bed_filenames]
        else:
            var_bed_filenames = ['snps_p.bed', 'ins_p.bed', 'del_p.bed', 'dups_p.bed', 'inv_p.bed']
            var_bed_filenames = [variants_bed_dir + i for i in var_bed_filenames]

    snps, ins, dels, dups, inv \
        = var_bed_filenames[0], var_bed_filenames[1], var_bed_filenames[2], var_bed_filenames[3], var_bed_filenames[4]


    # Overlap with chromHMM annotations
    if chromhmm_segway_combined is True:
        peaks_overlap_all = subprocess.run([bedtools_path, 'intersect', '-a', chromhmm_segway_bed_subset_path, '-b', snps, ins, dels, dups, inv, '-wa'],
                                               stdout=subprocess.PIPE)
    else:
        peaks_overlap_all = subprocess.run([bedtools_path, 'intersect', '-a', chromhmm_bed_subset_path, '-b', snps, ins, dels, dups, inv, '-wa'],
                                               stdout=subprocess.PIPE)
    peaks_overlap_decoded = peaks_overlap_all.stdout.decode('utf-8').split('\n')
    peaks_overlap_l = [i for i in peaks_overlap_decoded if 'chr' in i]
    peaks_overlap_l = [i.split('\t') for i in peaks_overlap_l]
    states_l = [i[3] for i in peaks_overlap_l]

    if chromhmm_segway_combined is True:
        for annot in chromhmm_segway_annot:
            counts_states.append([annot, states_l.count(annot)])
    else:
        for annot in chromhmm_annot:
            counts_states.append([annot, states_l.count(annot)])

    if in_peaks is True:
        if chromhmm_segway_combined is True:
            print('Chromatin occupancy statistics (chromHMM & segway) computed (in peaks)!')
        else:
            print('Chromatin occupancy statistics (chromHMM) computed (in peaks)!')
    else:
        if chromhmm_segway_combined is True:
            print('Chromatin occupancy statistics (chromHMM & segway) computed (all)!')
        else:
            print('Chromatin occupancy statistics (chromHMM) computed (all)!')

    return counts_states


def chromhmm_instances_percentages(counts_states):
    percentages = []
    all_counts = sum([i[1] for i in counts_states])
    for i in counts_states:
        state, count = i[0], i[1]
        percentages.append([state, (count / all_counts) * 100])

    return percentages


def excel_chromhmm(percentages_states_l, output_file):
    to_excel(pd.DataFrame(percentages_states_l), output_file)



def num_archetypal_motifs_in_peaks(narrowpeak_dir, filenames_l, motifs_file_path, overlapping_motifs_bed_output_dir):
    """ Using motif data from https://www.nature.com/articles/s41586-020-2528-x#MOESM1 and
    https://www.vierstra.org/resources/motif_clustering, query each peak for archetypal motif overlap for each dataset
    """
    num_motifs_overlap_peaks = []

    for file in filenames_l:
        motifs_overlap_all = subprocess.run([bedtools_path, 'intersect', '-a', motifs_file_path, '-b', narrowpeak_dir + file, '-wa'],
                                           stdout=subprocess.PIPE)
        motifs_overlap_decoded = motifs_overlap_all.stdout.decode('utf-8').split('\n')
        motifs_overlap_l = [i for i in motifs_overlap_decoded if 'chr' in i]
        motifs_overlap_l = [i.split('\t') for i in motifs_overlap_l]
        num_motifs_overlap = len(motifs_overlap_l)
        motifs_overlap_l_arr = np.array(motifs_overlap_l)
        np.savetxt(overlapping_motifs_bed_output_dir + 'motifs_overlap_peak_%s.bed' % file, motifs_overlap_l_arr, fmt='%s', delimiter='\t')

        num_motifs_overlap_peaks.append(num_motifs_overlap)

    return num_motifs_overlap_peaks


def num_var_overlap_motifs(overlapping_motifs_bed_output_dir, variants_bed_dir, gt):
    """ Compute total number of variants overlapping / within motifs, and number of each variant type
    overlapping / within motifs

    Directory 'overlapping_motifs_bed_output_dir' contains all BED files with motifs overlapping / within peaks

    2 possible approaches (taking former)

    1) Intersect all variants (using all variant BED files) with peak-overlapping motif BED files
    2) Generate a peak-overlapping variant BED file for each narrowPeak and use those for motif-variant overlaps
    for every narrowPeak file
    """
    var_bed_filenames = []
    percent_motif_overlap_variants_in_peaks = []
    percent_motif_overlap_variants_in_peaks_per_type = []
    if gt == 'm':
        var_bed_filenames = ['snps_m.bed', 'ins_m.bed', 'del_m.bed', 'dups_m.bed', 'inv_m.bed']
        var_bed_filenames = [variants_bed_dir + i for i in var_bed_filenames]
    elif gt == 'p':
        var_bed_filenames = ['snps_p.bed', 'ins_p.bed', 'del_p.bed', 'dups_p.bed', 'inv_p.bed']
        var_bed_filenames = [variants_bed_dir + i for i in var_bed_filenames]

    snps, ins, dels, dups, inv \
        = var_bed_filenames[0], var_bed_filenames[1], var_bed_filenames[2], var_bed_filenames[3], var_bed_filenames[4]

    peaks_overlapping_motif_bed_files = subprocess.run('ls', stdout=subprocess.PIPE, cwd=overlapping_motifs_bed_output_dir)
    list_bed_l = peaks_overlapping_motif_bed_files.stdout.decode('utf-8').split('\n')        # files list
    list_motif_bed_filenames = [i for i in list_bed_l if 'motifs_overlap_peak_' in i]     # incl. extension

    for motif_bed in list_motif_bed_filenames:
        motif_file_path = overlapping_motifs_bed_output_dir + motif_bed

        # number of motifs that overlap with peaks for dataset
        with open(motif_file_path, 'r') as motif_file_i:
            num_motifs_overlap_peaks_i = len([i for i in motif_file_path.readlines() if 'chr' in i])

        # overlaps with any variant
        motifs_overlap_var_all = subprocess.run([bedtools_path, 'intersect', '-a', motif_file_path, '-b', snps, ins, dels, dups, inv, '-wa'],
                                               stdout=subprocess.PIPE)
        motifs_overlap_var_decoded = motifs_overlap_var_all.stdout.decode('utf-8').split('\n')
        motifs_overlap_var_l = [i for i in motifs_overlap_var_decoded if 'chr' in i]
        motifs_overlap_var_l = [i.split('\t') for i in motifs_overlap_var_l]
        percent_motif_overlap_variants_in_peaks.append((len(motifs_overlap_var_l) / num_motifs_overlap_peaks_i) * 100)

        # overlaps with each variant type
        motifs_snp_overlap = subprocess.run([bedtools_path, 'intersect', '-a', motif_file_path, '-b', snps, '-wa'], stdout=subprocess.PIPE)
        motifs_ins_overlap = subprocess.run([bedtools_path, 'intersect', '-a', motif_file_path, '-b', ins, '-wa'], stdout=subprocess.PIPE)
        motifs_dels_overlap = subprocess.run([bedtools_path, 'intersect', '-a', motif_file_path, '-b', dels, '-wa'], stdout=subprocess.PIPE)
        motifs_dups_overlap = subprocess.run([bedtools_path, 'intersect', '-a', motif_file_path, '-b', dups, '-wa'], stdout=subprocess.PIPE)
        motifs_inv_overlap = subprocess.run([bedtools_path, 'intersect', '-a', motif_file_path, '-b', inv, '-wa'], stdout=subprocess.PIPE)

        # process subprocess output
        num_snp_overlap_motif_decoded = motifs_snp_overlap.stdout.decode('utf-8').split('\n')
        num_ins_overlap_motif_decoded = motifs_ins_overlap.stdout.decode('utf-8').split('\n')
        num_dels_overlap_motif_decoded = motifs_dels_overlap.stdout.decode('utf-8').split('\n')
        num_dups_overlap_motif_decoded = motifs_dups_overlap.stdout.decode('utf-8').split('\n')
        num_inv_overlap_motif_decoded = motifs_inv_overlap.stdout.decode('utf-8').split('\n')

        # count total number of variant-containing peaks
        num_snp_overlap_motifs = len([i for i in num_snp_overlap_motif_decoded if 'chr' in i])
        num_ins_overlap_motifs = len([i for i in num_ins_overlap_motif_decoded if 'chr' in i])
        num_dels_overlap_motifs = len([i for i in num_dels_overlap_motif_decoded if 'chr' in i])
        num_dups_overlap_motifs = len([i for i in num_dups_overlap_motif_decoded if 'chr' in i])
        num_inv_overlap_motifs = len([i for i in num_inv_overlap_motif_decoded if 'chr' in i])

        num_var_overlap_motif_l = [num_snp_overlap_motifs,
                                num_ins_overlap_motifs,
                                num_dels_overlap_motifs,
                                num_dups_overlap_motifs,
                                num_inv_overlap_motifs]

        total_num_var_overlap_motif_dataset = sum([num_snp_overlap_motifs,
                                            num_ins_overlap_motifs,
                                            num_dels_overlap_motifs,
                                            num_dups_overlap_motifs,
                                            num_inv_overlap_motifs])

        percentages = [(i / total_num_var_overlap_motif_dataset) * 100 for i in num_var_overlap_motif_l]

        percent_motif_overlap_variants_in_peaks_per_type.append(percentages)

    return percent_motif_overlap_variants_in_peaks_per_type     # every item a list 'percentages', for one dataset (I think)


def generate_ExPecto_snps_bed(expecto_snps_path, cols_path, output_bed_path):
    """ SNPs presumably from GTEx, 1000 Genomes proejct, and GWAS Catalog, from this paper:
    https://www.nature.com/articles/s41588-018-0160-6#Abs1

    Make first 3 columns BED format and output modified file
    """
    # snps_df = pd.read_table(expecto_snps_path)

    with open(cols_path, 'r') as cols:
        cols_lines_l = cols.readlines()[0].rstrip().split('\t')

    snps_df = pd.read_table(expecto_snps_path, names=cols_lines_l)
    coor_0_based = snps_df['pos'] - 1     # I think the provided coordinates are 1-based, based on a docstring from this code (line 149):
    # https://github.com/FunctionLab/ExPecto/blob/master/chromatin.py
    snps_df.insert(1, '0_based_coor', list(coor_0_based))
    np.savetxt(output_bed_path, np.array(snps_df), fmt='%s', delimiter='\t')


def overlap_expecto_snps(bedtools_path, var_bed_dir_path, expecto_bed_path, expecto_cols_path, vcf_paths, col_names_vcf_l, gt):
    """ Note: should check if ExPecto SNP matches VCF SNP (base, not just coordinate)

    :param bedtools_path:
    :param var_bed_dir_path:
    :param expecto_bed_path:
    :param gt:
    :return:
    """
    snps = 0
    if gt == 'm':
        snps = 'snps_m.bed'
    elif gt == 'p':
        snps = 'snps_p.bed'

    snps_overlap_all = subprocess.run([bedtools_path, 'intersect', '-a', expecto_bed_path, '-b', var_bed_dir_path + snps, '-wa', '-wb'],
        stdout=subprocess.PIPE)
    snps_overlap_decoded = snps_overlap_all.stdout.decode('utf-8').split('\n')
    snps_overlap_df = pd.DataFrame([i.split('\t') for i in snps_overlap_decoded if 'chr' in i])

    with open(expecto_cols_path, 'r') as cols:
        cols_lines_l = cols.readlines()[0].rstrip().split('\t')
    cols_lines_l.insert(1, 'sta')

    cols_lines_l = cols_lines_l + ['snp_chr', 'snp_sta', 'snp_end']    # of overlapped region (SNPs from VCF, not from ExPecto)

    snps_overlap_df.columns = cols_lines_l

    k562_only = pd.concat([snps_overlap_df['chr'], snps_overlap_df['sta'], snps_overlap_df['pos'], snps_overlap_df['K562.1'],
                           snps_overlap_df['K562']], axis=1)

    snp_df = get_snp_df(vcf_paths, col_names_vcf_l)
    series_coor_overlaps = pd.Series(snps_overlap_df['sta'], dtype='int64')
    snp_df_filt = snp_df[snp_df['POS'].isin(list(set(series_coor_overlaps)))]
    num_rows_snp_df_filt = snp_df_filt.shape[0]
    # for i in range(num_rows_snp_df_filt):
    #     row_snp_df_filt, row_expecto_snp_df = snp_df_filt[i, :], k562_only[i, :]
    #     expecto_ref, vcf_ref = row_snp_df_filt[3], row_expecto_snp_df[3]
    #     expecto_alt, vcf_alt = row_snp_df_filt[4], row_expecto_snp_df[4]
    #     if expecto_ref == vcf_ref and expecto_alt == vcf_alt:
    #         continue
    #     else:
    #         print('ExPecto ref/alt allele does not match VCF ref/alt allele!')


def overlap_epimap_enhancers(bedtools_path, variants_bed_dir, epimap_bed_path, epimap_file_path, gt):
    """ Identify how many EpiMap-linked enhancers contain variants (all, per type)

    :param bedtools_path:
    :param var_bed_dir_path:
    :param epimap_file_path:
    :return:
    """
    var_bed_filenames = []
    if gt == 'm':
        var_bed_filenames = ['snps_m.bed', 'ins_m.bed', 'del_m.bed', 'dups_m.bed', 'inv_m.bed']
        var_bed_filenames = [variants_bed_dir + i for i in var_bed_filenames]
    elif gt == 'p':
        var_bed_filenames = ['snps_p.bed', 'ins_p.bed', 'del_p.bed', 'dups_p.bed', 'inv_p.bed']
        var_bed_filenames = [variants_bed_dir + i for i in var_bed_filenames]

    snps, ins, dels, dups, inv \
        = var_bed_filenames[0], var_bed_filenames[1], var_bed_filenames[2], var_bed_filenames[3], var_bed_filenames[4]

    with gzip.open(epimap_file_path, 'r') as epimap:
        epimap_lines = epimap.readlines()
        epimap_lines_l = [i.decode('utf-8') for i in epimap_lines]
        epimap_lines_l = [i.rstrip().split('\t') for i in epimap_lines_l]
        np.savetxt(epimap_bed_path, np.array(epimap_lines_l), fmt='%s', delimiter='\t')
        epimap_lines_df = pd.DataFrame(epimap_lines_l)
        epimap_lines_df.columns = ['chr', 'sta', 'end', 'gene', 'corr', 'type']
        epimap_lines_df = epimap_lines_df.drop_duplicates(subset=['chr', 'sta', 'end']).reset_index(drop=True)
        num_enh_w_linked_gene = epimap_lines_df.shape[0]

    epimap_overlap_var_all = subprocess.run(
        [bedtools_path, 'intersect', '-a', epimap_bed_path, '-b', snps, ins, dels, dups, inv, '-u'],
        stdout=subprocess.PIPE)
    epimap_overlap_var_decoded = epimap_overlap_var_all.stdout.decode('utf-8').split('\n')
    epimap_overlap_var_l = [i for i in epimap_overlap_var_decoded if 'chr' in i]
    epimap_overlap_var_df = pd.DataFrame([i.split('\t') for i in epimap_overlap_var_l])
    epimap_overlap_var_df.columns = ['chr', 'sta', 'end', 'gene', 'corr', 'type']

    # some enhancers are linked to multiple genes
    epimap_overlap_var_df = epimap_overlap_var_df.drop_duplicates(subset=['chr', 'sta', 'end']).reset_index(drop=True)

    # number of variants overlapping with enhancers, each with one or more linked genes
    num_var_in_enh = epimap_overlap_var_df.shape[0]

    print('overlap_epimap_enhancers(): Percent of variants overlapping with EpiMap enhancers for %s: %.5f' % (gt, (num_var_in_enh / num_enh_w_linked_gene) * 100))

    # return (num_var_in_enh / num_enh_w_linked_gene) * 100   # return percentage


def overlap_abc_enhancers(bedtools_path, variants_bed_dir, abc_bed_path, abc_file_path, gt):
    """ Identify variants in "ABC enhancers", defined in this paper: https://www.nature.com/articles/s41586-021-03446-x

    Some entries in ABC file doesn't seem to be biological enhancers, just element-gene links
    """
    var_bed_filenames = []
    if gt == 'm':
        var_bed_filenames = ['snps_m.bed', 'ins_m.bed', 'del_m.bed', 'dups_m.bed', 'inv_m.bed']
        var_bed_filenames = [variants_bed_dir + i for i in var_bed_filenames]
    elif gt == 'p':
        var_bed_filenames = ['snps_p.bed', 'ins_p.bed', 'del_p.bed', 'dups_p.bed', 'inv_p.bed']
        var_bed_filenames = [variants_bed_dir + i for i in var_bed_filenames]

    snps, ins, dels, dups, inv \
        = var_bed_filenames[0], var_bed_filenames[1], var_bed_filenames[2], var_bed_filenames[3], var_bed_filenames[4]

    with gzip.open(abc_file_path, 'r') as abc:
        abc_lines = abc.readlines()
        abc_lines_l = [i.decode('utf-8') for i in abc_lines]
        abc_lines_l = [i.rstrip().split('\t') for i in abc_lines_l]
        print('overlap_abc_enhancers() for %s: abc_lines_l generated' % gt)
        np.savetxt(abc_bed_path, np.array(abc_lines_l), fmt='%s', delimiter='\t')
        print('overlap_abc_enhancers() for %s: ABC BED file generated' % gt)
        columns = abc_lines_l[0]
        abc_lines_df = pd.DataFrame(abc_lines_l[1:])
        abc_lines_df.columns = columns
        abc_lines_df = abc_lines_df.drop_duplicates(subset=columns[:3]).reset_index(drop=True)
        num_element_w_linked_gene = abc_lines_df.shape[0]

    print('overlap_abc_enhancers() for %s: Running bedtools intersect' % gt)
    abc_overlap_var_all = subprocess.run(
        [bedtools_path, 'intersect', '-a', abc_bed_path, '-b', snps, ins, dels, dups, inv, '-u'],
        stdout=subprocess.PIPE)
    abc_overlap_var_decoded = abc_overlap_var_all.stdout.decode('utf-8').split('\n')
    abc_overlap_var_l = [i for i in abc_overlap_var_decoded if 'chr' in i]
    abc_overlap_var_df = pd.DataFrame([i.split('\t') for i in abc_overlap_var_l])
    abc_overlap_var_df.columns = columns

    abc_overlap_var_df = abc_overlap_var_df.drop_duplicates(subset=columns[:3]).reset_index(drop=True)

    # number of variants overlapping with enhancers, each with one or more linked genes
    num_var_in_element = abc_overlap_var_df.shape[0]
    percent_var_in_element = (num_var_in_element / num_element_w_linked_gene) * 100

    print('%.5f percent of variants overlapping with ABC enhancers / elements (with linked genes) for %s' % (percent_var_in_element, gt))

    # return (num_var_in_element / num_element_w_linked_gene) * 100   # return percentage


def compare_peaks_custom_hg19(bedtools_path, path_narrowpeak_custom, path_narrowpeak_hg19):
    """ Compare peak occupancy for maternal and paternal custom genome and hg19

    Check overlaps using bedtools intersect:
    1) How many retained peaks?
    2) How many new peaks?
    3) How many lost peaks?

    Would help to visualize on IGV.

    Of the new peaks, how many have signal change, shape change? (can use BPNet to confirm any
    changes in contribution scores after parsing raw data)
    """


def run_functions(dir_narrowpeak, path_metadata, excel_dir_path, var_bed_dir_path, dir_in_peak_var_bed, gt, cols, paths, bedtools_path, chromhmm_bed_subset_path, chromhmm_data_path,
                  num_total_peaks, epimap_bed_path, epimap_file_path, abc_bed_path, abc_file_path, chromhmm_segway_path, chromhmm_bed_segway_subset_path, offset_score_plots_dir_path):
    """

    :param dir_narrowpeak:
    :param coor_l: contains coordinates per variant (for either maternal or paternal)
    :param excel_dir_path: directory to output excel files, including final '/'
    :return:
    """
    keep_id, id_factor_dict = parse_metadata(path_metadata)
    narrowpeak_filenames = filenames(dir_narrowpeak, keep_id)

    print('Generating BED files for each variant for %s...' % gt)
    var_all = get_var_all(paths, cols)
    print('Generating coordinates for each variant for %s...' % gt)
    snps_coor_all, ins_coor_all, del_coor_all, dups_coor_all, inv_coor_all = coor_all_chr(var_all, gt)
    print('Generating coordinate dataframes for each variant for %s...' % gt)
    snps_df_final, ins_df_final, del_df_final, dups_df_final, inv_df_final = gen_bed_df_var(snps_coor_all, ins_coor_all,
                                                                                            del_coor_all, dups_coor_all,
                                                                                            inv_coor_all)
    print('Generating BED files for each variant for %s...' % gt)
    gen_bed_files([snps_df_final, ins_df_final, del_df_final, dups_df_final, inv_df_final],
                  var_bed_dir_path, gt)
    print('BED files generated for each variant for %s' % gt)

    num_peaks_overlap_per_dataset, num_per_variant_overlap_per_dataset = \
        run_bedtools(dir_narrowpeak, bedtools_path, narrowpeak_filenames, var_bed_dir_path, num_total_peaks, id_factor_dict, gt)

    percentages_l = compute_percent_peaks_per_dataset(num_peaks_overlap_per_dataset)
    percent_per_var = compute_percent_overlapping_var_per_dataset(num_per_variant_overlap_per_dataset)
    
    # Variant BED files for in-peak variants
    bed_var_in_peaks(dir_narrowpeak, bedtools_path, narrowpeak_filenames, var_bed_dir_path, num_total_peaks, dir_in_peak_var_bed, gt)
    
    plot_sample_num_peaks_var_in_peaks(num_total_peaks, percentages_l, excel_dir_path)

    if gt == 'm':
        to_excel(pd.DataFrame(percentages_l), excel_dir_path + 'percent_peaks_overlap_per_dataset_m.xlsx')
        to_excel(pd.DataFrame(percent_per_var), excel_dir_path + 'percent_per_variant_overlap_per_dataset_m.xlsx')
        print('Excel files (maternal) for num_peaks_overlap_per_dataset and'
              'num_per_variant_overlap_per_dataset generated!')
    elif gt == 'p':
        to_excel(pd.DataFrame(percentages_l), excel_dir_path + 'percent_peaks_overlap_per_dataset_p.xlsx')
        to_excel(pd.DataFrame(percent_per_var), excel_dir_path + 'percent_per_variant_overlap_per_dataset_p.xlsx')
        print('Excel files (paternal) for num_peaks_overlap_per_dataset and'
              'num_per_variant_overlap_per_dataset generated!')


    # SNP position analysis
    offset_scores_l = compute_snp_positions(dir_narrowpeak, bedtools_path, narrowpeak_filenames, var_bed_dir_path, gt)
    enriched_offset_scores_samples = parse_offset_scores(offset_scores_l, counts_cutoff=500)
    
    plot_top_samples_offset_scores(enriched_offset_scores_samples, offset_score_plots_dir_path, gt)
    
    # TODO: visually check values in these lists (> 30 variants in peak is a lot) (*turns out result of including all, incl. overlapping variants)
    num_var_in_peaks_counts_sorted_l = compute_num_var_per_peak(narrowpeak_filenames, bedtools_path, dir_narrowpeak, var_bed_dir_path, gt)
    plot_num_var_per_peak_all_datasets(num_var_in_peaks_counts_sorted_l, excel_dir_path, gt)

    Get subset chromHMM, chromHMM + segway BED files
    load_edit_chromhmm(chromhmm_data_path, chromhmm_bed_subset_path)
    load_edit_chromhmm(chromhmm_segway_path, chromhmm_bed_segway_subset_path)
    
    # All variants (chromHMM)
    counts_states_all_chromhmm = classify_var_chromhmm(bedtools_path, chromhmm_bed_subset_path, chromhmm_bed_segway_subset_path, var_bed_dir_path, gt,
                                              in_peaks=False, chromhmm_segway_combined=False)
    counts_states_in_peaks_chromhmm = classify_var_chromhmm(bedtools_path, chromhmm_bed_subset_path, chromhmm_bed_segway_subset_path, dir_in_peak_var_bed, gt,
                                                   in_peaks=True, chromhmm_segway_combined=False)
    percentages_all_chromhmm = chromhmm_instances_percentages(counts_states_all_chromhmm)
    percentages_in_peaks_chromhmm = chromhmm_instances_percentages(counts_states_in_peaks_chromhmm)
    
    excel_chromhmm(percentages_all_chromhmm, excel_dir_path + 'chromhmm_percentages_all_%s.xlsx' % gt)
    
    # Only variants in peaks (chromHMM)
    excel_chromhmm(percentages_in_peaks_chromhmm, excel_dir_path + 'chromhmm_percentages_in_peaks_%s.xlsx' % gt)

    # All variants (chromHMM + segway)
    counts_states_all_chromhmm_segway = classify_var_chromhmm(bedtools_path, chromhmm_bed_subset_path, chromhmm_bed_segway_subset_path, var_bed_dir_path, gt,
                                              in_peaks=False, chromhmm_segway_combined=True)
    counts_states_in_peaks_chromhmm_segway = classify_var_chromhmm(bedtools_path, chromhmm_bed_subset_path, chromhmm_bed_segway_subset_path, dir_in_peak_var_bed, gt,
                                                   in_peaks=True, chromhmm_segway_combined=True)
    percentages_all_chromhmm_segway = chromhmm_instances_percentages(counts_states_all_chromhmm_segway)
    percentages_in_peaks_chromhmm_segway = chromhmm_instances_percentages(counts_states_in_peaks_chromhmm_segway)

    excel_chromhmm(percentages_all_chromhmm_segway, excel_dir_path + 'chromhmm_segway_percentages_all_%s.xlsx' % gt)

    # Only variants in peaks (chromHMM + segway)
    excel_chromhmm(percentages_in_peaks_chromhmm_segway, excel_dir_path + 'chromhmm_segway_percentages_in_peaks_%s.xlsx' % gt)


    # Motif archetype analysis
    num_archetypal_motifs_in_peaks(dir_narrowpeak, narrowpeak_filenames,
                                   motifs_file_path, overlapping_motifs_bed_output_dir)
    
    num_var_overlap_motifs(overlapping_motifs_bed_output_dir,
                           var_bed_dir_path, gt)

    # EpiMap
    overlap_epimap_enhancers(bedtools_path,
                             var_bed_dir_path,
                             epimap_bed_path,
                             epimap_file_path,
                             gt)

    # ABC
    overlap_abc_enhancers(bedtools_path,
                          var_bed_dir_path,
                          abc_bed_path,
                          abc_file_path,
                          gt)



def run_all():
    vcf_dir = '/path/to/vcf/'
    f1, f2, f3, f4 = 'ENCFF574MDJ_varied.vcf', 'ENCFF752OAX_SNP.vcf', 'ENCFF785JVR_DEL.vcf', 'ENCFF863MPP_BND.vcf'
    paths = [vcf_dir + f1, vcf_dir + f2, vcf_dir + f3, vcf_dir + f4]
    cols = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'K562']
    bed_files_dir = '/path/to/bed/'
    excel_files_dir = '/path/to/excel_output/'
    metadata_path = '/path/to/metadata.tsv'
    var_bed_dir_path = '/path/to/variants_bed_files/'
    bedtools_path = '/path/to/bedtools'
    chromhmm_data_path = '/path/to/BSS00762_18_CALLS_segments.bed.gz'
    chromhmm_bed_subset_path = '/path/to/chromhmm_bed_subset.bed'
    chromhmm_segway_path = '/path/to/wgEncodeAwgSegmentationCombinedK562.bed.gz'
    chromhmm_bed_segway_subset_path = '/path/to/chromhmm_segway_bed_subset.bed'
    # expecto_snps_cols_path = '/path/to/colnames.final'
    # expecto_snps_path = '/path/to/combined_snps.0.3.final.sorted'
    # expecto_cols_path = '/path/to/colnames.final'
    epimap_file_path = '/path/to/BSS00762_collated_pred.tsv.gz'
    epimap_bed_path = '/path/to/k562_epimap.bed'
    abc_bed_path = '/path/to/abc.bed'
    abc_file_path = '/path/to/AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz'
    dir_in_peak_var_bed = '/path/to/variants_in_peaks_bed_files/'
    offset_score_plots_dir_path = '/path/to/offset_scores_plots/'
    # motifs_file_path = 'NA'
    # overlapping_motifs_bed_output_dir = 'NA'

    num_total_peaks = count_peaks(metadata_path, bed_files_dir, excel_files_dir)

    print('Running functions for maternal genome')
    run_functions(bed_files_dir, metadata_path, excel_files_dir, var_bed_dir_path, dir_in_peak_var_bed, 'm', cols, paths, bedtools_path, chromhmm_bed_subset_path, chromhmm_data_path,
                  num_total_peaks, epimap_bed_path, epimap_file_path, abc_bed_path, abc_file_path, chromhmm_segway_path, chromhmm_bed_segway_subset_path, offset_score_plots_dir_path)

    print('Running functions for paternal genome')
    run_functions(bed_files_dir, metadata_path, excel_files_dir, var_bed_dir_path, dir_in_peak_var_bed, 'p', cols, paths, bedtools_path, chromhmm_bed_subset_path, chromhmm_data_path,
                  num_total_peaks, epimap_bed_path, epimap_file_path, abc_bed_path, abc_file_path, chromhmm_segway_path, chromhmm_bed_segway_subset_path, offset_score_plots_dir_path)

    
# Run the functions
run_all()
