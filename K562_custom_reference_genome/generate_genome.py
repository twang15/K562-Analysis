""" CLI to run k562_custom_genome.py
"""
import k562_custom_genome as k
import argparse

# Program summary
parser = argparse.ArgumentParser(description='Generate maternal and paternal custom reference genomes (in FASTA format) '
                                             'from VCF files. Written for K562 at this time.')

# Required arguments
parser.add_argument('vcf1', help='First input VCF filename (incl. extension) / file path')
parser.add_argument('vcf2', help='Second input VCF filename (incl. extension) / file path')
parser.add_argument('vcf3', help='Third input VCF filename (incl. extension) / file path')
parser.add_argument('vcf4', help='Fourth input VCF filename (incl. extension) / file path')

parser.add_argument('reference_path', help='Directory containing all chromosomal reference FASTA files (i.e. chr1.fa, '
                                           'chr2.fa, etc.). Only hg19 supported at this time.')
parser.add_argument('output_dir', help='Directory to store output FASTA files')
args = parser.parse_args()

# Inputs to run_functions_dict() in vcf_encode.py
fa_dir = args.reference_path
outdir = args.output_dir
paths_vcf = [args.vcf1, args.vcf2, args.vcf3, args.vcf4]

# Generate maternal custom genome
print('\n=================================\n')
print('Generating maternal custom genome')
k.run_functions_dict(paths_vcf, fa_dir, 'm', outdir)

# Generate paternal custom genome
print('\n=================================\n')
print('Generating paternal custom genome')
k.run_functions_dict(paths_vcf, fa_dir, 'p', outdir)

