# prepare the vcf input
pudb motif_deletion_snps.py

# SeqEnumerator
pudb SeqEnumerator.py -g /oak/stanford/scg/prj_ENCODE/ChIP/basepairmodel/basepairmodels/tutorial/reference/hg38.genome.fa -c 0 -s FORWARD -wd 1056 -wu 1057 -q BOTH --independent yes -i /oak/stanford/scg/prj_ENCODE/ChIP/basepairmodel/basepairmodels/tutorial/Output/CTCF_deletion.vcf -o ins_containing.DNA.fasta -r ins_neighbor_report.report.txt

# verify output:
sed -n "2p" CTCF_deletion.vcf.tsv | awk '{split($9, chars, ""); for (i=1; i<=length($9); i++) {print chars[i]} }' > ref
sed -n "2p" CTCF_deletion.vcf.tsv | awk '{split($10, chars, ""); for (i=1; i<=length($10); i++) {print chars[i]} }' > var
vimdiff ref var

# predict
srun --nodes=1 --ntasks=1 --cpus-per-task=10 --partition=interactive --account=default --time=24:00:00 --pty /bin/bash
source ~/conda_setup
conda activate basepairmodels
pudb predict.py --model Models_from_Anshul/ENCSR000EGM/models/model_split000.h5 -g reference/hg38.genome.fa -s reference/hg38.chrom.sizes -f Output/CTCF_deletion.vcf.tsv -o Output
