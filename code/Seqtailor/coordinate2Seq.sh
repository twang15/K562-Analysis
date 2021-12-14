# bash coordinate2Seq.sh 10 100 chr1

## 0-based or 1-based? https://www.biostars.org/p/84686/
# 1. the web-interface of the UCSC browser is 1-based, 
#   as it meant to be used by biologists, but all internal representations 
#   (text files, database tables, binary files) are 0-based, as they're mostly used by programmers.
#   Unfortunately, for historical reasons, there two exceptions, two formats 
#   (the wiggle and bigWig text file and database formats) are 1-based.

# 2. While BAM is 0-based, once you pull it out to something human-readable it can get turned into 1-based data (SAM).

begin=$1
end=$2
chrom=$3

wget http://togows.org/api/ucsc/hg38/${chrom}:${begin}-${end} -O sequence.txt &> /dev/null

sequence=`cat sequence.txt`
echo "$sequence"

rm sequence.txt
