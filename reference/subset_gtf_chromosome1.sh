
gtf=$1

zcat $gtf | grep -e '^#'  > ${gtf%.gtf.gz}.chr1.gtf
zcat $gtf | grep -P '^1\t'  >> ${gtf%.gtf.gz}.chr1.gtf