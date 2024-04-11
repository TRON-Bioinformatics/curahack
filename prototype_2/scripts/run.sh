!#/bin/bash

export POLARS_MAX_THREADS=4
for v1 in {86..108}
do
    for v2 in {87..109}
    do
        python diff_gtf.py --output output/$v2.vs.$v1.csv \
          ../reference/Homo_sapiens.GRCh38.$v2.chr1.gtf.gz \
          ../reference/Homo_sapiens.GRCh38.$v1.chr1.gtf.gz > output/$v2.vs.$v1.log
    done
done