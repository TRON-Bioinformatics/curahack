````
##################
### gtf to bed ###
##################

### Human ###
# 86
awk -v OFS="\t" '$7 == "+" {print $1, $4, $5, $10, $14}' Homo_sapiens.GRCh38.86_chr.gtf \
| sed 's/"//g' | sed 's/ //g' | sed 's/;//g' \
| grep -hv 'ENST' \
| sort -k 6,6 -k2,2n \
| awk -v OFS="\t" '{print $1, $2, $3, $5, $6, "+"}' | awk '{if ($2>0) print $0}' | grep -hv 'chrRP*' | grep -hv 'chrKI*' | grep -hv 'chrGL*' | sortBed -i - > genes_86_plus.bed

awk -v OFS="\t" '$7 == "-" {print $1, $4, $5, $10, $14}' Homo_sapiens.GRCh38.86_chr.gtf \
| sed 's/"//g' | sed 's/ //g' | sed 's/;//g' \
| grep -hv 'ENST' \
| sort -k 6,6 -k2,2n \
| awk -v OFS="\t" '{print $1, $2, $3, $5, $6, "-"}' | awk '{if ($2>0) print $0}' | grep -hv 'chrRP*' | grep -hv 'chrKI*' | grep -hv 'chrGL*' | sortBed -i - > genes_86_minus.bed

# merge
cat genes_86_minus.bed genes_86_plus.bed \
| sortBed -i - \
| uniq  > genes_hg38_86_both.bed

# 109
```
awk -v OFS="\t" '$7 == "+" {print $1, $4, $5, $10, $14}'
Homo_sapiens.GRCh38.109_chr.gtf \
| sed 's/"//g' | sed 's/ //g' | sed 's/;//g' \
| grep -hv 'ENST' \
| sort -k 6,6 -k2,2n \
| awk -v OFS="\t" '{print $1, $2, $3, $5, $6, "+"}' | awk '{if ($2>0) print $0}' | grep -hv 'chrRP*' | grep -hv 'chrKI*' | grep -hv 'chrGL*' | sortBed -i - > genes_109_plus.bed

awk -v OFS="\t" '$7 == "-" {print $1, $4, $5, $10, $14}' Homo_sapiens.GRCh38.109_chr.gtf \
| sed 's/"//g' | sed 's/ //g' | sed 's/;//g' \
| grep -hv 'ENST' \
| sort -k 6,6 -k2,2n \
| awk -v OFS="\t" '{print $1, $2, $3, $5, $6, "-"}' | awk '{if ($2>0) print $0}' | grep -hv 'chrRP*' | grep -hv 'chrKI*' | grep -hv 'chrGL*' | sortBed -i - > genes_109_minus.bed

# merge
cat genes_109_minus.bed genes_109_plus.bed \
| sortBed -i - \
| uniq  > genes_hg38_109_both.bed


####################
### bed to fasta ###
####################
bedtools getfasta -fi genome.fa -bed genes_hg38_86_both.bed -fo hg38_86_gtf.fa

bedtools getfasta -fi genome.fa -bed genes_hg38_109_both.bed -fo hg38_109_gtf.fa
````