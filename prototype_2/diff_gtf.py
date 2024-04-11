from typing import Tuple, Dict

import gtfparse
import polars as pl
import argparse


# open and load a gtf file with gtfparse
def load_gtf(file_path):
    return gtfparse.read_gtf(file_path)


# for every gene in the gtf file, it creates a unique key  with chr-start-end and calculates a hash over all exons
# coordinates
def hash_genes_in_gtf(gtf) -> Dict[str, Tuple[str, str]]:
    gene_dict = {}
    # iterate over genes
    for gene in gtf.filter(pl.col("feature") == "gene").iter_rows(named=True):
        gene_id = gene['gene_id']

        # define a gene unique key
        # NOTE: it does not use the gene id in case they change between versions
        chromosome = gene['seqname']
        start = gene['start']
        end = gene['end']
        strand = gene['strand']
        key = chromosome + '_' + str(start) + '_' + str(end) + '_' + strand

        exons_coordinates = []

        # iterate over exons for this gene_id in genomic order
        exons_by_gene_id = gtf.filter((pl.col("gene_id") == gene_id) & (pl.col("feature") == "exon")) \
            .sort(["start", "end"])
        for exon in exons_by_gene_id.iter_rows(named=True):
            exons_coordinates.append(str(exon['start']) + "_" + str(exon['end']))

        # calculate a hash over all exons coordinates
        exons_hash = hash(','.join(exons_coordinates))

        # store the gene_id and the hash in a dictionary
        gene_dict[key] = (gene_id, exons_hash)

    return gene_dict


def compare_two_gtfs(gene_dict_1, gene_dict_2):
    # find genes that are in gene_dict_1 but not in gene_dict_2
    diff_1 = diff_gene_dicts(gene_dict_1, gene_dict_2)

    # find genes that are in gene_dict_2 but not in gene_dict_1
    diff_2 = diff_gene_dicts(gene_dict_2, gene_dict_1)

    common_gene_ids_modified_coords = set(diff_1.values()).intersection(set(diff_2.values()))
    print("There are {} shared gene ids with modified coordinates in both GTFs".format(
        len(common_gene_ids_modified_coords)))
    print(common_gene_ids_modified_coords)

    unique_gene_ids_1 = set(diff_1.values()).difference(set(diff_2.values()))
    print("There are {} unique gene ids in GTF 1".format(len(unique_gene_ids_1)))
    print(unique_gene_ids_1)

    unique_gene_ids_2 = set(diff_2.values()).difference(set(diff_1.values()))
    print("There are {} unique gene ids in GTF 2".format(len(unique_gene_ids_2)))
    print(unique_gene_ids_2)

    # find genes that are in both gene_dict_1 and gene_dict_2
    common = {k: (i, h) for k, (i, h) in gene_dict_1.items() if k in gene_dict_2}
    print("There are {} common genes with equal gene coordinates in both GTFs".format(len(common)))

    identical = [(i, gene_dict_2[k][0]) for k, (i, h) in gene_dict_1.items() if k in gene_dict_2 and h == gene_dict_2[k][1]]
    print("There are {} common genes with identical coordinates across all exons in both GTFs".format(len(identical)))

    identical = pl.DataFrame({"gene_id_1": [i for i, i2 in identical], "gene_id_2": [i2 for i, i2 in identical]})

    jaccard_index = len(identical) / (len(gene_dict_1) + len(gene_dict_2) - len(identical))
    print("Jaccard index: {}".format(jaccard_index))

    return identical


def diff_gene_dicts(gene_dict_1, gene_dict_2):
    return {k: i for k, (i, _) in gene_dict_1.items() if k not in gene_dict_2}


def main():
    # create a command line interface with argparse
    parser = argparse.ArgumentParser(description='Diff two GTF files')
    parser.add_argument('gtf_file', type=str, help='GTF file 1')
    parser.add_argument('gtf_file_2', type=str, help='GTF file 2')
    parser.add_argument('--output', type=str, help='Output file', default='identical_genes.csv')
    args = parser.parse_args()

    # load the GTF files
    gtf = load_gtf(args.gtf_file)
    gtf_2 = load_gtf(args.gtf_file_2)

    # hash the genes in the GTF files
    gene_dict = hash_genes_in_gtf(gtf)
    gene_dict_2 = hash_genes_in_gtf(gtf_2)

    # compare the GTF files
    identical_df = compare_two_gtfs(gene_dict, gene_dict_2)

    # save the result to a file
    identical_df.write_csv(args.output)


if __name__ == "__main__":
    main()


