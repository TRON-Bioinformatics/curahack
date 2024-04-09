#!/usr/bin/env python

import sys

def get_exon_coordinates(gtf_file):
    count = 0
    exon_coordinates = []
    with open(gtf_file) as inf:
        for line in inf:
            if line[0] == "#":
                continue
            elements = line.rstrip().split("\t")
            if elements[1] == "ensembl" and elements[2] == "exon":
                

                chrom = elements[0]
                if chrom != "1":
                    continue
                start = elements[3]
                end = elements[4]
                info_field = elements[8]
                gene_id = ""
                for ele in info_field.split(";"):
                    if ele.strip().startswith("gene_id"):
                        gene_id = ele.split()[1].strip("\"")
                
                count +=1
                exon_coordinates.append((chrom, start, end, gene_id))
    
    print("Number of exons: {}".format(count))
    return exon_coordinates


def write_to_file(exon_coordinates, output_file):
    with open(output_file, "w") as outf:
        for (chrom, start, end, gene_id) in exon_coordinates:
            outf.write("{};{};{};{}\n".format(chrom, start, end, gene_id))




# Example usage:
input_filename = sys.argv[1]
output_filename = sys.argv[2]
exon_coord = get_exon_coordinates(input_filename)
write_to_file(exon_coord, output_filename)
