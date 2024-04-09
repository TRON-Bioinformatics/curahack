
import csv
import sys


def overlap_coordinates(genomic_coords, gene_coords):
    overlaps = []
    count = 0
    total = len(genomic_coords)
    for g_coord in genomic_coords:
        count += 1
        if count % 10000 == 0:
            print(f"{count}/{total} kmers processed.")
        for gene_coord in gene_coords:
            if (g_coord[0] == gene_coord[0] and
                g_coord[1] >= gene_coord[1] and
                g_coord[2] <= gene_coord[2]):

                overlaps.append((g_coord, gene_coord, g_coord[3]))
    return overlaps


def parse_exon_coordinates(coord_csv):
    exon_coord = []
    # Open the CSV file
    with open(coord_csv, 'r') as file:
        # Create a CSV reader object
        reader = csv.reader(file, delimiter = ";")
        
        # Iterate over each row in the CSV file
        for row in reader:
            # Process each row
            exon_coord.append((row[0], int(row[1]), int(row[2]), row[3]))

    return exon_coord


def parse_kmer_counts(kmer_counts):
    kmer_coord = []
    # Open the CSV file
    with open(kmer_counts, 'r') as file:
        # Create a CSV reader object
        reader = csv.reader(file, delimiter = ";")
        
        # Iterate over each row in the CSV file
        for row in reader:
            # Process each row
            coord_str = row[0]
            count = float(row[1])
            coord_split = coord_str.split("_")
            for i in range(0, len(coord_split)-1, 2):
                chrom = coord_split[i]
                pos = int(coord_split[i+1])
            
                kmer_coord.append((chrom, pos, pos+31, count))
    return kmer_coord
        

input_counts = sys.argv[1]
input_coordinates = "exons_chr1.csv"

genomic_coordinates = parse_kmer_counts(input_counts)
print("Genomic coordinates parsed.")
gene_coordinates = parse_exon_coordinates(input_coordinates)
print("Exonic coordinates parsed.")


# Example genomic coordinates
# genomic_coordinates = [
#     ('chr1', 1180, 1211),
#     ('chr1', 1500, 1531),
#     ('chr2', 3000, 4000)
# ]

# Example gene coordinates
# gene_coordinates = [
#     ('chr1', 1200, 1800, 'GeneA'),
#     ('chr2', 3200, 3800, 'GeneB')
# ]

# Find overlaps between genomic coordinates and gene coordinates
overlaps = overlap_coordinates(genomic_coordinates, gene_coordinates)

# Process and analyze the overlapping regions
count_dict = {}
for overlap in overlaps:
    
    genomic_region = overlap[0]
    gene_region = overlap[1]
    count = overlap[2]
    #print(f"Overlap found between genomic region {genomic_region} and gene region {gene_region}: kmer_count={count}")
    gene_id = gene_region[3]
    if gene_id in count_dict:
        count_dict[gene_id] += count
    else:
        count_dict[gene_id] = count

with open("results.csv", "w") as outf:
    for key, val in count_dict.items():
        outf.write("{};{}\n".format(key, val))