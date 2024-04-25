#!/usr/bin/env python

from argparse import ArgumentParser

from Bio import SeqIO


def load_genomic_kmers(filename):
    """
    Read sequences from a FASTA file.
    Args:
        filename (str): Path to the FASTA file.
    Returns:
        list: List of tuples where each tuple contains the sequence ID and sequence.
    """
    kmer_dict = {}
    kmer_size = 0
    with open(filename, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            kmer_dict[str(record.seq)] = record.id
            kmer_size = len(str(record.seq))
    return kmer_dict, kmer_size


def read_fastq(filename):
    """
    Read sequences from a FASTQ file.
    Args:
        filename (str): Path to the FASTQ file.
    Returns:
        list: List of tuples where each tuple contains the sequence ID and sequence.
    """
    sequences = []
    with open(filename, "r") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            yield (record.id, str(record.seq))
            #sequences.append((record.id, str(record.seq)))
    #return sequences


def generate_kmers(sequence, k):
    """Generates non-overlapping 31-mers from a sequence."""
    for i in range(0, len(sequence) - k + 1):
        yield sequence[i:i+k]


def compare_sequences_with_mismatch(seq1, seq2):
    """
    Compare two sequences allowing for one mismatch.
    
    Args:
        seq1 (str): First sequence.
        seq2 (str): Second sequence.
        
    Returns:
        bool: True if sequences match with at most one mismatch, False otherwise.
    """
    if len(seq1) != len(seq2):
        return False
    
    mismatches = 0
    for base1, base2 in zip(seq1, seq2):
        if base1 != base2:
            mismatches += 1
            if mismatches > 1:
                return False
    
    return True


def count_kmers_in_sequences(fastq_file, kmer_dict, k):
    """Counts occurrences of each entry in the dictionary within all sequences."""
    n = 0
    counts = {}
    counts_norm = {}
    for identifier, sequence in read_fastq(fastq_file):
        n += 1
        if n % 100000 == 0:
            print("Processed {} reads".format(n))
        count = 0
        kmers = []
        for query_kmer in generate_kmers(sequence, k):
            # Exact matching of kmers
            # TODO: Check if allowing mismatches would improve results
            if query_kmer in kmer_dict:
            #for ref_kmer in kmer_dict:
                #if compare_sequences_with_mismatch(query_kmer, ref_kmer):
                count += 1
                kmers.append(kmer_dict[query_kmer])
        for kmer in kmers:
            # Only count each read once: if multiple kmers are present in one read, just add the fraction to the count table
            try:
                counts_norm[kmer] += 1/count
            except:
                counts_norm[kmer] = 1/count
            try:
                counts[kmer] += 1
            except:
                counts[kmer] = 1
    return counts, counts_norm


def write_read_count_fasta(output_filename, counts):
    """Writes counts to a FASTA file."""
    print(f"Number of exported counts: {len(counts)}")
    with open(output_filename, 'w') as file:
        for identifier, count in counts.items():
            file.write(f'>{identifier}\n')
            file.write(f'{count}\n')


def write_read_count_csv(output_filename, counts):
    """Writes counts to a csv file."""
    print(f"Number of exported counts: {len(counts)}")
    with open(output_filename, 'w') as file:
        for identifier, count in counts.items():
            file.write(f'{identifier};{count}\n')


def main():
    parser = ArgumentParser(description="Counts kmers from FASTQ reads in genomic kmer DB")
    parser.add_argument("-i", "--input_fastq", dest="input_fastq", help="Specify input FASTQ file")
    parser.add_argument("-g", "--genomic_kmers", dest="genomic_kmers", help="Specify genomic kmers FASTA")
    #parser.add_argument("-k", "--kmer_size", dest="kmer_size", type=int, help="Specify kmer size", default=31)
    parser.add_argument("-o", "--output_csv", dest="output_csv", help="Specify output CSV file")

    args = parser.parse_args()

    kmer_dict, k = load_genomic_kmers(args.genomic_kmers)
    print("Loaded genomic kmers.")
    counts, counts_norm = count_kmers_in_sequences(args.input_fastq, kmer_dict, k)
    write_read_count_csv(args.output_csv, counts)
    write_read_count_csv(args.output_csv + ".norm.csv", counts_norm)


if __name__ == '__main__':
    main()
    