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
    with open(filename, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            kmer_dict[str(record.seq)] = record.id
    return kmer_dict


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


def count_kmers_in_sequences(fastq_file, kmer_dict, k):
    """Counts occurrences of each entry in the dictionary within all sequences."""
    n = 0
    #num_sequences = len(sequences)
    counts = {} #{kmer: 0 for kmer in kmer_dictionary}
    for identifier, sequence in read_fastq(fastq_file):
        n += 1
        if n % 100000 == 0:
            print("Processed {} reads".format(n))
        for kmer in generate_kmers(sequence, k):
            if kmer in kmer_dict:
                try:
                    counts[kmer_dict[kmer]] += 1
                except:
                    counts[kmer_dict[kmer]] = 1
    return counts


def count_kmers_in_sequences_with_fractions(sequences, kmer_dict, k):
    """Counts occurrences of each entry in the dictionary within all sequences."""
    counts = {}
    for identifier, sequence in sequences:
        count = 0
        kmers = []
        for kmer in generate_kmers(sequence, k):
            if kmer in kmer_dict:
                count += 1
                kmers.append(kmer_dict[kmer])
        for kmer in kmers:
            # Only count each read once: if multiple kmers are present in one read, just add the fraction to the count table
            try:
                counts[kmer] += 1/count
            except:
                counts[kmer] = 1/count
    return counts


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
    parser.add_argument("-k", "--kmer_size", dest="kmer_size", type=int, help="Specify kmer size", default=31)
    parser.add_argument("-o", "--output_csv", dest="output_csv", help="Specify output CSV file")

    args = parser.parse_args()

    kmer_dict = load_genomic_kmers(args.genomic_kmers)
    print("Loaded genomic kmers.")
    counts = count_kmers_in_sequences(args.input_fastq, kmer_dict, args.kmer_size)
    write_read_count_csv(args.output_csv, counts)


if __name__ == '__main__':
    main()
    