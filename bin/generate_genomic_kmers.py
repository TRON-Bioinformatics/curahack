#!/usr/bin/env python

from argparse import ArgumentParser

from Bio import SeqIO


def read_fasta(filename):
    """
    Read sequences from a FASTA file.
    Args:
        filename (str): Path to the FASTA file.
    Returns:
        list: List of tuples where each tuple contains the sequence ID and sequence.
    """
    sequences = []
    with open(filename, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            yield (record.id, str(record.seq))
            #sequences.append((record.id, str(record.seq)))
    #return sequences


def generate_kmers(sequence, k, step):
    """Generates non-overlapping k-mers from a sequence."""
    for i in range(0, len(sequence) - k + 1, step):
        yield sequence[i:i+k]


def generate_kmer_fasta(input_filename, output_filename, k, step):
    """Generates a new FASTA file with one entry per k-mer."""
    unique_kmers = {}  # Dictionary to store unique 31-mers and their occurrences
    duplicate_kmers = set()
    
    for identifier, sequence in read_fasta(input_filename):
        for i, kmer in enumerate(generate_kmers(sequence, k, step)):
            if 'N' in kmer:
                continue
            # Are positions correct?
            # Is there a way to find kmer in large sequence and return all positions?
            if kmer not in unique_kmers:
                pos = "{}_{:d}".format(identifier, i*step)
                unique_kmers[kmer] = pos
            else:
                duplicate_kmers.add(kmer)
                    
            if i % 1000000 == 0:
                print(i)
    print(f'Uniques:    {len(unique_kmers)}')
    print(f'Duplicates: {len(duplicate_kmers)}')
    with open(output_filename, 'w') as file:
        for kmer, pos in unique_kmers.items():
            if kmer in duplicate_kmers:
                continue
            file.write(f'>{pos}\n')
            file.write(f'{kmer}\n')



def main():
    parser = ArgumentParser(description="Generates kmers from genomic FASTA")
    parser.add_argument("-i", "--input_fasta", dest="input_fasta", help="Specify input FASTA file")
    parser.add_argument("-s", "--step_size", dest="step_size", type=int, help="Specify step size during kmer generation", default=15)
    parser.add_argument("-k", "--kmer_size", dest="kmer_size", type=int, help="Specify kmer size", default=31)
    parser.add_argument("-o", "--output_fasta", dest="output_fasta", help="Specify output FASTA file")

    args = parser.parse_args()

    k = args.kmer_size
    generate_kmer_fasta(args.input_fasta, args.output_fasta, k, args.step_size)
    print(f'Generated new FASTA file with {k}-mers: {args.output_fasta}')


if __name__ == '__main__':
    main()