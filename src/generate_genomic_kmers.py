def read_fasta(filename):
    """Reads a FASTA file and returns a dictionary with sequence identifiers as keys and sequences as values."""
    sequences = {}
    with open(filename, 'r') as file:
        identifier = None
        sequence = ''
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if identifier is not None:
                    sequences[identifier] = sequence
                identifier = line[1:]
                sequence = ''
            else:
                sequence += line
        if identifier is not None:
            sequences[identifier] = sequence
    return sequences

def generate_kmers(sequence, k, step):
    """Generates non-overlapping k-mers from a sequence."""
    for i in range(0, len(sequence) - k + 1, step):
        yield sequence[i:i+k]

def generate_kmer_fasta(input_filename, output_filename, k, step):
    """Generates a new FASTA file with one entry per k-mer."""
    sequences = read_fasta(input_filename)
    unique_kmers = {}  # Dictionary to store unique 31-mers and their occurrences
    duplicates = 0
    uniques = 0
    for identifier, sequence in sequences.items():
        for i, kmer in enumerate(generate_kmers(sequence, k, step)):
            if all(char == 'N' for char in kmer):
                continue
            if kmer not in unique_kmers:
                unique_kmers[kmer] = [(identifier, i*k+1)]
                uniques += 1
            else:
                unique_kmers[kmer].append((identifier, i*k+1))
                duplicates += 1
            if i % 1000000 == 0:
                print(i)
    print(f'Uniques:    {uniques}')
    print(f'Duplicates: {duplicates}')
    new_sequences = {}  # Dictionary to store sequences for output
    for kmer, occurrences in unique_kmers.items():
        identifier = '_'.join([f'{id}_{pos}' for id, pos in occurrences])
        new_sequences[identifier] = kmer

    write_fasta(output_filename, new_sequences)

def write_fasta(output_filename, sequences):
    """Writes sequences to a FASTA file."""
    print(len(sequences))
    with open(output_filename, 'w') as file:
        for identifier, sequence in sequences.items():
            file.write(f'>{identifier}\n')
            file.write(f'{sequence}\n')


if __name__ == '__main__':
    input_filename = 'Homo_sapiens.GRCh38.dna.chromosome.1.fa'
    k = 31
    step = 15
    output_filename = f'kmers_chr1_k{k}_step{step}.fasta'
    generate_kmer_fasta(input_filename, output_filename, k, step)
    print(f'Generated new FASTA file with {k}-mers: {output_filename}')
