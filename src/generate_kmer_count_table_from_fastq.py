
def read_kmer_fasta(filename):
    """Reads a FASTA file and returns a list of tuples containing sequence IDs and sequences."""
    sequences = {}
    with open(filename, 'r') as file:
        identifier = None
        sequence = ''
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if identifier is not None:
                    #sequences.append((identifier, sequence))
                    sequences[sequence] = identifier
                identifier = line[1:]
                sequence = ''
            else:
                sequence += line
        if identifier is not None:
            sequences[sequence] = identifier
            #sequences.append((identifier, sequence))
    return sequences

def read_fastq(filename):
    """Reads a FastQ file and returns a list of tuples containing sequence IDs and sequences."""
    sequences = {}
    with open(filename, 'r') as file:
        identifier = None
        sequence = ''
        for line in file:
            line = line.strip()
            if line.startswith('@'):  # FastQ header
                if identifier is not None:
                    #sequences.append((identifier, sequence))
                    sequences[identifier] = sequence
                identifier = line[1:]
                sequence = ''
            elif line.startswith('+'):  # Quality header, skip
                continue
            else:
                sequence += line
        if identifier is not None:
            #sequences.append((identifier, sequence))
            sequences[identifier] = sequence
    return sequences

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

def generate_kmers(sequence, k):
    """Generates non-overlapping 31-mers from a sequence."""
    for i in range(0, len(sequence) - k + 1):
        yield sequence[i:i+k]

def count_kmers_in_sequences(sequences, kmer_dictionary, k):
    """Counts occurrences of each entry in the dictionary within all sequences."""
    counts = {}#{kmer: 0 for kmer in kmer_dictionary}
    for sequence in sequences.values():
        for kmer in generate_kmers(sequence, k):
            if kmer in kmer_dictionary:
                try:
                    counts[kmer_dictionary[kmer]] += 1
                except:
                    counts[kmer_dictionary[kmer]] = 1
    return counts

def count_kmers_in_sequences_with_fractions(sequences, kmer_dictionary, k):
    """Counts occurrences of each entry in the dictionary within all sequences."""
    counts = {}
    for sequence in sequences.values():
        count = 0
        kmers = []
        for kmer in generate_kmers(sequence, k):
            if kmer in kmer_dictionary:
                count = count + 1
                kmers.append(kmer_dictionary[kmer])
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

if __name__ == '__main__':
    kmer_file = 'kmers_chr1_k31_step15.fasta'
    reads_file = 'Control_siRNA_1.fastq'
    k = 31
    output_filename = 'kmer_read_counts_control_siRNA_1_k{k}.csv'
    kmer_dictionary = read_kmer_fasta(kmer_file)
    sequences = read_fastq(reads_file)
    counts = count_kmers_in_sequences(sequences, kmer_dictionary, k)
    write_read_count_csv(output_filename, counts)