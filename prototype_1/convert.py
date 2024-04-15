def convert(path: str):
    converted_file = ""
    with open(path, "r") as output_file:
        for i, line in enumerate(output_file):
            kmer = line.split(" ")[0]
            converted_file += f">kmer{i}\n{kmer}\n"
    with open("converted_output_stat5b3.fasta", "w") as converted_output:
        converted_output.write(converted_file)


if __name__ == "__main__":
    test_path = "C:/Users/pasca/OneDrive/Desktop/output_stat5b3.txt"
    convert(test_path)
