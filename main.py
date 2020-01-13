import argparse

def parse_fasta(fname):
    with open(fname, "r") as fh:

        # Create variables for storing the identifiers and the sequence.
        identifier = None
        sequence = []

        for line in fh:
            line = line.strip()  # Remove trailing newline characters.
            if line.startswith(">"):
                if identifier is None:
                    # This only happens when the first line of the
                    # FASTA file is parsed.
                    identifier = line[1:]
                else:
                    # This happens every time a new FASTA record is
                    # encountered.

                    # Start by yielding the entry that has been built up.
                    yield identifier, sequence

                    # Then reinitialise the identifier and sequence
                    # variables to build up a new record.
                    identifier = line[1:]
                    sequence = []
            else:
                # This happens every time a sequence line is encountered.
                seq_on_line = list(line)
                sequence.extend(seq_on_line)
        # yielding last fasta entry
        yield identifier, sequence

def parse_fastq(fname):
    with open(fname, "r") as fh:

        # Create variables for storing the identifiers and the sequence.
        identifier = None
        sequence = []
        quality = []
        readseq = None

        for line in fh:
            line = line.strip()  # Remove trailing newline characters.
            if line.startswith("@"):
                readseq = True
                if identifier is None:
                    # This only happens when the first line of the
                    # FASTQ file is parsed.
                    identifier = line[1:]
                else:
                    # This happens every time a new FASTQ record is
                    # encountered.

                    # Start by yielding the entry that has been built up.
                    yield identifier, sequence, quality

                    # Then reinitialise the identifier, sequence and quality
                    # variables to build up a new record.
                    identifier = line[1:]
                    sequence = []
                    quality = []
            elif line.startswith("+"):
                readseq = False
            else:
                if readseq:
                    # This happens every time a sequence line is encountered.
                    seq_on_line = list(line)
                    sequence.extend(seq_on_line)
                else:
                    # This happens every time a quality line is encountered.
                    qual_on_line = list(line)
                    quality.extend(qual_on_line)
        yield identifier, sequence, quality

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fastq", type=argparse.FileType('r'),
                    help="input fastq file containing sequenced reads")
    parser.add_argument("--ref_fasta", dest="ref", help="fasta file containing reference sequences")
    args = parser.parse_args()
    return args

def main():
    
    args = argparser()
    print(args.ref)
    


if __name__ == '__main__':
    main()