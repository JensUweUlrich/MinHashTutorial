import argparse
import matplotlib.pyplot as plt

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
    parser.add_argument("--fastq", dest="fastq",
                    help="input fastq file containing sequenced reads")
    parser.add_argument("--ref_fasta", dest="ref", help="fasta file containing reference sequences")
    args = parser.parse_args()
    return args

def calculate_read_metrics(fastq_file):
    # initialize counting variables
    read_number = 0
    bp_number = 0
    min_readlength = None
    max_readlength = 0
    read_n50 = None
    read_lengths = []
    mean_qual = None
    sum_qual = 0
    # iterate over all reads in the given fastq file
    for identifier, sequence, quality in parse_fastq(fastq_file):
        # increase read counter by 1 for every read
        read_number += 1
        # add bp count of current read to sum of all read base pairs
        bp_number += len(sequence)
        # set length of first read to minimum read length
        # or update minimum read length if current read is smaller than minimum read length
        if min_readlength == None:
            min_readlength = len(sequence)
        elif len(sequence) < min_readlength:
            min_readlength = len(sequence)

        # update maximum read length if current read length is bigger than maximum read length
        if len(sequence) > max_readlength:
            max_readlength = len(sequence)
        
        read_lengths.append(len(sequence))
        # increase counter for current read length in read length dictionary
        #if len(sequence) in read_lengths:
        #    read_lengths[len(sequence)] += 1
        # else:
        #    read_lengths[len(sequence)] = 1

        
        # first sum up all quality values of all bases of the read
        for q in quality:#
            # use Sanger encoding of quality values
            sum_qual += ord(q) - 33
        # add mean quality of the read to the sum of mean read qualities

    # get mean quality score of all reads
    mean_qual = sum_qual / bp_number
    
    # set counter to 0 for not occuring read lengths
    #for i in range(1, max_readlength, 1):
        #if not i in read_lengths:
            #read_lengths[i] = 0

    n50_bp = 0;
    # iterate read length counter dict in descending order
    for i in sorted(read_lengths, reverse=True):
        # add bp in reads with that length to the counter
        n50_bp += i
        # if more than half of all bp is in reads longer than i => stop iteration
        # and set read_n50 to i
        if n50_bp >= bp_number/2:
            read_n50 = i
            break


    print("Number Reads : ", read_number)
    print("Number of Base Pairs : ", bp_number)
    print("Maximum Read Length : ", max_readlength)
    print("Minimum Read Length : ", min_readlength)
    print("Mean Read Length : ", bp_number/read_number)
    print("Read N50 : ", read_n50)
    print("Mean Quality Score : ", mean_qual)
    return read_lengths

def main():
    
    args = argparser()
    lengths = calculate_read_metrics(args.fastq)

    #print(lengths.values())
    plt.hist(lengths, bins=50)
    plt.show()
    
    


if __name__ == '__main__':
    main()