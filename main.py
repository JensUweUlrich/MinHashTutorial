import argparse
import random
import bisect
import matplotlib.pyplot as plt

nts = {}
rand = random.randint(1,1000)

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
    parser.add_argument("--kmer-size", dest="K", type=int, help="size of the kmers used for minhashing")
    parser.add_argument("--sketch-size", dest="S", type=int, help="size of the minhash sketches")
    parser.add_argument("--output", dest="out", help="fastq output file")
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
    read_qualities = []
    mean_qual = None
    sum_qual = 0
    # iterate over all reads in the given fastq file
    for identifier, sequence, quality in parse_fastq(fastq_file):
        # increase read counter by 1 for every read
        read_number += 1
        # add bp count of current read to sum of all read base pairs
        bp_number += len(sequence)
        # set length of first read to minimum read length
        # or update minimum read length if current read is smaller than minimum
        # read length
        if min_readlength == None:
            min_readlength = len(sequence)
        elif len(sequence) < min_readlength:
            min_readlength = len(sequence)

        # update maximum read length if current read length is bigger than
        # maximum read length
        if len(sequence) > max_readlength:
            max_readlength = len(sequence)
        
        read_lengths.append(len(sequence))
        
        current_quality = 0
        # first sum up all quality values of all bases of the read
        for q in quality:#
            # use Sanger encoding of quality values
            sum_qual += ord(q) - 33
            current_quality += ord(q) - 33

        if len(quality) > 0:
            read_qualities.append(current_quality / len(quality))
    # get mean quality score of all reads
    mean_qual = sum_qual / bp_number
    
    n50_bp = 0
    # iterate read length counter dict in descending order
    for i in sorted(read_lengths, reverse=True):
        # add bp in reads with that length to the counter
        n50_bp += i
        # if more than half of all bp is in reads longer than i => stop
        # iteration
        # and set read_n50 to i
        if n50_bp >= bp_number / 2:
            read_n50 = i
            break


    print("Number Reads : ", read_number)
    print("Number of Base Pairs : ", bp_number)
    print("Maximum Read Length : ", max_readlength)
    print("Minimum Read Length : ", min_readlength)
    print("Mean Read Length : ", bp_number / read_number)
    print("Read N50 : ", read_n50)
    print("Mean Quality Score : ", mean_qual)

    plt.hist(read_lengths, bins=500)
    plt.show()

    plt.hist(read_qualities)
    plt.show()
    
def initialize_nts():
    nts['A'] = 0
    nts['C'] = 1
    nts['G'] = 2
    nts['T'] = 3
    nts['N'] = 4

# calculates the hash value for a given kmer
def hash_kmer(kmer):
    # reverse the kmer first
    kmer_reverse = kmer[::-1]
    hash_v = 0
    # iterate over all characters of the kmer for hash value computation
    for i in range(0, len(kmer)):
        hash_v += int(nts[kmer_reverse[i]]) * (4 ** i)

    # XOR hash value with a random integer in order to avoid lexicographical
    # sorting of kmers
    # otherwise kmers beginning with A would always have the smallest hash
    # value
    # could interfere jaccard similarity for small sketch sizes
    hash_v ^= rand
    return hash_v

# calculate reverse complement of a given kmer
def reverse_complement(kmer):
    kmer_reverse = kmer[::-1]
    rev_comp = []
    for n in kmer_reverse:
        if n == 'A':
            rev_comp.append('T')
            continue
        if n == 'C':
            rev_comp.append('G')
            continue
        if n == 'G':
            rev_comp.append('C')
            continue
        if n == 'T':
            rev_comp.append('A')
            continue
        if n == 'N':
            rev_comp.append('N')

    return rev_comp


# builds the minhash sketch for a give dna sequence
# sequence: dna sequence
# k : kmer size
# s : sketch size
def build_sketch(sequence, k, s):
    # initialize resulting sketch of hash values
    sketch = []
    # iterate over all kmers of the given sequence
    # store hash values of the kmers within a list
    max_h = 0
    for i in range(0, len(sequence) - k):
        current_kmer = sequence[i:i + k]
        h = hash_kmer(current_kmer)
        # take only thge smaller hash value of current kmer abnd its reverse complement
        h_rc = hash_kmer(reverse_complement(current_kmer))
        if h_rc < h:
            h = h_rc

        # add first hash value to the sketch 
        if max_h == 0:
            sketch.append(h)
            max_h = h
            continue

        # add further hash values to the sketch if h is not already in the sketch
        # and either the sketch is not full or hash value is smaller than the maximum hash value in the sketch
        if not h in sketch and (len(sketch) < s or h < max_h):
            # insert hash value and keep list sorted
            bisect.insort(sketch, h)
            if h > max_h:
                max_h = h
            if len(sketch) > s:
                sketch = sketch[:s]

    return sketch

# calculate jaccard similarity for given sketches
def estimate_jaccard(sketch1, sketch2):
    # initialize counter for hash values contained in both sketches
    counter = 0
    # if sketches are of different size, iterate over the elements of the
    # smaller sketch
    # increment the counter if hash value is part of both sketches
    s = len(sketch1)
    if (len(sketch1) <= len(sketch2)):
        for h in sketch1:
            if h in sketch2:
                counter += 1
    else:
        s = len(sketch2)
        for h in sketch2:
            if h in sketch1:
                counter += 1
    if s == 0:
        return 0
    else:
        return counter / s


def write_fastq_record(fh, id, seq, qual):
    # sequence identifier line
    fh.write("@" + id + "\n")
    # divide sequence in chunks of 70 chars and print lines seperately
    for i in range(0, int(len(seq) / 70) + 1):
        s = seq[(i * 70) : ((i + 1) * 70)]
        fh.write("".join(s) + "\n")
    
    # quality identifier line
    fh.write("+" + id + "\n")
    # divide qualities in chunks of 70 chars and print lines seperately
    for i in range(0, int(len(qual) / 70) + 1):
        s = qual[(i * 70) : ((i + 1) * 70)]
        fh.write("".join(s) + "\n")
    
def add_entry(dict, identifier, sequence, quality, jaccard):
    entry = {}
    entry['jaccard'] = jaccard
    entry['sequence'] = sequence
    entry['quality'] = quality
    dict[identifier] = entry
   
def main():
    
    args = argparser()
    # initialize nt translation
    initialize_nts()

    # build sketch of positive control
    ref_sketches = {}
    for id, seq in parse_fasta(args.ref):
        ref_sketches[id] = build_sketch(seq, args.K, args.S)
    
    # iterate over all reads in fastq file
    
    result = {}
    min_jacc = 1.0
    for identifier, sequence, quality in parse_fastq(args.fastq):
        # build sketch for current read
        read_sketch = build_sketch(sequence, args.K, args.S)
        for ref_sketch in ref_sketches.values():
            # calculate jaccard similarity for reference and read sketch
            jacc = estimate_jaccard(ref_sketch, read_sketch)
            # add first 10 entries to the dictionary
            if len(result) < 10:
                add_entry(result, identifier, sequence, quality, jacc)
                if float(jacc) < float(min_jacc):
                    min_jacc = jacc
            else:
                if float(jacc) > float(min_jacc):
                    d_key = None
                    for k, v in result.items():
                        if float(v['jaccard']) == float(min_jacc):
                            d_key = k
                            break
                    del result[d_key]
                    add_entry(result, identifier, sequence, quality, jacc)
                    min_jacc = 1.0
                    for k, v in result.items():
                        if float(v['jaccard']) < float(min_jacc):
                            min_jacc = v['jaccard']

    fh = open(args.output, "w")
    for k,v in result.items():
        write_fastq_record(fh, k, v['sequence'], v['quality'])
    
    fh.close()
    
    
    

if __name__ == '__main__':
    main()