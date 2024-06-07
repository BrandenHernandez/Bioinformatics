"""============================================================================
Title : FASTAsorter.py

Description : Python module for reading and writing FASTA files
which will be used in most of the assignments in this course.

Author: Branden Hernandez

Date: 9-18-23

Version: 0.1

Usage: Read through the entire FASTA file taken as input (-i)
then sort the sequences from the longest sequence down to the shortest sequence.

Python Version: 3.9
============================================================================= """

import sys

# About:
# FASTA File Specifications
# The FASTA file format is a widely used bioinformatics text format for representing nucleotide sequences
# (such as DNA or RNA) or protein sequences. The file is entirely text based using standard ASCII characters.
# The format is as follows:
# 
# 1. Header Line:
#   a. Each sequence in the FASTA format begins with a single-line description, also known as the
# sequence identifier. This line is distinguished from the sequence data by a greater-than (">")
# symbol at the beginning.
# 
#   b. It is recommended, but not required, to keep the first word of the description line as the unique
# identifier for the sequence.
# 
#   c. Example:
# >sequence_identifier Description of the sequence
# 
# 2. Sequence data:
#   a. Following the header line, the lines of sequence data start.
# 
#   b. These lines contain the actual nucleotide (A, T, C, G, U for DNA/RNA) or protein (represented by
# standard one-letter codes) sequences.
# 
#   c. There should be no whitespace characters within the sequence. The sequence can span multiple
# lines, but no line should be longer than 80 characters, which is a rule to make the file easily
# readable.
# 
#   d. Example:
# AGCTTTTCATTCTGACTGCAACGGGCAATATGTCT
# CTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC
# 
# 3. The header line and the sequence data are tied together as a single unit, where the header line acts as the
# metadata for the sequence and the sequence data is the data.
# 
# 4. Multiple sequences:
#   a. In the case of multiple sequences, each sequence must follow the two-step format of a single
# header line followed by sequence lines.
# 
#   b. The header line of the next sequence marks the end of the previous sequence.
# 
#   c. Example:
# >sequence1 AGCTTTTCATTCTGACTGCAACGGGCAATATGTCT
#            CTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC
# >sequence2 CTAGCTAGCTAGCTGACGATGCGATTACGTATCGT
#            ATCGGCTAGGCTAGCTAGCTCGATCGAGCTAGCTA

# Building Blocks
gene = []

blah = sys.argv


def i(input_path):
    input_file_path = input_path
    # input_file_path = "100_sequences.fna"  # test
    return input_file_path


def o(output_path):
    output_file_path = output_path
    return output_file_path


def read_from_file(i_file_path):  # accept file path as argument
    with open(i_file_path, 'r') as file:  # open file path in read mode
        data = file.read()  # stuff data into variable 'data'
    return data  # return variable


def write_to_file(o_file_path, data):  # take file path and data as arguments
    with open(o_file_path, 'w') as file:  # open the file path in write mode
        file.write(data)  # write data contents to file at file path


def reorganize(data):
    for num in range(1, len(data)):
        key_item = data[num]
        j = num - 1
        while j >= 0 and len(data[j]) < len(key_item):
            temp = data[j]
            data[j] = data[j + 1]
            data[j + 1] = temp
            j -= 1
    return data


def parser3000(data):
    parsed = []
    parsed = data.split(">")
    ordered = reorganize(parsed)
    return ordered


# print("Should be file path " + blah[0])  # test
# print("Should be i " + blah[1])  # test
# print("Should be input " + blah[2])  # test
# print("Should be o " + blah[3])  # test
# print("Should be output " + blah[4])  # test

if blah[2].__contains__('.fna'):  # extension check
    buffer = ""
    # take input file, feed to reader function, feed to parser function, then build list
    gene = parser3000(read_from_file(i(blah[2])))

    # go through each item in list
    for i in range(1, len(gene)):
        # print each gene sequence
        buffer += (">" + gene[i] + "\n")

    write_to_file(o(blah[4]), buffer)

else:
    sys.exit("\nERROR: Please upload a .fna file path\n")

print('\n' + read_from_file(o(blah[4])))

# Thanks for looking
