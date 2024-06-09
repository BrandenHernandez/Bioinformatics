import sys
import numpy

"""============================================================================
Title : Smith-Waterman_matrix.py

Description : Python module for reading FASTA file and 
using Smith-Waterman matrix to determine best fit for localized sequence alignment.

Author: Branden Hernandez

Date: 10-13-23

Version: 0.3

Usage: Read through a FASTA file containing two sequences taken as input (-i) then align the first
sequence to the second sequence using the local alignment algorithm Smith-Waterman. The resulting aligned
subsequences will then be written to an output file as given by the output argument (-o)

Organization: Texas Tech University

Python Version: 3.9
============================================================================= """

# Building Blocks
gene = []
blah = sys.argv
DIAG = 0
UP = 1
LEFT = 2
X = -1
arrows = [u"\u2196", u"\u2191", u"\u2190"]  # diagonal, up, left
output_file = ""
global sub_file


# From Assignment I
# --------------------------------------------------------------------------

def read_from_file(i_file_path):  # accept file path as argument
    with open(i_file_path, 'r') as file:  # open file path in read mode
        file_data = file.read()  # stuff data into variable 'data'
    return file_data  # return variable


def write_to_file(o_file_path, input_data):  # take file path and data as arguments
    with open(o_file_path, 'w') as file:  # open the file path in write mode
        file.write(input_data)  # write data contents to file at file path


def reorganize(input_data):
    for num in range(1, len(input_data)):
        key_item = input_data[num]
        j = num - 1
        while j >= 0 and len(input_data[j]) < len(key_item):
            temp = input_data[j]
            input_data[j] = input_data[j + 1]
            input_data[j + 1] = temp
            j -= 1
    return input_data


def parser3000(input_data):
    parsed = []
    parsed = input_data.split(">")
    ordered = reorganize(parsed)
    return ordered


def reader(genes, user_input, output):
    # print("Should be file path " + user_input[0])  # test
    # print("Should be -i " + user_input[1])  # test
    # print("Should be input file " + user_input[2])  # test
    # print("Should be -o " + user_input[3])  # test
    # print("Should be output file " + user_input[4])  # test
    # print("Should be -s " + user_input[5])  # test
    # print("Should be substitution file " + user_input[6])  # test
    # print(" test" + read_from_file(blah[6])) # test


    if user_input[1].__contains__('i') and user_input[2].__contains__('.fna'):  # extension check
        buffer = ""  # create a string to combine output
        # take input file, feed to reader function, feed to parser function, then build list
        genes = parser3000(read_from_file(user_input[2]))

        if user_input[3].__contains__('o'):
            output = user_input[4]

        if user_input[5].__contains__('-s'):
            sub_file = read_from_file(user_input[6])

        # go through each item in list
        for i in range(len(genes)):
            # print each gene sequence
            buffer += (">" + genes[i] + "\n")

    else:
        sys.exit("\nERROR: Please upload a .fna file path\n")

    # print('\n' + read_from_file(output_file))  # test
    # print('\n ' + read_from_file(sub_file))  # test
    # return 2 sequences to compare
    return output, buffer, sub_file


# From Assignment II
# --------------------------------------------------------------------------
def transmute_negatives(num):
    if num < 0:
        num = 0
    return num


def matrix_builder(seq1, seq2, indel):
    # Building Blocks
    match = 1
    miss = -1
    weight = 0
    n = len(seq1)
    m = len(seq2)
    s = numpy.zeros((n + 1, m + 1))  # matrix
    ptr = numpy.zeros((n + 1, m + 1), dtype=int)

    # initialize outer values
    for i in range(1, n + 1):  # side row
        s[i, 0] = 0 * i
    for j in range(1, m + 1):  # top row
        s[0, j] = 0 * j

    # initialize traceback outer arrows
    ptr[0, 0] = X
    ptr[0, 1:] = LEFT
    ptr[1:, 0] = UP

    # Grader =======================================
    for i in range(1, n + 1):  # row
        for j in range(1, m + 1):  # column
            # match
            if seq1[i - 1] == seq2[j - 1]:
                s[i, j] = s[i - 1, j - 1] + match  # enter value into matrix
                ptr[i, j] = DIAG  # enter value into pointer for traceback
                weight = s[i, j]
            # mismatch
            else:
                weight = miss + s[i - 1, j - 1]
                s[i, j] = transmute_negatives(s[i - 1, j - 1] + miss)  # enter value into matrix
                ptr[i, j] = DIAG  # enter value into pointer for traceback
                if s[i, j] == 0.0:
                    ptr[i, j] = X
            # compare gap penalty row up
            if s[i - 1, j] + indel > weight or s[i - 1, j] + indel > s[i, j]:
                s[i, j] = transmute_negatives(s[i - 1, j] + indel)  # enter value into matrix
                ptr[i, j] = UP  # enter value into pointer for traceback
                if s[i, j] == 0.0:
                    ptr[i, j] = X
            # compare gap penalty column up
            if s[i, j - 1] + indel > weight or s[i, j - 1] + indel > s[i, j]:
                s[i, j] = transmute_negatives(s[i, j - 1] + indel)  # enter value into matrix
                ptr[i, j] = LEFT  # enter value into pointer for traceback

    return s, ptr


# start with the highest value and trace until a 0 is met

def matrix_trace(seq1, seq2, ptr, i, j):
    align1 = ""
    align2 = ""
    n, m = (i, j)
    curr = ptr[i, j]
    while n > 0 or m > 0:
        if curr == DIAG:  # print both
            align1 = seq1[n - 1] + align1
            align2 = seq2[m - 1] + align2
            n -= 1
            m -= 1
        elif curr == LEFT:
            align1 = "-" + align1  # gap one
            align2 = seq2[m - 1] + align2  # print the other
            m -= 1
        elif curr == UP:
            align1 = seq1[n - 1] + align1  # print the one
            align2 = "-" + align2  # gap the other
            n -= 1
        elif curr == X:
            return align1, align2

        curr = ptr[n, m]

    return align1, align2


def print_nw_matrix(s, seq1, seq2):
    max_val = 0, 0, 0
    print("\n" + "~~~" * 25)
    print("\tSmith-Waterman matrix:\n")
    print("\t    " + " \t".join(seq2))  # Distance from corner + space between letters. join(letters)
    for i in range(len(s)):
        if i == 0:  # initial gap
            print("  ", end=" ")
        if i > 0:  # row letters
            print(" " + seq1[i - 1], end=" ")
        for j in range(len(s[i])):
            if s[i, j] >= max_val[0]:
                max_val = s[i, j], i, j  # highest value and position
            # space before the number to account for the lack of negative
            print("{val}".format(val=s[i, j]), end=" ")
        if i >= 0:  # new line
            print("")
    return max_val


def print_ptr_matrix(ptr, seq1, seq2):
    print("\n" + "~~~" * 25)
    print("\tTraceback:\n")
    print("\t  " + "  ".join(seq2))
    for i in range(len(ptr)):
        if i == 0:
            print(" ", end=" ")
        if i > 0:
            print(seq1[i - 1], end=" ")  # print letters
        for j in range(len(ptr[i])):
            if ptr[i, j] == -1:
                print(" " + "x", end=" ")
            elif ptr[i, j] >= 3:
                print(" " + arrows[ptr[i, j] - 3], end=" ")
            else:
                # if ptr[i, j] <= 3:
                print(" " + arrows[ptr[i, j]], end=" ")

        if i >= 0:
            print("  ")


def smith_waterman(seq1, seq2, gapped):
    # Build the value matrix
    s, ptr = matrix_builder(seq1, seq2, gapped)
    # Print the matrix and get best score
    score, i, j = print_nw_matrix(s, seq1, seq2)
    # Build the trace back matrix and get the alignment
    alignment = matrix_trace(seq1, seq2, ptr, i, j)
    # Print the pointer matrix
    print_ptr_matrix(ptr, seq1, seq2)

    return alignment, score


# --------------------------------------------------------------------------

if __name__ == '__main__':
    blah = ['Smith-Waterman_matrix.py', '-i', 'local_example.fna',
            '-o', 'test.fna', '-s', 'nucleotide.mtx']  # test

    # read_and_write
    data = reader(gene, blah, output_file)
    #output = data[0]
    #buffer = data[1]
    #sub = data[2]
    # print("test " + data[1])

    sequences = data[1].splitlines()
    sequencesA = sequences[0]  # first input
    sequence1 = sequences[1]  # first sequence
    sequencesB = sequences[3]  # second input
    sequence2 = sequences[4]  # second sequence

    # print("sequence 1: ", sequence1)  # test the first sequence
    # print("sequence 2: ", sequence2)  # test the second sequence

    if data[2].__contains__("-"):
        gap = int(data[2][-5][-1] + data[2][-4][-1])

    # feed both sequence into the N_W algo and spit out the aligned seq + the score
    align, sequences = smith_waterman(sequence1, sequence2, gap)

    # print aligned sequences with score
    print("\n" + "~~~" * 25)
    print("\tAlignment:\n")
    # first string
    print(sequencesA + "; score = {score}\n {align}".format(score=sequences, align=align[0]))
    # second string gapped
    print(sequencesB + "; score = {score}\n {align}".format(score=sequences, align=align[1]))

    write_to_file(data[0], "{Sequence1}; score = {score1}\n{align1}\n"
                           "{Sequence2}; score = {score2}\n{align2}".format(
        Sequence1=sequencesA, score1=sequences, align1=align[0],
        Sequence2=sequencesB, score2=sequences, align2=align[1]))

# Thanks for looking
