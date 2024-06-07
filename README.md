FASTA File Specifications
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
