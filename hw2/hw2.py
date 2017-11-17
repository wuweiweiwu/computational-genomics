import sys, os
import argparse
import numpy as np
import re
import matplotlib.pyplot as plt

'''
Using argparser to get arguments for query and ref and matches files
'''
def make_arg_parser():
    parser = argparse.ArgumentParser(prog='hw2.py',
                          version="%prog 1.0",
                          formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-q","--query",
                      default=argparse.SUPPRESS,
                      required=True,
                      help="Path to query fasta [required]")
    parser.add_argument("-r","--ref",
                      default=argparse.SUPPRESS,
                      required=True,
                      help="Path to reference fasta [required]")
    parser.add_argument("-m","--matches",
                      default=None,
                      required=False,
                      help="matches.txt file [optional]")
    return parser

'''
basic needleman_wunsch
scoring system:
    -3 for mismatch
    1 for matches
    -2 for gap
    ignore affine gaps / indels
'''
MISMATCH = -3
MATCH = 1
GAP = -2
# regular needleman_wunsch algorithm
def needleman_wunsch(query, ref):
    # create a matrix of (len(ref) + 1) x (len(query) + 1) for the grid
    # it is 1 extra so we account for the top row and leftmost column
    matrix = [[0] * (len(query)+1) for _ in range(len(ref)+ 1)]
    # for loop to iterate throught the matrix
    for i in range(1,len(ref)+ 1):
        for j in range(1,len(query)+ 1):
            # assign a penalty based on if it is a match or not for the diagonal score
            penalty = MISMATCH
            if ref[i-1] == query[j-1]:
                penalty = MATCH
            diag = matrix[i-1][j-1] + penalty
            # assign gap penalties from the up and left
            up = matrix[i-1][j] + GAP
            left = matrix[i][j-1] + GAP
            # assign the current value to the maximum of all 3
            matrix[i][j] = max([diag, up, left])

    # start to reconstuct alignment from bottom right of matrix
    # score is that value at the bottom right
    startingRow, startingColumn = len(ref), len(query)
    score = matrix[startingRow][startingColumn]
    # initially the aligments are empty
    alignment1 = alignment2 = ""

    # indices for iterating through the matrix
    # i = row, j = column
    i, j = startingRow, startingColumn

    # loop until the indices are invalid
    while i > 0 and j > 0:
        # get the value for diagonal, up, and left
        diag = matrix[i-1][j-1]
        up = matrix[i-1][j]
        left = matrix[i][j-1]
        # get the max value. thats the direction we want to go
        maxVal = max([diag, up, left])
        # if we are going diagonal we either have a mismatch or match
        if maxVal == diag:
            # if there is match append the character to each alignment
            if ref[i-1] == query[j-1]:
                alignment1 = ref[i-1] + alignment1
                alignment2 = ref[i-1] + alignment2
            # if there is not a match append a blank "-" to the second alignment
            else:
                alignment1 = ref[i-1] + alignment1
                alignment2 = "-" + alignment2
            # decrement both pointers since we are going diagonal
            i -= 1
            j -= 1
        # if we are going left that means that there is a blank in the first alignment
        elif maxVal == left:
            alignment1 = "-" + alignment1
            alignment2 = query[j-1] + alignment2
            # decrement the column pointer
            j -= 1
        # if we are going up that means there is a blank in the second alignment
        elif maxVal == up:
            alignment1 = ref[i-1] + alignment1
            alignment2 = "-" + alignment2
            # decrement row pointer
            i -= 1
    # print out the results
    print "ref:"
    print ref
    print
    print "query:"
    print query
    print
    print "Alignment:"
    print alignment1
    print
    print alignment2
    print
    print "score:"
    print score
    # return the score and alignments
    return score, alignment1, alignment2

# anchored_needleman_wunsch based on the regular needleman_wunsch
def anchored_needleman_wunsch(query, ref, matches):
    total_score = 0
    # for each match get the substring from query and ref and run needleman_wunsch on the result substring
    # increment total score everytime
    for match in matches:
        start1, end1, start2, end2 = match
        q = query[start1-1:end1]
        r = ref[start2-1:end2]
        score, alignment1, alignment2 = needleman_wunsch(query=q,ref=r)
        total_score += score
    # print score
    print "total score:"
    print total_score
    # return the total score
    return total_score

# helper function to read the sequence from the .fa files returns a string
def read_sequence(filename):
    with open(filename) as f:
        content = f.read().splitlines()
    return "".join(content[1:])

# helper function to read the match file and returns a list of list of indices
def read_matches(filename):
    ret = []
    with open(filename) as f:
        content = f.read().splitlines()
    for line in content:
        ret.append(map(int, re.split(r'\t+', line)))
    return ret

# helper function using numpy to permute the input sequence. returns numpy.array
def permute(sequence):
    return np.random.permutation(list(sequence))

# function to run anchored_needleman_wunsch 10000 times and uses matplotlib to show a histogram
def run_10k_times(query, ref, matches):
    results = []
    for _ in xrange(10000):
        # permute query and ref
        rand_q = permute(query)
        rand_r = permute(ref)
        # the the score from anchored_needleman_wunsch
        score = anchored_needleman_wunsch(rand_q, rand_r, matches)
        results.append(score)
    # convert result array to numpy.array
    arr = np.array(results)
    # create and show histogram and auto bins
    plt.hist(arr, bins='auto')
    plt.show()

if __name__ == '__main__':
    parser = make_arg_parser()
    args = parser.parse_args()
    query = read_sequence(args.query)
    ref = read_sequence(args.ref)
    # if there is no args then call regular needleman_wunsch
    if not args.matches:
        needleman_wunsch(query=query, ref=ref)
    else:
        matches = read_matches(args.matches)
        #uncomment to run the 10000 times random permutations
        # run_10k_times(query=query, ref=ref, matches=matches)
        anchored_needleman_wunsch(query=query, ref=ref, matches=matches)
