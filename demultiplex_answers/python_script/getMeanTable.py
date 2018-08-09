#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import argparse

def get_arguments():
	'''Takes a FASTQ file and reports the mean for each quality score for each position. File is required,
    and an output name is required.'''
	parser = argparse.ArgumentParser(description = "Stats for each quality score position")
	parser.add_argument("-f", "--filename", help= "Designates the FASTQ file",\
		required=True,type=str)
	parser.add_argument("-o", "--output", help= "Designates the output TSV name",\
		required=True,type=str)
	return parser.parse_args()

args = get_arguments()

file = args.filename
outputFile = args.output

scores = np.zeros(101, dtype=float)

def convert_Phred(char):
    '''Converts a character to a Phred score.'''
    return ord(char) - 33

def totalScores(file, scores):
    '''Adds up all scores based on basepair position.'''
    with open(file, "r") as fh:
        LN = 0
        for line in fh:
            LN += 1
            if LN % 4 == 0:
                line = line.strip("\n")
                for char in range(len(line)):
                    scores[char] += convert_Phred(line[char])
                if LN % 100000 == 0:
                    print(LN)
        return scores, LN

def getMean(scores, LN):
    '''Calculates the mean for each position.'''
    reads = LN / 4
    for value in range(len(scores)):
        scores[value] = scores[value] / reads
    return scores

def printStats(scores):
    '''Prints the TSV table with mean for each position in the file.'''
    with open(outputFile, "w") as fh:
        fh.write("#_Base_Pairs\tMean_Quality_Score\n")
        for position in range(len(scores)):
            resultString1 = str(position) + "\t" + str(scores[position]) + "\n"
            fh.write(resultString1)

scores, LN = totalScores(file, scores)
mean_scores = getMean(scores, LN)
printStats(mean_scores)