#!/usr/bin/env python

'''
Title: add_bq_from_fastq.py
Date: 20221221
Author: Adam Nunn
Description:
	This program takes an input BAM file alongside corresponding FASTQ files
	and appends the original base quality information to the appropriate
	alignments. To be used for example in the case of BS-Seeker2 which neglects
	to provide the appropriate base quality information in output BAM files.  

	++ R1 forward -> same as R1 file
	+- R1 reverse -> revcomp R1 file
	-+ R2 forward -> same as R2 file
	-- R2 reverse -> revcomp R2 file

List of functions:
	main()
	build_fastq()


Procedure:
	1. Try to build the fastq dictionary if 'MATE1' and 'MATE2' are given
	2. Open pysam.AlignmentFile objects for 'BAM' and 'OUT'
	3. Iterate through each alignment in 'BAM' pysam.AlignmentFile object
	4. Determine whether the sequence is forward or reverse strand
	5. Append original base qualities from 'MATE1' and/or 'MATE2'
	6. Write the new, modified alignment to the 'OUT' pysam.AlignmentFile object


Usage:
		./add_bq_from_fastq.py [-h, --help] \
		[-1, --mate1] mate_1.fastq \
		[-2, --mate2] mate_2.fastq \
		BAM OUT

eg. ./add_bq_from_fastq.py -1 mate_1.fastq -2 mate_2.fastq input.bam output.bam
'''

###################
## INIT ENVIRONMENT

import argparse
import pysam
from array import array


##################
## DEFINE __MAIN__
def main(BAM,OUT,MATE1=None,MATE2=None):

	# Build dictionaries from fastq files
	mate1 = build_fastq(MATE1)
	try: mate2 = build_fastq(MATE2)
	except (TypeError):
		mate2 = None
		pass

	# Open pysam.AlignmentFile objects for reading and writing
	with pysam.AlignmentFile(BAM, "rb") as original, pysam.AlignmentFile(OUT, "wb", header=original.header) as modified:

		# iterate over each alignment from 'original'
		for alignment in original:

			# filter out unmapped, secondary alignments, or qcfailed alignments
			if (alignment.is_unmapped) or (alignment.is_secondary) or (alignment.is_qcfail): continue
			if (alignment.query_alignment_length == 0): continue

			# get alignment variables
			qnam = alignment.query_name
			qseq = alignment.query_sequence
			qual = alignment.query_qualities

			# skip erroneous reads
			if (qnam == None) or (qseq == None): continue

			# base quality string not present in bam file
			if qual == None:

				# reverse complement
				if alignment.is_reverse:

					if (alignment.is_paired and alignment.is_read1) or not alignment.is_paired:
						if qnam in mate1: alignment.query_qualities = mate1[qnam][::-1]
						else: print("WARN: {} present in bam but missing from mate1 file".format(qnam))
						
					else:
						if qnam in mate2: alignment.query_qualities = mate2[qnam][::-1]
						else: print("WARN: {} present in bam but missing from mate2 file".format(qnam))	

				# not reverse complement
				else:

					if (alignment.is_paired and alignment.is_read1) or not alignment.is_paired:
						if qnam in mate1: alignment.query_qualities = mate1[qnam]
						else: print("WARN: {} present in bam but missing from mate1 file".format(qnam))
					else:
						if qnam in mate2: alignment.query_qualities = mate2[qnam][::-1]
						else: print("WARN: {} present in bam but missing from mate2 file".format(qnam))	
				

			# write new read to 'modified' file
			modified.write(alignment)

	###################
	# Goodbye message
	print("\n----------------")
	print("Success!\n")


## END OF __MAIN__
##################


###################
## DEFINE FUNCTIONS

####### Function to build 'fastq' dictionary from 'FASTQ' file object 
def build_fastq(FASTQ):

	######################################################################################################
	## FASTQ = string of path to input fastq eg. "/path/to/mate_1.fastq"                                ##
	######################################################################################################

	# declare an empty dictionary and stage the 'FASTQ' file
	dfastq = dict()
	with open(FASTQ, 'r') as fastq:
		count = 0

		# iterate through the lines of the 'FASTQ' file
		for line in fastq:
			line = line.rstrip()
			count += 1

			# identify the sequence ID
			if count % 4 == 1:
				line = line.split(" ")
				ID = line[0][1:]

			# identify base quality string
			elif count % 4 == 0:
				qual = array('B',[])
				for bq in line.upper(): qual.append(ord(bq)-33)
				dfastq[ID] = qual

			# skip other lines
			else: continue

	# return the constructed genome dictionary
	return dfastq

## END OF FUNCTIONS
###################

#############
## RUN SCRIPT

# define argparse
usage = ''' This program takes an input BAM file alongside corresponding FASTQ files
	and appends the original base quality information to the appropriate
	alignments. To be used for example in the case of BS-Seeker2 which neglects
	to provide the appropriate base quality information in output BAM files. '''

parser = argparse.ArgumentParser(description=usage)

parser.add_argument('inbam', metavar='<BAM>', help='[REQUIRED] Path to input BAM file')
parser.add_argument('outbam', metavar='<OUT>', help='[REQUIRED] Path to output BAM file')
parser.add_argument('-1','--mate1', metavar='<MATE1>', help='[REQUIRED] Path to input MATE1 fastq file.')
parser.add_argument('-2','--mate2', metavar='<MATE2>', help='[OPTIONAL] Path to input MATE2 fastq file. Not required in case of single-end sequencing.')

args = parser.parse_args()
parser.parse_args()

# call main()
if __name__ == '__main__':
	main(args.inbam,args.outbam,args.mate1,args.mate2)

## END OF SCRIPT
################