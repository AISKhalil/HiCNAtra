import os, sys, getopt
import glob
import logging

import numpy as np

import mirnylib.genome
import hiclib.mapping
import mirnylib.h5dict

def Genome_Mapping(argv):
	inputfile = ''
	outputfile = ''
	try:
		opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
	except getopt.GetoptError:
		print('Mapping.py -i <inputfile> -o <outputfile>')
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print('Mapping.py -i <inputfile> -o <outputfile>')
	 		sys.exit()
		elif opt in ("-i", "--ifile"):
	 		inputfile = arg
		elif opt in ("-o", "--ofile"):
	 		outputfile = arg
	logging.basicConfig(level=logging.DEBUG)

	if not os.path.exists('tmp'):
	    os.mkdir('tmp')

	# bowtie2 index file
	bowtie_bin_path   = '.../bowtie2'
	bowtie_index_path = '../HiCNAtraTool/Annotations/hg19/UCSC_chromFa/chrIndex'

	# hg19 reference genome
	fasta_path = '../HiCNAtraTool/Annotations/hg19/UCSC_chromFa'
	genome_db  = mirnylib.genome.Genome(fasta_path, readChrms=['#', 'X'])

	# input & output files
	firstFastq  = inputfile  + "_1.fastq"
	secondFastq = inputfile  + "_2.fastq"	
	firstBam    = outputfile + "_1.bam"
	secondBam   = outputfile + "_2.bam"
	outputDict  = outputfile + ".hdf5"


	############################################
	# 1. Aligning the raw reads with bowtie2 ###
	print("Align the raw reads with bowtie2...")

	hiclib.mapping.iterative_mapping(
	    bowtie_path=bowtie_bin_path,
	    bowtie_index_path=bowtie_index_path,
	    fastq_path= firstFastq,
	    out_sam_path= firstBam,
	    min_seq_len=20,
	    len_step=5,
	    seq_start=0,
	    seq_end=100,	
	    nthreads=12,
	    tmp_dir='./tmp',
	    bowtie_flags='--very-sensitive')

	hiclib.mapping.iterative_mapping(
	    bowtie_path=bowtie_bin_path,
	    bowtie_index_path=bowtie_index_path,
	    fastq_path= secondFastq,
	    out_sam_path= secondBam,
	    min_seq_len=20,
	    len_step=5,
	    seq_start=0,
	    seq_end=100,	
	    nthreads=12,
	    tmp_dir='./tmp',
	    bowtie_flags='--very-sensitive')


	############################################
	# 2. Parsing the bowtie2 output          ###
	print("Parsing the bowtie2 generated BAMs...")

	lib = mirnylib.h5dict.h5dict(outputDict)

	hiclib.mapping.parse_sam(
		sam_basename1=firstBam,
		sam_basename2=secondBam,
		out_dict=lib,
		genome_db=genome_db,
		keep_ids=False,
		enzyme_name='HindIII')

	print("Done!")

if __name__ == "__main__":
	Genome_Mapping(sys.argv[1:])
