#!/usr/bin/python

'''
    File name: get_metadata_from_gb.py
    Author: Jose G. Macia-Vicente
    Date created: 2019-08-16
    Python Version: 2.7
'''

from Bio import SeqIO
import re
import sys
import argparse

parser = argparse.ArgumentParser(
    description='Retrieve protein sequences from GenBank records.')

parser.add_argument(
    "gb",
    help = 'input file with GenBank records')

parser.add_argument(
    "-o", help = 'output file',
    type=str,
    default='outSeq.fasta')

args = parser.parse_args()

of = args.o

if1 = args.gb
of1 = args.o

count = 0

handle = open(if1)
seqOut = open(of1, "w")

for seq_record in SeqIO.parse(handle, "genbank"):
	features = seq_record.features[3]
	accno = seq_record.name
	organism = "_".join(seq_record.annotations["organism"].split())
	if features.qualifiers.get("translation"):
		prot = str(features.qualifiers.get("translation")[0].lower())
		seqOut.write(">"+accno+"_"+organism+"\n")
		seqOut.write(prot+"\n")
		count +=1
	else:
		print "No protein sequence for "+accno+"_"+organism

print str(count)+" sequences processed"

handle.close()
seqOut.close()
