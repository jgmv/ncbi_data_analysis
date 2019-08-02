#!/usr/bin/python

'''
    File name: get_metadata_from_gb.py
    Author: Jose G. Macia-Vicente
    Date created: 2014-11-18
    Date last modified: 2019-07-25
    Python Version: 2.7
'''

from Bio import SeqIO
import re
import sys
import argparse

parser = argparse.ArgumentParser(
    description='Retrieve fasta sequences from GenBank records.')

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
	accno = seq_record.name
	organism = "_".join(seq_record.annotations["organism"].split())
	seqOut.write(">"+accno+"_"+organism+"\n")
	seqOut.write(str(seq_record.seq.lower())+"\n")		

	count +=1

print str(count)+" sequences processed"

handle.close()
seqOut.close()
