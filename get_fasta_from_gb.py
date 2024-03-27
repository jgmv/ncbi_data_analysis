#!/usr/bin/python

'''
    File name: get_metadata_from_gb.py
    Author: Jose G. Macia-Vicente
    Date created: 2014-11-18
    Date last modified: 2024-03-27
    Python Version: 3.10.1
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

parser.add_argument(
    "-accession", help = 'output sequence accession ("acc", default) or strain ("str")',
    type=str,
    default='acc')

args = parser.parse_args()

of = args.o

if1 = args.gb
of1 = args.o
seqname = args.accession

count = 0

handle = open(if1)
seqOut = open(of1, "w")

if seqname != 'acc' and seqname != 'str':
	print("'accession' must be 'acc' or 'str'")
	sys.exit()

for seq_record in SeqIO.parse(handle, "genbank"):
	accno = seq_record.name
	organism = "_".join(seq_record.annotations["organism"].split())
	if seqname == 'acc':
		seqOut.write(">"+accno+"_"+organism+"\n")
	else:
		features = seq_record.features[0]
		strain = "_".join(str(features.qualifiers.get("strain"))[2:-2].split())
		if strain == "":
			strain = "NA"
		seqOut.write(">"+organism+"_"+strain+"\n")		
	seqOut.write(str(seq_record.seq.lower())+"\n")		

	count +=1

print(str(count)+" sequences processed")

handle.close()
seqOut.close()
