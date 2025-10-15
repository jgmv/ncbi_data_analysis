#!/usr/bin/python

'''
    File name: get_metadata_from_gb.py
    Author: Jose G. Macia-Vicente
    Date created: 2014-11-18
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
    "-accession", help = 'output sequence accession ("acc", default), strain ("str"), or isolate ("isol")',
    type=str,
    default='acc')

parser.add_argument(
    "-taxon", help = 'add organism name to sequence id ("yes", default; "no")',
    type=str,
    default='yes')

parser.add_argument(
    "-length", help = 'add sequence lengthto sequence id ("yes", default; "no")',
    type=str,
    default='yes')

args = parser.parse_args()

of = args.o

if1 = args.gb
of1 = args.o
seqname = args.accession
taxon = args.taxon
length = args.length

count = 0

handle = open(if1)
seqOut = open(of1, "w")

if seqname != 'acc' and seqname != 'str' and seqname != 'isol':
	print("'accession' must be 'acc', 'str' or 'isol'")
	sys.exit()

for seq_record in SeqIO.parse(handle, "genbank"):
	accno = seq_record.name
	organism = "_".join(seq_record.annotations["organism"].split())
	features = seq_record.features[0]
	if seqname == 'acc':
		label = accno
	elif seqname == 'isol':
		label = "_".join(str(features.qualifiers.get("isolate"))[2:-2].split())
	else:
		label = "_".join(str(features.qualifiers.get("strain"))[2:-2].split())
	if label == "":
		label = "NA"
	if taxon == 'yes':
		label = label+"_"+organism
	if length == 'yes':
		label = label+" "+str(len(seq_record.seq))
	seqOut.write(">"+label+"\n")	
	seqOut.write(str(seq_record.seq.lower())+"\n")		

	count +=1

print(str(count)+" sequences processed")

handle.close()
seqOut.close()
