#!/usr/bin/python

'''
    File name: parse_taxids.py
    Author: Jose G. Macia-Vicente
    Date created: 2020-04-10
    Date last modified: 2020-04-10
    Python Version: 3.7
'''


from ete3 import NCBITaxa
import re
import sys
import argparse


parser = argparse.ArgumentParser(
    description='Obtains lineages from a list of taxids.')

parser.add_argument(
    "i",
    help = 'input file with taxid list')

parser.add_argument(
    "-o", help = 'output file',
    type=str,
    default='output_taxa.csv')


args = parser.parse_args()


if1 = args.i
of1 = args.o

count = 0

handle = open(if1, "r")
taxOut = open(of1, "w")
ncbi = NCBITaxa()

for record in handle:
	lineage = ncbi.get_lineage(record)
	names = ncbi.get_taxid_translator(lineage)
	taxOut.write(record.strip())
	for taxid in lineage:
		taxOut.write(";"+names[taxid])
	taxOut.write("\n")
	count +=1

print(str(count))+" taxids read"

handle.close()
taxOut.close()
