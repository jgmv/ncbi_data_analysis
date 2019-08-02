#!/usr/bin/python

'''
    File name: fetch_gb.py
    Version: 2.0
    Author: Jose G. Macia-Vicente
    Date created: 2016-07-01
    Date last modified: 2019-08-02
    Python Version: 2.7
'''

__author__ = "Jose G. Macia-Vicente"

from Bio import Entrez
from Bio import SeqIO
import sys
import argparse

# input arguments parameters
parser = argparse.ArgumentParser(
    description='Retrieve GenBank records from NCBI from GI list.')

parser.add_argument(
    "gi",
    help = 'input file with GI list')

parser.add_argument(
    "-o", help = 'output file',
    type=str,
    default='gbList.txt')

parser.add_argument(
    "-email", help = 'tell NCBI who you are!',
	type=str,
    default = 'your@email')

args = parser.parse_args()

giList = args.gi

of = args.o

# process query
Entrez.email = args.email

id_list = set(open(giList, 'rU'))

id_list = list(id_list) 

# save results to output file
giOut = open(of, "w")

count = 0
for gi in id_list:
	handle = Entrez.efetch(db="nuccore", id=gi, rettype="gb", retmode="text")
	giOut.write(handle.read())
	count +=1
	if count % 100 == 0:
		print str(count)+" records processed..." 

giOut.close()
handle.close()

