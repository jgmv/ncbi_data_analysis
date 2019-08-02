#!/usr/bin/python

'''
    File name: search_ncbi_gi_by_term.py
    Author: Jose G. Macia-Vicente
    Date created: 2016-06-29
    Date last modified: 2019-07-26
    Python Version: 2.7
'''

__author__ = "Jose G. Macia-Vicente"

from Bio import Entrez
import sys
import argparse

# input arguments parameters

parser = argparse.ArgumentParser(
    description='Retrieve GI list from NCBI from Entrez terms.')

parser.add_argument(
    "query",
    help = 'input string with Entrez query (with quotations)',
    type=str)

parser.add_argument(
    "-o", help = 'output file',
    type=str,
    default='giList.txt')

args = parser.parse_args()

query = args.query

of = args.o

# process query
Entrez.email = "MaciaVicente@em.uni-frankfurt.de"
handle = Entrez.esearch(
    db="nucleotide",
    term=query,
    retmax=2000000)

record = Entrez.read(handle)

# save results to output file
giOut = open(of, "w")

for i in record["IdList"]:
    giOut.write(i+"\n")

print record["Count"] + " records extracted."

giOut.close()
handle.close()
