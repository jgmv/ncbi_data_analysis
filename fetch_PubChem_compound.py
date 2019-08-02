#!/usr/bin/python

'''
    File name: fetch_PubChem_compound.py
    Author: Jose G. Macia-Vicente
    Date created: 2017-02-28
    Date last modified: 2017-02-28
    Python Version: 2.7
'''

__author__ = "Jose G. Macia-Vicente"

import argparse
from pubchempy import *

# input arguments parameters

parser = argparse.ArgumentParser(
    description='Retrieve PubChem records from CID list.')

parser.add_argument(
    "cid",
    help = 'input file with CID list')

parser.add_argument(
    "-o", help = 'output file',
    type=str,
    default='output.csv')

args = parser.parse_args()

# process query
handle = open(args.cid)

# save results to output file
of = open(args.o, "w")

of.write("cid|name|molecular_formula|molecular_weight|canonical_smiles|isomeric_smiles|iupac_name|exact_mass\n")
for line in handle:
    try:
      c = Compound.from_cid(line)
      of.write(str(c.cid) + "|"
        + str(c.synonyms[0]) + "|"
        + str(c.molecular_formula) + "|"
        + str(c.molecular_weight) + "|"
        + str(c.canonical_smiles) + "|"
        + str(c.isomeric_smiles) + "|"
        + str(c.iupac_name) + "|"
        + str(c.exact_mass) + "\n")
    except Exception:
      pass

of.close()
handle.close()

