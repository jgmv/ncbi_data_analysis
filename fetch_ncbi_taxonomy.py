#!/usr/bin/python3

'''
    Author: Jose G. Macia-Vicente
    Date created: 2023-05-10
    Date last modified: 2023-05-10
    Function: fetches taxids and taxonomic lineages for taxa 
'''

from ete3 import NCBITaxa
import argparse
import os

# arguments
parser = argparse.ArgumentParser(
    description='Obtains NCBI taxids and lineages from a list of taxa.')
parser.add_argument(
    "i",
    help = 'input file with taxonomic names')
parser.add_argument(
    "-o", help = 'output file',
    type=str,
    default='output.tsv')
parser.add_argument(
    "-f", help = 'proportion of word similarity for inexact name search',
    type=float,
    default=1.0)
parser.add_argument(
    "-c", help = 'clean-up taxonomic database?',
    action=argparse.BooleanOptionalAction)

# parse arguments
args = parser.parse_args()
handle = open(args.i, "r")
out = open(args.o, "w")

# reset number of entries
count = 0

# fetch taxonomic database from NCBI
ncbi = NCBITaxa()

# get taxonomic info for each entry
if args.f == 1:
    out.write("query\ttaxid\tkingdom\tphylum\tclass\torder\tfamily\n")
else:
    out.write("query\tfuzzy_match\ttaxid\tkingdom\tphylum\tclass\torder\tfamily\n")

for record in handle:
    record = record.strip()
    if args.f == 1:
        try:
            taxid = ncbi.get_name_translator([record])
            taxid = str(taxid[record]).replace("]", "").replace("[", "")
        except:
            print(record+" not found in database.")
            taxid = 1
    else:
        fuzzy = ncbi.get_fuzzy_name_translation(record, args.f)
        taxid = fuzzy[0]
        if taxid == None:
            taxid = 1
        match = str(fuzzy[1])+"|"+str(fuzzy[2])
    lineage = ncbi.get_lineage(taxid)
    ranks = ncbi.get_rank(lineage)
    kingdom = "unknown"
    phylum = "unknown"
    classs = "unknown"
    order = "unknown"
    family = "unknown"
    for i in ranks:
        if ranks[i] == "kingdom":
            kingdom = ncbi.get_taxid_translator([i])
            kingdom = next(iter(kingdom.values()))
        elif ranks[i] == "phylum":
            phylum = ncbi.get_taxid_translator([i])
            phylum = next(iter(phylum.values()))
        elif ranks[i] == "class":
            classs = ncbi.get_taxid_translator([i])
            classs = next(iter(classs.values()))
        elif ranks[i] == "order":
            order = ncbi.get_taxid_translator([i])
            order = next(iter(order.values()))
        elif ranks[i] == "family":
            family = ncbi.get_taxid_translator([i])
            family = next(iter(family.values()))
    if args.f == 1:
        out.write(record+"\t"+str(taxid)+"\t"+kingdom+"\t"+phylum+"\t"+classs+"\t"+order+"\t"+family+"\n")
    else:
        out.write(record+"\t"+match+"\t"+str(taxid)+"\t"+kingdom+"\t"+phylum+"\t"+classs+"\t"+order+"\t"+family+"\n")
    count +=1
print(str(count)+" queries processed.")

out.close()

# clean-up environment
if args.c:
    os.system("rm taxdump.tar.gz")

# end
