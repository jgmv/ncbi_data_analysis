#!/usr/bin/python

'''
    File name: get_metadata_from_BioSample.py
    Author: Jose G. Macia-Vicente
    Date created: 2019-08-06
    Python Version: 2.7
'''


from Bio import SeqIO
import re
import sys
import argparse

parser = argparse.ArgumentParser(
    description='Retrieve metadata from BioSample records as a table.')

parser.add_argument(
    "gb",
    help = 'input file with BioSample records')

parser.add_argument(
    "-o", help = 'output file',
    type=str,
    default='bsTable.csv')

args = parser.parse_args()

of = args.o

if1 = args.gb
of1 = args.o



handle = open(if1)
out = open(of1, "w")

# write header for output file
out.write("description\tbiosample\torganism\thost\tdate\tlocation\televation \
  \tlat_lon\n")

count = 0
for line in handle:
	if line.startswith("1: "):
		count +=1
		description = line.split(": ")[1].rstrip()
		host = "NA"
		date = "NA"
		location = "NA"
		elevation = "NA"
		lat_lon = "NA"
		lat = "NA"
		lon = "NA"
	if line.startswith("Identifiers: "):
		biosample = line.split(": ")[2]
		biosample = biosample.replace("; Sample name", "")
	if line.startswith("Organism: "):
		organism = line.split(": ")[1].rstrip()
	if line.startswith("    /host="):
		host = line.split("=")[1].replace("\"", "").rstrip()
	if line.startswith("    /collection date="):
		date = line.split("=")[1].replace("\"", "").rstrip()
	if line.startswith("    /geographic location="):
		location = line.split("=")[1].replace("\"", "").rstrip()
	if line.startswith("    /elevation="):
		elevation = line.split("=")[1].replace("\"", "").rstrip()
	if line.startswith("    /latitude and longitude="):
		lat_lon = line.split("=")[1].replace("\"", "").rstrip()
		if len(lat_lon.split(" ")) == 4:
			lat = lat_lon.split(" ")[0]
			lon = lat_lon.split(" ")[2]
			if lat_lon.split(" ")[1] == "S":
				lat = float(lat) * -1
			if lat_lon.split(" ")[3] == "W":
				lon = float(lon) * -1
			lat_lon = str(lat)+" "+str(lon)
	if  line.startswith("Accession: "):
		out.write(description+"\t"+biosample+"\t"+organism+"\t"+host+"\t"+date \
		  +"\t"+location+"\t"+elevation+"\t"+lat_lon+"\n")


print str(count)+" records processed"

handle.close()
out.close()
