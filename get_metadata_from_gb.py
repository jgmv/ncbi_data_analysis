#!/usr/bin/python

'''
    File name: get_metadata_from_gb.py
    Author: Jose G. Macia-Vicente
    Date created: 2015-05-26
    Date last modified: 2019-07-25
    Python Version: 2.7
'''


from Bio import SeqIO
import re
import sys
import argparse

parser = argparse.ArgumentParser(
    description='Retrieve metadata from GenBank records.')

parser.add_argument(
    "gb",
    help = 'input file with GenBank records')

parser.add_argument(
    "-o", help = 'output file',
    type=str,
    default='gbTable.csv')

args = parser.parse_args()

of = args.o

if1 = args.gb
of1 = args.o

count = 0

handle = open(if1)
seqOut = open(of1, "w")

# write header for output file
seqOut.write("accession|sequence_version|organism|kingdom|phylum|Class|order|family|genus|species|db_xref|strain|isolate|specimen_voucher|clone|country|lat_lon|host|isolation_source|collection_date|PCR_primers|environmental_sample|note|authors|title|journal|pubmed_id|ref_comment"+"\n")

for seq_record in SeqIO.parse(handle, "genbank"):	
	accno = seq_record.name
	try: seq_v = accno+"."+str(seq_record.annotations["sequence_version"])
	except:	seq_v = "NA"

	# extract taxonomy
	organism = "_".join(seq_record.annotations["organism"].split())
	if organism == "":
		organism = "NA"

	tax = seq_record.annotations["taxonomy"]

	try: kingdom = tax[1]
	except:	kingdom = "NA"

	try: phylum = tax[3]
	except:	phylum = "NA"

	try: Class = tax[5]
	except:	Class = "NA"

	try: order = tax[7]
	except:	order = "NA"

	try: family = tax[8]
	except:	family = "NA"

	try: genus = tax[9]
	except:	genus = "NA"

	try: species = tax[10]
	except:	species = "NA"

	# extract reference data
	ref = seq_record.annotations["references"][0]

	authors = ref.authors
	if authors == "":
		authors = "NA"

	title = ref.title
	if title == "":
		title = "NA"

	journal = ref.journal
	if journal == "":
		journal = "NA"

	pubmed = ref.pubmed_id
	if pubmed == "":
		pubmed = "NA"

	ref_comment = ref.comment
	if ref_comment == "":
		ref_comment = "NA"

	# extract features
	features = seq_record.features[0]

	db_xref = str(features.qualifiers.get("db_xref"))[2:-2]
	if db_xref == "":
		db_xref = "NA"

	strain = "_".join(str(features.qualifiers.get("strain"))[2:-2].split())
	if strain == "":
		strain = "NA"

	isolate = "_".join(str(features.qualifiers.get("isolate"))[2:-2].split())
	if isolate == "":
		isolate = "NA"

	specimen_voucher = str(features.qualifiers.get("specimen_voucher"))[2:-2]
	if specimen_voucher == "":
		specimen_voucher = "NA"

	clone = str(features.qualifiers.get("clone"))[2:-2]
	if clone == "":
		clone = "NA"

	country = str(features.qualifiers.get("country"))[2:-2]
	if country == "":
		country = "NA"

	lat_lon = str(features.qualifiers.get("lat_lon"))[2:-2]
	if lat_lon == "":
		lat_lon = "NA"

	host = str(features.qualifiers.get("host"))[2:-2]
	if host == "":
		host = "NA"

	isolation_source = str(features.qualifiers.get("isolation_source"))[2:-2]
	if isolation_source == "":
		isolation_source = "NA"

	collection_date = str(features.qualifiers.get("collection_date"))[2:-2]
	if collection_date == "":
		collection_date = "NA"

	primers = str(features.qualifiers.get("PCR_primers"))[2:-2]
	if primers == "":
		primers = "NA"

	environmental_sample = str(features.qualifiers.get("environmental_sample"))
	if environmental_sample == "None":
		environmental_sample = "FALSE"
	else:
		environmental_sample = "TRUE"

	note = str(features.qualifiers.get("note"))[2:-2]
	if note == "":
		note = "NA"

	# write to output file
	seqOut.write(
		accno+"|"+
		seq_v+"|"+
		organism+"|"+
		kingdom+"|"+
		phylum+"|"+
		Class+"|"+
		order+"|"+
		family+"|"+
		genus+"|"+
		species+"|"+
		db_xref+"|"+
		strain+"|"+
		isolate+"|"+
		specimen_voucher+"|"+
		clone+"|"+
		country+"|"+
		lat_lon+"|"+
		host+"|"+
		isolation_source+"|"+
		collection_date+"|"+
		primers+"|"+
		environmental_sample+"|"+
		note+"|"+
		authors+"|"+
		title+"|"+
		journal+"|"+
		pubmed+"|"+
		ref_comment
		+"\n")

	count +=1

print str(count)+" sequences processed"

handle.close()
seqOut.close()
