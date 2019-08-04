#!/usr/bin/python

'''
    File name: draw_chem.py
    Author: Jose G. Macia-Vicente
    Date created: 2017-03-03
    Date last modified: 2017-03-03
    Python Version: 2.7
'''

import argparse
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Draw

# input arguments parameters

parser = argparse.ArgumentParser(
    description='Draw chemical structures from provided file.')

parser.add_argument(
    "cid",
    help = 'input CID from list')

parser.add_argument(
    "-i", help = 'input file',
    type=str,
    default='smiles.txt')

parser.add_argument(
    "-o", help = 'output file',
    type=str,
    default='output.pdf')

args = parser.parse_args()

# process query
handle = open(args.i)

# for black and white (comment for color):
opt = Draw.DrawingOptions()
from collections import defaultdict
dd = defaultdict(lambda:(0,0,0))  # <- use a default color of black
opt.elemDict = dd

# draw compound
for line in handle:
  if line.split()[0] in args.cid:
    smile = Chem.MolFromSmiles(line.split()[1])
    #smile = Chem.AddHs(smile)
    Chem.Compute2DCoords(smile)
    #Chem.EmbedMolecule(smile)
    #Chem.MMFFOptimizeMolecule(smile)
    #Chem.UFFOptimizeMolecule(smile)
    #smile = Chem.RemoveHs(smile)
    Draw.MolToFile(smile, args.o, options = opt)


