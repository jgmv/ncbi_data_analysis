#!/usr/bin/python

'''
    File name: draw_chem.py
    Author: Jose G. Macia-Vicente
    Date created: 2017-03-03
    Date last modified: 2025-03-20
'''

import argparse
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem import Draw
from rdkit.Chem.Draw import MolDrawing
from rdkit.Chem.Draw.MolDrawing import DrawingOptions
import rdkit

# input arguments parameters
parser = argparse.ArgumentParser(
    description='Draw simple chemical structures from SMILEs.')

parser.add_argument(
    "-s",
    type=str,
    default="",
    help = 'input single SMILEs in-line')

parser.add_argument(
    "-f", help = 'input file with names and SMILEs (tab-separated) in rows',
    type=str,
    default='smiles.txt')

parser.add_argument(
    "-o", help = 'output file',
    type=str,
    default='output.svg')

parser.add_argument(
    "-size", help = 'size in pixels',
    type=int,
    default=300)

parser.add_argument(
    "--stereo", help = 'show stereo bonds',
    action=argparse.BooleanOptionalAction)

args = parser.parse_args()

# set options (not working)
#opt = Draw.MolDrawing.DrawingOptions()

# define drawings
def draw_mol(smile):
    d = Chem.MolFromSmiles(smile)
    rdDepictor.Compute2DCoords(d)
    rdDepictor.StraightenDepiction(d)
    return(d)
    
# process query
if args.s == "":
    handle = open(args.f)
    for line in handle:
        name = line.split()[0]
        smile = draw_mol(line.split()[1])
        Draw.MolToFile(smile, name + '.svg', size = (args.size, args.size), \
            wedgeBonds = args.stereo)
else:
    smile = draw_mol(args.s)
    Draw.MolToFile(smile, args.o, size = (args.size, args.size), \
        wedgeBonds = args.stereo)

# end
