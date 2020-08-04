#!/usr/bin/python
import sys

#This program takes three inputs: the alignment in phylip format, a treefile in Newick format, and outputfile name
#It returns a nested dictionary with various results depending on analysis
from Bio.Phylo.PAML import codeml

cml = codeml.Codeml()
cml.alignment = sys.argv[1]
cml.tree = sys.argv[2]
cml.out_file = sys.argv[7]
cml.set_options(verbose=1)
cml.set_options(seqtype=1) # codon 
cml.set_options(model=sys.argv[3]) # *0:one, 1:b, 2:2 or more dN/dS ratios for branches
cml.set_options(fix_omega=sys.argv[4]) # * fix omega
cml.set_options(omega=sys.argv [5]) # *if fixed, omega
cml.set_options(cleandata = 0)
cml.set_options(CodonFreq=sys.argv[6]) # codon_freq
cml.set_options(NSsites= [2]) #1/:neutral;2:selection; 3:discrete;4:freqs;
 # 5:gamma;6:2gamma;7:beta;8:beta&w;9:beta&gamma;
 # 10:beta&gamma+1; 11:beta&normal>1; 12:0&2normal>1;
 # 13:3normal>0 
cml.working_dir = "."

results = cml.run(verbose=True)

#`results = yn.read(sys.argv[2])

