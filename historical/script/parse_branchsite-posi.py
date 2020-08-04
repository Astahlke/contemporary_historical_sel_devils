#!/usr/bin/python
import sys
import os
import csv

from Bio.Phylo.PAML import codeml

#runparse.sh supplies the following args
#sysargv1 is directory
#sysargv2 is geneID

fname = os.path.join(sys.argv[1])
gene = sys.argv[2]
filename = sys.argv[3]

#with open(filename, 'a') as outfile:
#	outfile.write(str(gene) + '-null' +'\t' + 'lnL' + '\t' +  'branch length' + '\t' +  'prop0' + '\t ' + 'background0' + '\t' + 'foreground0'+ '\t' + '\t '+ 'prop1' + 'background1' + '\t' + 'foreground1'+  '\t' + 'prop2' + '\t '+ 'background2' + '\t' + 'foreground2'+  '\t' + 'prop3' + '\t ' 'background3' + '\t' + 'foreground3' + '\n')

results = codeml.read(fname)
lnl = results['NSsites'][2]['lnL']
treel = results['NSsites'][2]['tree length']

p0 = results['NSsites'][2]['parameters']['site classes'][0]['proportion']
b0 = results['NSsites'][2]['parameters']['site classes'][0]['branch types']['background']
f0 = results['NSsites'][2]['parameters']['site classes'][0]['branch types']['foreground']

p1 = results['NSsites'][2]['parameters']['site classes'][1]['proportion']
b1 = results['NSsites'][2]['parameters']['site classes'][1]['branch types']['background']
f1 = results['NSsites'][2]['parameters']['site classes'][1]['branch types']['foreground']

p2 = results['NSsites'][2]['parameters']['site classes'][2]['proportion']
b2 = results['NSsites'][2]['parameters']['site classes'][2]['branch types']['background']
f2 = results['NSsites'][2]['parameters']['site classes'][2]['branch types']['foreground']

p3 = results['NSsites'][2]['parameters']['site classes'][3]['proportion']
b3 = results['NSsites'][2]['parameters']['site classes'][3]['branch types']['background']
f3 = results['NSsites'][2]['parameters']['site classes'][3]['branch types']['foreground']

with open(filename, 'a') as outfile:
        outfile.write(str(gene) + '-posi' +'\t' + str(lnl) + '\t' +  str(treel) + '\t' + str(p0) + '\t' + str(b0) + '\t' + str(f0) + '\t' + str(p1) + '\t'+ str(b1) + '\t' + str(f1)+ '\t' + str(p2) + '\t'+ str(b2) + '\t' + str(f2) + '\t' + str(p3) + '\t'+ str(b3) + '\t' + str(f3) + '\n')
