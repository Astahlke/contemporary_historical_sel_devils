import sys
import os

inputfile = open(sys.argv[1], 'r')
fasta = inputfile.read()

length = fasta.split(">")
nuc = fasta.split("\n")

i = 1
seqlength = {}
while i< len(length):
    #print(len(pci_nuc[i]))
    #print(pci_len[i])
    start = int(length[i].split()[1].split('[')[1])
    end = int(length[i].split()[3].split(']')[0])
    #print(end-start)
    seqlength[i] = (end-start)
    i +=1

k=len(length[sorted(seqlength, key=seqlength.__getitem__, reverse=True)[0]].split())
aa_seq=''.join(length[sorted(seqlength, key=seqlength.__getitem__, reverse=True)[0]].split()[4:k])
loc_name = length[sorted(seqlength, key=seqlength.__getitem__, reverse=True)[0]].split('\n')[0].split()[0]
print('>gene' + loc_name[6:9] + '\n' + aa_seq)

# filename=''.join()

# with open(filename, 'a') as outfile:
#         outfile.write('>' + loc_name '\n' + aa_seq)
