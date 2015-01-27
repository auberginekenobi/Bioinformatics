# -*- coding: utf-8 -*-
"""
Created on Mon Jan 26 14:46:03 2015

@author: chapmano
"""
import copy
import numpy

#Problem 1
def variance(vals):
    """
    input: list of numerical vals.
    output: the variance.
    """
    mean=numpy.mean(vals)
    sdsum=0
    for val in vals:
        sub=val-mean
        sd=sub**2
        sdsum+=sd
    return sdsum/float(len(vals))

def write():    
    f=open("ChoData.txt",'r')
    w=open("ChoVariance.txt",'w')
    w.write(f.readline())#write the title
    
    cutoff = float(raw_input("Pick a number between 0 and 100: "))
    
    nameline={}
    namevar=[]
    
    for line in f:
        line=line.strip()
        words=line.split('\t')
        name=words[0]
        nameline[name]=line
        vals=map(int,words[1:])
        var=variance(vals)
        namevar.append((var,name))
    namevar.sort(reverse=True)
    numtokeep=int((100-cutoff)*len(namevar)/100)
    print numtokeep
    namevar=namevar[:numtokeep]
    for (var,name) in namevar:
        w.write(nameline[name])
        print nameline[name]
    f.close()
    w.close()
    
#Problem 2
def translate(nucSeq, genCode):
    if type(nucSeq) != str:
        return "Error: Input sequence must be a string.\n"
    if len(nucSeq) % 3 != 0:
        return "Error: Input sequence must be multiple of 3.\n"
    aaSeq = ''
    nucSeq = nucSeq.upper()
    for i in range(0,len(nucSeq),3):
        codon = nucSeq[i:i+3]
        if codon in genCode:
            aaSeq = aaSeq + genCode[codon]
        else:
            return "Error: Sequence must only contain A, C, G, or T.\n"
    return aaSeq   
    
    
infile=raw_input('Enter the file with nucleotide sequences: ')
outfile=raw_input('Enter the output file name:')
code=raw_input("""Select the genetic code (default=code 1):
1. Standard Code
2. The Vertebrate Mitochondrial Code
3. The Yeast Mitochondrial Code
6. The Ciliate Code
Your Selection: """)
standardCode = {'AAA': 'K', 'AAC': 'N', 'AAG': 'K', 'AAT': 'N',
           'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
           'AGA': 'R', 'AGC': 'S', 'AGG': 'R', 'AGT': 'S',
           'ATA': 'I', 'ATC': 'I', 'ATG': 'M', 'ATT': 'I',
           'CAA': 'Q', 'CAC': 'H', 'CAG': 'Q', 'CAT': 'H',
           'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
           'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
           'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
           'GAA': 'E', 'GAC': 'D', 'GAG': 'E', 'GAT': 'D',
           'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
           'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
           'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
           'TAA': '*', 'TAC': 'Y', 'TAG': '*', 'TAT': 'Y',
           'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
           'TGA': '*', 'TGC': 'C', 'TGG': 'W', 'TGT': 'C',
           'TTA': 'L', 'TTC': 'F', 'TTG': 'L', 'TTT': 'F'}
vertMitoCode = copy.deepcopy(standardCode)
vertMitoCode['AGA'] = '*'
vertMitoCode['AGG'] = '*'
vertMitoCode['ATA'] = 'M'
vertMitoCode['TGA'] = 'W'
yeastMitoCode = copy.deepcopy(standardCode)
yeastMitoCode['ATA'] = 'M'
yeastMitoCode['CTT'] = 'T'
yeastMitoCode['CTC'] = 'T'
yeastMitoCode['CTA'] = 'T'
yeastMitoCode['CTG'] = 'T'
yeastMitoCode['TGA'] = 'W'
yeastMitoCode['CGA'] = 'X'
yeastMitoCode['CGG'] = 'X'
ciliateCode = copy.deepcopy(standardCode)
ciliateCode['TAA'] = 'Q'
ciliateCode['TAG'] = 'Q'





           
codeD = {1: standardCode, 2: vertMitoCode, 3: yeastMitoCode, 6: ciliateCode}
