#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 30 23:23:59 2020

@author: korolaxyw
"""
import random
import collections
import matplotlib.pyplot as plt

Nucleotides = ["A","T","C","G"]
DNAComplement = {'A':'T', 'T':'A', 'C':'G','G':'C'}
Amino_Acid_Codons = {
    "TTT": "Phe", "TTC": "Phe",
    "TTA": "Leu", "TTG": "Leu", "CTT": "Leu", "CTC": "Leu", "CTA": "Leu", "CTG": "Leu",
    "ATT": "Ile", "ATC": "Ile", "ATA": "Ile",
    "GTT": "Val", "GTC": "Val", "GTA": "Val", "GTG": "Val",
    "TCT": "Ser", "TCC": "Ser", "TCA": "Ser", "TCG": "Ser", "AGT": "Ser", "AGC": "Ser",
    "CCT": "Pro", "CCC": "Pro", "CCA": "Pro", "CCG": "Pro",
    "ACT": "Thr", "ACC": "Thr", "ACA": "Thr", "ACG": "Thr",
    "GCT": "Ala", "GCC": "Ala", "GCA": "Ala", "GCG": "Ala",
    "TAT": "Tyr", "TAC": "Tyr",
    "CAT": "His", "CAC": "His",
    "CAA": "Gln", "CAG": "Gln",
    "AAT": "Asn", "AAC": "Asn",
    "AAA": "Lys", "AAG": "Lys",
    "GAT": "Asp", "GAC": "Asp",
    "GAA": "Glu", "GAG": "Glu",
    "TGT": "Cys", "TGC": "Cys",
    "TGG": "Trp",
    "CGT": "Arg", "CGC": "Arg", "CGA": "Arg", "CGG": "Arg", "AGA": "Arg", "AGG": "Arg",
    "GGT": "Gly", "GGC": "Gly", "GGA": "Gly", "GGG": "Gly",
    "ATG": "Met",
    "TAA": "Stop", "TAG": "Stop", "TGA": "Stop",
    }

fasta={}

with open('sequence.fasta') as f:
      sequence = ""  
      for line in f:
          if line.startswith(">"):
            name=line[1:].rstrip()
            fasta[name]=""
            continue
          fasta[name]+=line.strip().upper()

Covid_19_Sequence = fasta['COVID']

rndDNAStr = ''.join([random.choice(Nucleotides)
                     for nuc in range(50)])

# Check the sequence to make sure it is DNA string containing only the four nucleosides
def validateSeq(dna_seq):
    tmpseq = dna_seq.upper()
    for nuc in tmpseq:
        if nuc not in Nucleotides:
            return False
    return tmpseq

DNAStr = validateSeq(rndDNAStr)
covid19_DNAStr = validateSeq(Covid_19_Sequence)
#print(DNAStr)

# Count the frequency of occurence of each of the four bases in the DNA sequence string
def countNucFrequency(seq):
    tmpFreqDict = {"A": 0, "T": 0, "C": 0, "G": 0}
    for nuc in seq:
        tmpFreqDict[nuc] += 1
    return tmpFreqDict
    
#print(countNucFrequency(DNAStr))

# Transcription from DNA to RNA
def transcription(seq):
    transcript = seq.replace("T","U")
    return transcript

#print(transcription(DNAStr))

# Reverse Complement of the DNA single strand via base pairing
def DNAcomplement(seq):
    Basepairing_strand = ''.join([DNAComplement[nuc] for nuc in seq])[::1]
    return Basepairing_strand

rDNAStr = DNAcomplement(DNAStr)
r_covid19_DNAStr = DNAcomplement(covid19_DNAStr)

# Each of the four bases in DNA sequence will be assigned a unique color when the result is printed in terminal
def coloring(seq):
    colors = {
        'A': '\033[92m',
        'C': '\033[94m',
        'G': '\033[93m',
        'T': '\033[91m',
        'U': '\033[91m',
        'reset': '\033[0;0m'
    }
    
    tmpStr = ""
    for nuc in seq:
        if nuc in colors:
            tmpStr += colors[nuc] + nuc
        else:
            tmpStr += colors['reset'] + nuc
            
    return tmpStr + '\033[0;0m'

# Count the proportion of numbers of G and C in the DNA sequence
def gc_portion(seq):
    gccontent = round((seq.count('C')+seq.count('G'))/len(seq)*100)
    return gccontent

#print(gc_portion(DNAStr))

# Translate DNA sequence to Amino Acid Codons so the protein sequence can be analyzed
def translation(seq, starting_position):
    amino_acids = [Amino_Acid_Codons[seq[position:position + 3]] for position in range(starting_position, len(seq) - 2, 3)]
    return amino_acids

# Count the number of the amino acid(s) of interest appeared in the DNA sequence
def aminoacidFreq(seq, aminoacid):
    freqlist = []
    for i in range(0, len(seq) - 2, 3):
        if Amino_Acid_Codons[seq[i:i + 3]] == aminoacid:
            freqlist.append(seq[i:i +3])
    return freqlist

#print(aminoacidFreq(covid19_DNAStr, "Lys"))

# Analyze the frequency of different nucleotides combinations that contribute to the same amino acid codon of interest
def countCodonFrequency(seq, aminoacid):
    tempCodonlist = []
    # the codon is translated from position 0, and every 3 nucleotides form 1 amino acid
    for i in range(0, len(seq) - 2, 3):
        if Amino_Acid_Codons[seq[i:i + 3]] == aminoacid:
            tempCodonlist.append(seq[i:i +3])
    aafreqDict = dict(collections.Counter(tempCodonlist))
    total_weight = sum(aafreqDict.values())
    for seq in aafreqDict:
        aafreqDict[seq] = round(aafreqDict[seq]/ total_weight, 2)
    return aafreqDict

# Read the translation of DNA sequence in 3 different fram shifts to find out the best ORF
def reading_frame(seq):
    frame = []
    frame.append(translation(seq, 0))
    # then shift the reading frame by 1 nucleotide
    frame.append(translation(seq, 1))
    # then shift the reading frame by 2 nucleotides
    frame.append(translation(seq, 2))
    return frame

# List of proteins that translate between 'Met'(included) and 'Stop'(not included) because we only care about ORF
def ORF_protein(aa_seq): #aa represents amino_acid
    current_protein = []
    orf_proteins = []
    for aa in aa_seq:
        if aa == "Stop":
            if current_protein:
                for s in current_protein:
                    orf_proteins.append(s)
                current_protein = []
        else:
            if aa == "Met":
                current_protein.append("")
            for i in range(len(current_protein)):
                current_protein[i] += aa
    return orf_proteins

        
#print(f'\nSequence: {coloring(covid19_DNAStr)}\n')
#print(f'[1] : Sequence Length: {len(covid19_DNAStr)}\n')
#print(f'[2] : Nucleotide Frequencies: {countNucFrequency(covid19_DNAStr)}\n')
#print(f'[3] : DNA Single Strand DNA Complement:{DNAcomplement(covid19_DNAStr)}\n')
#print(f"[4] : Covid19 DNA Complentary Sequence: 3' {DNAcomplement(covid19_DNAStr)} 5'\n")
#print(f'[5] : Transcription Result: {transcription(r_covid19_DNAStr)}\n')
#print(f'[6] : GC_Portions of DNA Sequence: {gc_portion(covid19_DNAStr)}%\n')
#print(f'[7] : Translation Result: {translation(covid19_DNAStr, 0)}\n')
#print(f'[8] : Lys Count: {len(aminoacidFreq(covid19_DNAStr, "Lys"))}\n')
#print(f'[9] : Amino Acid Frequency (Lys): {countCodonFrequency(covid19_DNAStr, "Lys")}\n')
#print(f'[10] : Cys Count: {len(aminoacidFreq(covid19_DNAStr, "Cys"))}\n')  
#print(f'[11] : Amino Acid Frequency (Cys): {countCodonFrequency(covid19_DNAStr, "Cys")}\n')
#print(f'[12] : Tyr Count: {len(aminoacidFreq(covid19_DNAStr, "Tyr"))}\n')
#print(f'[13] : Amino Acid Frequency (Tyr): {countCodonFrequency(covid19_DNAStr, "Tyr")}\n')
#print('[14] : Reading Frames:')    
#for frame in reading_frame(covid19_DNAStr):
   # print(frame)
#print('[15] : All Open Reading Frames:')    
#for frame in reading_frame(covid19_DNAStr):
    #print(ORF_protein(frame))

    
           
