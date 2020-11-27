#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 30 22:50:16 2020

@author: korolaxyw
"""

from DNAToolkit import *
import random

rndDNAStr = ''.join([random.choice(Nucleotides)
                     for nuc in range(50)])

DNAStr = validateSeq(rndDNAStr)

print(countNucFrequency(DNAStr))

print(transcription(DNAStr))

rDNAStr = DNAcomplement(DNAStr)
print(rDNAStr)

print(f'\nSequence: {coloring(DNAStr)}\n')
print(f'[1] : Sequence Length: {len(DNAStr)}\n')
print(f'[2] : Nucleotide Frequencies: {countNucFrequency(DNAStr)}\n')
print(f'[3] : DNA Single Strand DNA Complement:{DNAcomplement(DNAStr)}\n')
print(f"[4] : Original DNA Strand + DNA Complemnt:\n5' {DNAStr} 3'")
print(f"   {''.join(['|' for a in range(len(DNAStr))])}")
print(f"3' {DNAcomplement(DNAStr)} 5'\n")
print(f'[5] : Transcription Result: {transcription(rDNAStr)}\n')
print(f"[6] : Original DNA Strand + Transcription Result:\n5' {DNAStr} 3'")
print(f"   {''.join(['|' for a in range(len(DNAStr))])}")
print(f"3' {transcription(rDNAStr)} 5'\n")
print(f'[7] : GC_Portions of DNA Sequence: {gc_portion(DNAStr)}%\n')
print(f'[8] : Translation Result: {translation(DNAStr, 0)}\n')
print(f'[9] : Lys Count: {len(aminoacidFreq(DNAStr, "Lys"))}\n')
print(f'[10] : Amino Acid Frequency (Lys): {countCodonFrequency(DNAStr, "Lys")}\n')
print(f'[11] : Cys Count: {len(aminoacidFreq(DNAStr, "Cys"))}\n')  
print(f'[12] : Amino Acid Frequency (Cys): {countCodonFrequency(DNAStr, "Cys")}\n')
print(f'[13] : Tyr Count: {len(aminoacidFreq(DNAStr, "Tyr"))}\n')
print(f'[14] : Amino Acid Frequency (Tyr): {countCodonFrequency(DNAStr, "Tyr")}\n')
print('[15] : Reading Frames:')    
for frame in reading_frame(DNAStr):
    print(frame)
print('[16] : All Open Reading Frames:')    
for frame in reading_frame(DNAStr):
    print(ORF_protein(frame))
