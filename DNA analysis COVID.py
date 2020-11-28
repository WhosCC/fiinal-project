from DNAToolkit import *
import matplotlib.pyplot as plt

fasta={}

with open('sequence.fasta') as f:
      sequence = ""  
      for line in f:
          if line.startswith(">"):
            name=line[1:].rstrip()
            fasta[name]=""
            continue
          fasta[name]+=line.strip().upper()
          
DNAStr=fasta['COVID']

#print(fasta['COVID'])

#Plot nucleotide frequency
base=countNucFrequency(DNAStr)
X=["A","T","C","G"]
Y=[int(base["A"]),int(base["T"]),int(base["C"]),int(base["G"])]
plt.bar(x=X, height=Y, color='steelblue')
for x, y in enumerate(Y):
    plt.text(x, y + 100, '%s' % y, ha='center', va='top')
plt.title("Base Frequency")
plt.xlabel("Nucleotide Type")
plt.ylabel("Nucleotide Number")
plt.show()


#print(transcription(DNAStr))

#rDNAStr = DNAcomplement(DNAStr)
#print(rDNAStr)

#print(f'\nSequence: {coloring(DNAStr)}\n')
#print(f'[1] : Sequence Length: {len(DNAStr)}\n')
#print(f'[2] : Nucleotide Frequencies: {countNucFrequency(DNAStr)}\n')
#print(f'[3] : DNA Single Strand DNA Complement:{DNAcomplement(DNAStr)}\n')

#print(f"[4] : Original DNA Strand + DNA Complemnt:\n5' {DNAStr} 3'")
#print(f"   {''.join(['|' for a in range(len(DNAStr))])}")
#print(f"3' {DNAcomplement(DNAStr)} 5'\n")
#TOO LONG TO PLOT?

#print(f'[5] : Transcription Result: {transcription(DNAStr)}\n')
#重复了吗？


#print(f"[6] : Original DNA Strand + Transcription Result:\n5' {DNAStr} 3'")
#print(f"   {''.join(['|' for a in range(len(DNAStr))])}")
#print(f"3' {transcription(rDNAStr)} 5'\n")
#print(f'[7] : GC_Portions of DNA Sequence: {gc_portion(DNAStr)}%\n')
#print(f'[8] : Translation Result: {translation(DNAStr, 0)}\n')

print(f'[9] : Lys Count: {len(aminoacidFreq(DNAStr, "Lys"))}\n')
print(f'[10] : Amino Acid Frequency (Lys): {countCodonFrequency(DNAStr, "Lys")}\n')
print(f'[11] : Cys Count: {len(aminoacidFreq(DNAStr, "Cys"))}\n')  
print(f'[12] : Amino Acid Frequency (Cys): {countCodonFrequency(DNAStr, "Cys")}\n')
print(f'[13] : Tyr Count: {len(aminoacidFreq(DNAStr, "Tyr"))}\n')
print(f'[14] : Amino Acid Frequency (Tyr): {countCodonFrequency(DNAStr, "Tyr")}\n')

lys=len(aminoacidFreq(DNAStr, "Lys"))
cys=len(aminoacidFreq(DNAStr, "Cys"))
tyr=len(aminoacidFreq(DNAStr, "Tyr"))

X=["Lys","Cys","Tyr"]
Y=[lys,cys,tyr]
plt.bar(x=X, height=Y, color='steelblue')
for x, y in enumerate(Y):
    plt.text(x, y +15, '%s' % y, ha='center', va='top')
plt.title("Base Frequency")
plt.xlabel("Nucleotide Type")
plt.ylabel("Nucleotide Number")
plt.show()

print('[15] : Reading Frames:')    
for frame in reading_frame(DNAStr):
   print(frame)
print('[16] : All Open Reading Frames:')    
for frame in reading_frame(DNAStr):
    print(ORF_protein(frame))