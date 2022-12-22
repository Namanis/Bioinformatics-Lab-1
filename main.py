from Bio import SeqIO
from Bio.Seq import Seq
from collections import Counter
from codons import DNA_Codons

for seq_record in SeqIO.parse(r"c:\Users\laury\Desktop\viruses\data\mamalian4.fasta", "fasta"):
    print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record))
    new_seq = seq_record.seq

def translate_seq(self, init_pos=0):
        """Translates a DNA sequence into an aminoacid sequence"""

        return [self[pos:pos + 3] for pos in range(init_pos, len(self) - 2, 3)]


count_nucleotides = {
    'ATG': new_seq.count('ATG'),
    'TAA': new_seq.count('TAA'),
    'TAG': new_seq.count('TAG'),
    'TGA': new_seq.count('TGA')
}
print(count_nucleotides)
print(new_seq.count("ATG"))

#print(translate_seq(new_seq))

def gen_reading_frames(self):
        """Generate the six reading frames of a DNA sequence, including reverse complement"""
        frames = []
        frames.append(translate_seq(new_seq, 0))
        frames.append(translate_seq(new_seq, 1))
        frames.append(translate_seq(new_seq, 2))
        tmp_seq = self.reverse_complement()
        frames.append(translate_seq(tmp_seq, 0))
        frames.append(translate_seq(tmp_seq, 1))
        frames.append(translate_seq(tmp_seq, 2))
        del tmp_seq
        return frames
#strrep = str(new_seq.seq)

#print(new_seq.reverse_complement())
def proteins_from_rf(aa_seq):
        """Compute all possible proteins in an aminoacid seq and return a list of possible proteins"""
        current_prot = []
        proteins = []
        for aa in aa_seq:
            if aa == "TAG" or aa == "TAA" or aa == "TGA":
                # STOP accumulating amino acids if _ - STOP was found
                if current_prot:
                    for p in current_prot:
                        proteins.append(p)
                    current_prot = []
          
            else:
                # START accumulating amino acids if M - START was found
                if aa == "ATG":
                    current_prot.append("")
                for i in range(len(current_prot)):
                    current_prot[i] += aa
        return proteins




array=[]
new_array=[]
# for rf in gen_reading_frames(new_seq):
#     array=rf
    #print(array)

for array in gen_reading_frames(new_seq):
    new_array = proteins_from_rf(array)
    #print(new_array)

# for array in proteins_from_rf(new_seq):
#     new_array = proteins_from_rf(new_seq)
#     print(new_array)

#print(new_array[0])
#print(new_array[1])
#print(new_array[2])
#print(proteins_from_rf(new_seq))

#print(gen_reading_frames(new_seq))



import re

pattern = re.compile(r'(?=(ATG(?:...)*?)(?=TAG|TGA|TAA))')

Seq= """CATGATGCTCAGCGAGGACAGCAAGGGCCCATTTACAGGAGCATAGTAA"""

revcompseq = Seq[::-1].maketrans("ATGC", "TACG") #reverse complement

print(Seq)
print (pattern.findall(Seq)) #forward search
print (pattern.findall(Seq[::-1].translate(revcompseq))) #backward search

#print(type(new_seq)) #this is the seq attribute of the seq1 object
fowardSeq = str(new_seq) #so cast the seq attribute to str

reverseSeq = str(new_seq.reverse_complement())

#print (pattern.findall(fowardSeq)) #forward search
#print (pattern.findall(reverseSeq))

longFowardSeq = []
longReverseSeq = []

for i in range(len(pattern.findall(fowardSeq))):
    if len(pattern.findall(fowardSeq)[i])>= 100:
        longFowardSeq.append((pattern.findall(fowardSeq)[i]))

        
for i in range(len(pattern.findall(reverseSeq))):
    if len(pattern.findall(reverseSeq)[i])>= 100:
        longReverseSeq.append((pattern.findall(reverseSeq)[i]))    
    

def codon_usage(self, aminoacid):
        """Provides the frequency of each codon encoding a given aminoacid in a DNA sequence"""
        tmpList = []

        for i in range(0, len(self) - 2, 3):
            if DNA_Codons[self[i:i + 3]] == aminoacid:
                tmpList.append(self[i:i + 3])


        freqDict = dict(Counter(tmpList))
        totalWight = sum(freqDict.values())
        for seq in freqDict:
            freqDict[seq] = round(freqDict[seq] / totalWight, 2)
        return freqDict


for i in DNA_Codons:
    currentCodon = 0
    if (currentCodon != str(DNA_Codons[i])):
        print(codon_usage(new_seq, DNA_Codons[i]))
        currentCodon = DNA_Codons[i]
        

    
    # print(currentCodon)
    # print(DNA_Codons)















