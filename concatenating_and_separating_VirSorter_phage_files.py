import os
from  os import walk
from Bio import SeqIO

#my VirSorter predicted phage were output in batches as there were too many files to do all at once
#so in order to have all the predicted phage in one file, the large directory needed to be
#walked through and the batch files' contents concatenated together

path = "C://Users//genev//Documents//PUTONTI_LAB//PATRIC_VirSorter_results//all_batch_contigs//"

#this creates a list of the batch files within the larger directory
file_list = []
for (dirpath, dirnames, filenames) in walk(path):
    for y in filenames:
        if y.endswith(".contigs.fasta"):
            file_list.append(y)

#the fasta contents of each batch file are then written out to a new file with all of the predicted phages
total_phage_contigs = open("total_pa_phage_contigs.fasta","w")
for i in file_list:
    file = open(path+i,"r").read()
    total_phage_contigs.write(file)

total_phage_contigs.close()

#for another part of my workflow (Anvi'o) all of the predicted phages needed to be separated into their own files
#so after opening the large concatenated file of phages, each phage has its own file created
#and the fasta ID and sequence for that phage is added to the file
file_sequences = SeqIO.parse(open(path+"total_pa_phage_contigs.fasta"),"fasta")
file_sequences = list(file_sequences)

for i in range(len(file_sequences)):
    new_file = open(path + "all_separate_contigs\\" + str(i) + "_" + file_sequences[i].description[0:19] + ".fasta","w")
    new_file.write(">" + str(file_sequences[i].description))
    new_file.write(str(file_sequences[i].seq))
    new_file.close()
    
