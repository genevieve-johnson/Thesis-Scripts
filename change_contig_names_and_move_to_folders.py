import os
from os import walk

#VirSorter outputs predicted phages with confidences of 1-3 and 4-6,
#1 & 4 being most confident and 3 & 6 being least confident
#this script goes through all of the VirSorter predicted phage files,
#it then pulls out the 1 & 4 phages and rewrites them to new files in a separate folder
#it also changes the contig names to include the phage name (the name of the file)
#it also repeats this process for 2 & 5 phages

path = "C://Users//genev//Documents//PUTONTI_LAB//fake_virsorter//"
new_path = "C://Users//genev//Documents//PUTONTI_LAB//fake_virsorter//predicted_fasta//"
for (directpath, directnames, files) in walk(new_path):
    for x in files:
        #pull out 1s and 4s
        if x.endswith("1.fasta") or x.endswith("4.fasta"):
            #name of file (phage name)
            fname = str(x[0:-6]+"_")
            #opens that phage file
            with open(new_path+x,'r') as f:
                #creates new file to write out to in new folder
                fout = open(path+"//renamed_seqs_1_4//"+fname+"new.fasta","w")
                #changes contig name to include phage name
                for line in f:
                    if line[0] == ">":
                        fout.write(">" + fname + line[1:])
                    else:
                        fout.write(line)
            fout.close()
        #and then separately pull out 2s and 5s
        elif x.endswith("2.fasta") or x.endswith("5.fasta"):
            #name of file (phage name)
            fname = str(x[0:-6]+"_")
            #opens that phage file
            with open(new_path+x,'r') as f:
                #creates new file to write out to in new folder
                fout = open(path+"//renamed_seqs_2_5//"+fname+"new.fasta","w")
                #changes contig name to include phage name
                for line in f:
                    if line[0] == ">":
                        fout.write(">" + fname + line[1:])
                    else:
                        fout.write(line)
            fout.close()

                    
                    
