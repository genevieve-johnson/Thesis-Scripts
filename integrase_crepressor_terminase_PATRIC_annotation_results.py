import os
from Bio import SeqIO

#List of all the PATRIC annotation result batch names -- these are the directory names
#List of all the PATRIC annotation file names
batches = ["batch_1A","batch_1B","batch_1C","batch_2","batch_3","batch_4","batch_5","batch_6","batch_7","batch_8","batch_9","batch_10","batch_11","batch_12","batch_13","batch_14","batch_15","batch_cp1","batch_cp2","batch_cp3","batch_101","batch_102","batch_103","batch_104","batch_105","batch_106","batch_107"]
file_names = ["Batch 1A","Batch 1B","Batch 1C","Batch 2","Batch 3","Batch 4","Batch 5","Batch 6","Batch 7","Batch 8","Batch 9","Batch 10","Batch 11","Batch 12","Batch 13","Batch 14","Batch 15","Batch cp1","Batch cp2","Batch cp3","Batch 101","Batch 102","Batch 103","Batch 104","Batch 105","Batch 106","Batch 107"]

#all the files to write out to
absolutely_everything = open("absolutely_all_genes.fasta","w")
all_integrase_fastas = open("all_integrase_fastas.fasta","w")
everything_besides_integrases_fastas = open("everything_besides_integrases_fastas.fasta","w")
all_c_repressor_fastas = open("all_c_repressor_fastas.fasta","w")
everything_besides_c_repressors_fastas = open("everything_besides_c_repressors_fastas.fasta","w")
all_terminase_fastas = open("all_terminase_fastas.fasta","w")
everything_besides_terminases_fastas = open("everything_besides_terminases_fastas.fasta","w")

#the final dictionaries for each gene
integrase_contigs = {}
other_than_integrase_contigs = {}
c_repressor_contigs = {}
other_than_c_repressor_contigs = {}
terminase_contigs = {}
other_than_terminase_contigs = {}

#going through all the batch files to find the integrases, c repressors, and terminases in the annotations
for i in range(len(batches)):
    batch = batches[i]
    file_name = file_names[i]
    path = "C://Users//genev//Documents//PUTONTI_LAB//PATRIC_VirSorter_results//"+batch+"//"

    fastas = SeqIO.parse(open(path+"Pseudomonas aeruginosa "+file_name+".feature_protein.fasta"),"fasta")

    #temporary dictionaries for the different genes 
    integrases = {}
    other_than_integrase = {}
    c_repressors = {}
    other_than_c_repressor = {}
    terminases = {}
    other_than_terminase = {}

    #if integrase is in the fasta description, split the description to get the feature ID name
    #feature ID name is key and the sequence is the value
    #otherwise put the other genes in the everything-except-integrase dictionary
    for fasta in fastas:
        if "integrase" in fasta.description:
            split_up = fasta.description.split(" ")
            integrases[split_up[0]] = str(fasta.seq)
        elif "Integrase" in fasta.description:
            split_up = fasta.description.split(" ")
            integrases[split_up[0]] = str(fasta.seq)
        else:
            split_up = fasta.description.split(" ")
            other_than_integrase[split_up[0]] = str(fasta.seq)

    #if c repressor is in the fasta description, split the description to get the feature ID name
    #feature ID name is key and the sequence is the value
    #otherwise put the other genes in the everything-except-crepressor dictionary
    for fasta in fastas:
        if "repressor" in fasta.description and "cI" in fasta.description:
            split_up = fasta.description.split(" ")
            c_repressors[split_up[0]] = str(fasta.seq)
        elif "repressor" in fasta.description and "CI" in fasta.description:
            split_up = fasta.description.split(" ")
            c_repressors[split_up[0]] = str(fasta.seq)
        elif "repressor" in fasta.description and "C1" in fasta.description:
            split_up = fasta.description.split(" ")
            c_repressors[split_up[0]] = str(fasta.seq)
        elif "repressor" in fasta.description and "c1" in fasta.description:
            split_up = fasta.description.split(" ")
            c_repressors[split_up[0]] = str(fasta.seq)
        else:
            split_up = fasta.description.split(" ")
            other_than_c_repressor[split_up[0]] = str(fasta.seq)

    #if terminase is in the fasta description, split the description to get the feature ID name
    #feature ID name is key and the sequence is the value
    #otherwise put the other genes in the everything-except-terminase dictionary
    for fasta in fastas:
        if "terminase" in fasta.description:
            split_up = fasta.description.split(" ")
            terminases[split_up[0]] = str(fasta.seq)
        elif "Terminase" in fasta.description:
            split_up = fasta.description.split(" ")
            terminases[split_up[0]] = str(fasta.seq)
        else:
            split_up = fasta.description.split(" ")
            other_than_terminase[split_up[0]] = str(fasta.seq)


    #this section just takes the feature ID name and goes through the features file to find the name of the phage
    #the phage name is then the key in the final gene dictionary and the sequence from the temporary dictionary is added as the value
    i_features = open(path+"Pseudomonas aeruginosa "+file_name+".features.txt","r")
    for line in i_features:
        i_feature_split = line.split("\t")
        if i_feature_split[0] in integrases.keys():
            integrase_contigs[i_feature_split[1]] = integrases[i_feature_split[0]]
        elif i_feature_split[0] in other_than_integrase.keys():
            other_than_integrase_contigs[i_feature_split[1]] = other_than_integrase[i_feature_split[0]]
    
    c_features = open(path+"Pseudomonas aeruginosa "+file_name+".features.txt","r")
    for line in c_features:
        c_feature_split = line.split("\t")
        if c_feature_split[0] in c_repressors.keys():
            c_repressor_contigs[c_feature_split[1]] = c_repressors[c_feature_split[0]]
        elif c_feature_split[0] in other_than_c_repressor.keys():
            other_than_c_repressor_contigs[c_feature_split[1]] = other_than_c_repressor[c_feature_split[0]]
    
    t_features = open(path+"Pseudomonas aeruginosa "+file_name+".features.txt","r")
    for line in t_features:
        t_feature_split = line.split("\t")
        if t_feature_split[0] in terminases.keys():
            terminase_contigs[t_feature_split[1]] = terminases[t_feature_split[0]]
        elif t_feature_split[0] in other_than_terminase.keys():
            other_than_terminase_contigs[t_feature_split[1]] = other_than_terminase[t_feature_split[0]]


#this section is just writing all of the dictionaries out to their respective files
for key,value in integrase_contigs.items():
    all_integrase_fastas.write(">"+key+"\n")
    all_integrase_fastas.write(value+"\n")
    absolutely_everything.write(">"+key+"\n")
    absolutely_everything.write(value+"\n")
all_integrase_fastas.close()
absolutely_everything.close()

for key,value in other_than_integrase_contigs.items():
    everything_besides_integrases_fastas.write(">"+key+"\n")
    everything_besides_integrases_fastas.write(value+"\n")
    absolutely_everything.write(">"+key+"\n")
    absolutely_everything.write(value+"\n")
everything_besides_integrases_fastas.close()
absolutely_everything.close()

for key,value in c_repressor_contigs.items():
    all_c_repressor_fastas.write(">"+key+"\n")
    all_c_repressor_fastas.write(value+"\n")
all_c_repressor_fastas.close()

for key,value in other_than_c_repressor_contigs.items():
    everything_besides_c_repressors_fastas.write(">"+key+"\n")
    everything_besides_c_repressors_fastas.write(value+"\n")
everything_besides_c_repressors_fastas.close()

for key,value in terminase_contigs.items():
    all_terminase_fastas.write(">"+key+"\n")
    all_terminase_fastas.write(value+"\n")
all_terminase_fastas.close()

for key,value in other_than_terminase_contigs.items():
    everything_besides_terminases_fastas.write(">"+key+"\n")
    everything_besides_terminases_fastas.write(value+"\n")
everything_besides_terminases_fastas.close()

