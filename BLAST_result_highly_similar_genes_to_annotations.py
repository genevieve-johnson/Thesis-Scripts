import pandas
from Bio import SeqIO

#we already have files for the annotated integrases, c repressors, terminases
#but what about genes with high BLAST similarity to these 3 genes that weren't annotated correctly?
#first a local BLAST compared all of the genes to each of the 3 genes to find similarities
#this script takes the BLAST result CSV files and finds the genes that have a percent identity above 70% and query coverage above 70%
#then it adds these highly similar genes to the integrase, c repressor, or terminase combined files

path = "C:\\Users\\genev\\Documents\\PUTONTI_LAB\\PATRIC_VirSorter_results\\"

#combined files to write out the old annotated and new highly similar genes
full_integrase_data = open(path + "combined_full_integrases.fasta","w")
full_c_repressor_data = open(path + "combined_full_c_repressors.fasta","w")
full_terminase_data = open(path + "combined_full_terminases.fasta","w")

#these are the known annotated integrase, c repressor, and terminase genes
known_integrases = open(path + "all_integrase_fastas.fasta","r")
known_c_repressors = open(path + "all_c_repressor_fastas.fasta","r")
known_terminases = open(path + "all_terminase_fastas.fasta","r")

#write all of the known annotated integrase, c repressor, and terminase genes to the new combined files
for line in known_integrases:
    full_integrase_data.write(line)
known_integrases.close()

for line in known_c_repressors:
    full_c_repressor_data.write(line)
known_c_repressors.close()

for line in known_terminases:
    full_terminase_data.write(line)
known_terminases.close()

#load the BLAST CSV files
integrase_data = path + "integrase_blast.csv"
integrase_df = pandas.read_csv(integrase_data)

c_repressor_data = path + "c_repressor_blast.csv"
c_repressor_df = pandas.read_csv(c_repressor_data)

terminase_data = path + "terminase_blast.csv"
terminase_df = pandas.read_csv(terminase_data)

#go through BLAST CSV files and if both the percent identity & query coverage are about 70%,
#add those highly similar genes to the corresponding list
integrase_pident_qcovs_70_up = []
for i in range(len(integrase_df.pident)):
    if integrase_df.pident[i] >= 70:
        if integrase_df.qcovs[i] >= 70:
            integrase_pident_qcovs_70_up.append(integrase_df.sseqid[i])

c_repressor_pident_qcovs_70_up = []
for i in range(len(c_repressor_df.pident)):
    if c_repressor_df.pident[i] >= 70:
        if c_repressor_df.qcovs[i] >= 70:
            c_repressor_pident_qcovs_70_up.append(c_repressor_df.sseqid[i])

terminase_pident_qcovs_70_up = []
for i in range(len(terminase_df.pident)):
    if terminase_df.pident[i] >= 70:
        if terminase_df.qcovs[i] >= 70:
            terminase_pident_qcovs_70_up.append(terminase_df.sseqid[i])

#this section just finds sequences for all of the highly similar genes from the BLAST results
#these genes are then added to the corresponding combined file with the known annotated genes
everything_besides_integrase = SeqIO.parse(open(path+"everything_besides_integrases_fastas.fasta"),"fasta")
for i in everything_besides_integrase:
    if i.description in integrase_pident_qcovs_70_up:
        full_integrase_data.write(">"+i.description+"\n")
        full_integrase_data.write(str(i.seq)+"\n")

everything_besides_c_repressor = SeqIO.parse(open(path+"everything_besides_c_repressors_fastas.fasta"),"fasta")
for i in everything_besides_c_repressor:
    if i.description in c_repressor_pident_qcovs_70_up:
        full_c_repressor_data.write(">"+i.description+"\n")
        full_c_repressor_data.write(str(i.seq)+"\n")

everything_besides_terminase = SeqIO.parse(open(path+"everything_besides_terminases_fastas.fasta"),"fasta")
for i in everything_besides_terminase:
    if i.description in terminase_pident_qcovs_70_up:
        full_terminase_data.write(">"+i.description+"\n")
        full_terminase_data.write(str(i.seq)+"\n")
        

full_integrase_data.close()
full_c_repressor_data.close()
full_terminase_data.close()
