from Bio import SeqIO
import numpy

#USEARCH clusters of c repressor sequences, clusters of integrase sequences, and clusters of terminase sequences
#we want to see how many phages each integrase cluster shares with each terminase cluster
#and how many phages each c repressor cluster shares with each terminase cluster

#create matrix of terminase clusters on x axis and integrase clusters on y axis
#create matrix of terminase clusters on x axis and c repressor clusters on y axis

#c repressor (c) has clusters 0-43
#integrase (i) has clusters 0-71
#terminase (t) has clusters 0-80

path = "C://Users//genev//Documents//PUTONTI_LAB//PATRIC_VirSorter_results//clusters//"

c_repressor_list = []
integrase_list = []
terminase_list = []

#make lists of cluster file names for easier access to files
for i in range(0,44):
    c_repressor_list.append("c_repressor_cluster"+str(i))

for i in range(0,72):
    integrase_list.append("integrase_cluster"+str(i))

for i in range(0,81):
    terminase_list.append("terminase_cluster"+str(i))


#make dictionaries for each of c/i/t with 0-43/0-71/0-80 as key and the phages in those clusters as the values
c_repressor_dict = {}
for i in range(len(c_repressor_list)):
    file = SeqIO.parse(open(path+c_repressor_list[i]),"fasta")
    c_repressor_dict[i] = []
    for j in file:
        #fix contig names so they just include phage name, not position of gene
        name = j.description[0:(j.description.find("cat")+27)]
        c_repressor_dict[i].append(name)

integrase_dict = {}
for i in range(len(integrase_list)):
    file = SeqIO.parse(open(path+integrase_list[i]),"fasta")
    integrase_dict[i] = []
    for j in file:
        #fix contig names so they just include phage name, not position of gene
        name = j.description[0:(j.description.find("cat")+27)]
        integrase_dict[i].append(name)

terminase_dict = {}
for i in range(len(terminase_list)):
    file = SeqIO.parse(open(path+terminase_list[i]),"fasta")
    terminase_dict[i] = []
    for j in file:
        #fix contig names so they just include phage name, not position of gene
        name = j.description[0:(j.description.find("cat")+27)]
        terminase_dict[i].append(name)


#matrix with c_repressor on y axis and terminase on x axis
#plus one extra row and column to get total number phages in each cluster for each
c_t_matrix = numpy.zeros((45,82))

#c_repressor and terminase - count number shared phages between each cluster
#i = c = rows and j = t = columns
#matrix[c][t]
for i in range(0,44):
    for j in range(0,81):
        count = 0
        for x in c_repressor_dict[i]:
            if x in terminase_dict[j]:
                count += 1
        c_t_matrix[i][j] = count
        c_t_matrix[i][81] = len(c_repressor_dict[i])
        c_t_matrix[44][j] = len(terminase_dict[j])


#matrix with integrase on y axis and terminase on x axis
#plus one extra row and column to get total number phages in each cluster for each
i_t_matrix = numpy.zeros((73,82))

#integrase and terminase - count number shared phages between each cluster
#i = i = rows and j = t = columns
#matrix[i][t]
for i in range(0,72):
    for j in range(0,81):
        count = 0
        for x in integrase_dict[i]:
            if x in terminase_dict[j]:
                count += 1
        i_t_matrix[i][j] = count
        i_t_matrix[i][81] = len(integrase_dict[i])
        i_t_matrix[72][j] = len(terminase_dict[j])


#now go through both matrices and if the number in the column/row equals the last row/column numbers
#then they share all the phage that is in one of two the clusters

c_t_cluster_tuples = []
for row in range(0,43):
    for column in range(0,80):
        number = c_t_matrix[row][column]
        if number > 0:
            if number == c_t_matrix[row][81]:
                c_t_cluster_tuples.append((row,column,"sharing all "+str(int(number))+" prophages in c_repressor cluster"))
            if number == c_t_matrix[44][column]:
                c_t_cluster_tuples.append((row,column,"sharing all "+str(int(number))+" prophages in terminase cluster"))

                
i_t_cluster_tuples = []
for row in range(0,71):
    for column in range(0,80):
        number = i_t_matrix[row][column]
        if number > 0:
            if number == i_t_matrix[row][81]:
                i_t_cluster_tuples.append((row,column,"sharing all "+str(int(number))+" prophages in integrase cluster"))
            if number == i_t_matrix[72][column]:
                i_t_cluster_tuples.append((row,column,"sharing all "+str(int(number))+" prophages in terminase cluster"))


#write this information out into a file
output = open(path+"cluster_sharing.txt","w")

output.write("C repressor and Terminase clusters with same prophages:"+"\n")
for i in c_t_cluster_tuples:
    output.write("c_repressor_cluster"+str(i[0])+", terminase_cluster"+str(i[1])+", "+str(i[2])+"\n")
output.write("\n")
output.write("Integrase and Terminase clusters with same prophages:"+"\n")
for i in i_t_cluster_tuples:
    output.write("integrase_cluster"+str(i[0])+", terminase_cluster"+str(i[1])+", "+str(i[2])+"\n")

output.close()
