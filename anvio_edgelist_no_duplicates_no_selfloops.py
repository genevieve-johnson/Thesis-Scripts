import pandas
import numpy as np

#the edgelist of all genes shared between every single phage is output from Anvi'o
#this edgelist does contain some self-loops and duplicate edges
#since there are over 6000 nodes and hundreds of thousands of edges, we need to cut down the size of the file
#so this script gets rid of any self-loops and duplicate edges that are unnecessary to include in the edgelist

path = "C:\\Users\\genev\\Documents\\PUTONTI_LAB\\"
cyto_data = path + "anvio_edgelist_for_cytoscape.csv"
cyto_df = pandas.read_csv(cyto_data,header=None)
cyto_array = np.array(cyto_df)

#each row of the array is [the name of first phage, the name of second phage, number genes shared]
#creates dictionary of key = phage names as tuple and value = number genes shared
cyto_dict = {}
for i in cyto_array:
    tuple_key = (i[0],i[1])
    cyto_dict[tuple_key] = i[2]

#if the first phage name is the same as the second phage name, this is a self-loop and the number of genes is set to 0
#if the two phages are already in the dictionary, we remove the duplicates by setting the number genes to 0
for i in cyto_dict:
    if i[0] == i[1]:
        cyto_dict[i] = 0
    reverse = (i[1],i[0])
    if reverse in cyto_dict.keys():
        if cyto_dict[reverse] != 0 and cyto_dict[i] != 0:
            cyto_dict[reverse] = 0

#this is just reshaping the array into a pandas dataframe to make things easier
new_cyto_df = pandas.DataFrame.from_dict(cyto_dict,orient='index',columns=["value"])
new_cyto_df.index.name = "names"
new_cyto_df.reset_index(inplace=True)

#if any of the number of genes are equal to 0 (like the self-loops or duplicates) they are removed from the dataframe
new_cyto_df = new_cyto_df[~(new_cyto_df["value"]==0)]

#this is just making the dataframe look nicer with fixed column names
new_cyto_df['first_name'], new_cyto_df['second_name'] = new_cyto_df.names.str
del new_cyto_df["names"]
new_cyto_df = new_cyto_df[["first_name","second_name","value"]]

#and the dataframe is written out to a new CSV file
new_cyto_df.to_csv(r'C:\\Users\\genev\\Documents\\PUTONTI_LAB\\new_edgelist_nozeros_noduplicates.csv',header=False,index=False)
