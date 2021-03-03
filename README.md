# Thesis-Scripts
Python scripts for thesis - prediction, annotation, analysis on bacteriophage genomes

Packages used: os, walking, numpy, Biopython, SeqIO, pandas

I have written several scripts to assist in sorting through and analyzing my sequence data for my thesis project. The main point of this project is to look at bacteriophage diversity and to compare a few specific phage genes.

I have uploaded some of the larger scripts written for this project including a script to concatenate or separate contigs into one large file or each into their own file. (concatenating_and_separating_VirSorter_phage_files.py)

I also have a script that finds the files of VirSorter predicted phage with higher confidence (category 1 and category 4), then it renames the contigs to include the name of the phages, and finally it moves the phages into a new folder to organize the phage by categories. (change_contig_names_and_move_to_folders.py)

There is a script that is used after annotating genes from PATRIC that goes through the annotation files, pulls out specific genes, and writes out the genes and their sequences to files to organize the genes. (integrase_crepressor_terminase_PATRIC_annotation_results.py)

Another script takes the BLAST CSV results of these specific annotated genes compared to all of the other genes and finds highly similar genes (percent identity and query coverage both above 70%) to add to the specific annotated gene files. (BLAST_result_highly_similar_genes_to_annotations.py)

One script takes in the cluster files for several different genes output from USEARCH and compares the number of phage shared in between different genes' clusters. So for gene-type-1 and gene-type-2, this script goes through all of the clusters of gene-type-1 and counts the number of phage in each cluster that are shared with each cluster of gene-type-2. These shared phage counts are then put into a matrix to find if there any clusters of the two gene types that share all of the same phages. (compare_integrase_crepressor_terminase_USEARCH_clusters.py)

Finally, the last script takes in an edgelist that has been produced by Anvi'o that consists of the number of genes shared between each of the phages. The script then goes through the edgelist to remove any self-loops or duplicate edges in order to reduce the size of the file and make it easier for downstream analysis in Cytoscape. (anvio_edgelist_no_duplicates_no_selfloops.py)
