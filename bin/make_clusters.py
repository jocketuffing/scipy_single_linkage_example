import pandas as pd
from scipy.cluster.hierarchy import single, fcluster,cut_tree
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import dendrogram, linkage,to_tree
from scipy.spatial.distance import squareform
import numpy as np
import json

#Specify the number of ADs for cut off when extracting single linkage clusters
ad_cutoff = 3
print('\nAD cut off set to: '+str(ad_cutoff))

#read the recordids of isolates of interest
isolate_annotation_external = {}
with open("data/annotation.json", "r") as infile:
    isolate_annotation_external = json.load(infile)

#transform keys to integers
isolate_annotation = {}
for id,recordid in isolate_annotation_external.items():
    isolate_annotation[int(id)] = recordid

#read the distance matrix of members in the cluster 2023-10.SALM.09.STRATHCONA
df = pd.read_excel('data/VisualisationDistanceMatrix__11_25_2024 18_33_37.xlsx',index_col='Isolates')

#print the AD distances compared to the german reference strain
print('Distances:\n')
for id,recordid in isolate_annotation.items():
    isolate_annotation[int(id)] = recordid
    if id == 7505:
        continue
    print(isolate_annotation[7505]+' <=> '+recordid+': '+str(df.loc[7505,id])+' ADs')
print('\n')

#create a dataframe to store cluster information for the isolates
df_isolates = pd.DataFrame(columns=['Isolate_id','cluster_id','annotation'])
df_isolates.set_index('Isolate_id',inplace=True)

#transform the distance matrix dataframe to a condensed distance matrix
condensed_dist_matrix = squareform(df)

#make a single linkage from the matrix
Z = linkage(condensed_dist_matrix,method='single',optimal_ordering=False) #'single'

#extract single linkage clusters using cut_tree with specified AD cut off
cluster_ids = cut_tree(Z, height=ad_cutoff)

#collection of results
unique_ids, counts = np.unique(cluster_ids, return_counts=True)
count_dict = dict(zip(unique_ids, counts))

#iterate over the elements in the distance matrix
for i in range(len(df.index)):
    
    #retrieve the isolate id
    isolate_id = df.index[i]

    #retrieve the cluster id and number of memebers created from the cut tree process
    cluster_id = cluster_ids[i][0]
    count = count_dict[cluster_id]

    #retrieve recordid if available
    annotation = ''

    if isolate_id in isolate_annotation.keys():
        annotation = isolate_annotation[isolate_id]

    if int(count) == 1:
        cluster_id = 'singleton'
    else:
        cluster_id = 'cluster_'+str(cluster_id)

    #save the results in the dataframe
    df_isolates.loc[isolate_id,'cluster_id'] = cluster_id
    df_isolates.loc[isolate_id,'annotation'] = annotation

#print the cluster result of isolates of interest
print('Single linkage clusters:\n')
print(df_isolates.loc[df_isolates['annotation']!=''])
print('\n')