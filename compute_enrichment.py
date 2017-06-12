# -*- coding: utf-8 -*-
"""
Created on Sat Jun 10 12:04:38 2017

@author: mtinti-x
"""

import pandas as pd
from string import strip
from goatools.base import download_go_basic_obo
from goatools.associations import read_associations
from goatools.go_enrichment import GOEnrichmentStudy
from goatools.obo_parser import GODag
import scipy.cluster.hierarchy as sch
import numpy as np
import os

#this function parse the trypdb tab file
#to create the input for the read_associations
#format = gene_id \t go_term \n
def format_file(trytrip_file='', res_file = ''):
    res_file = open(res_file,'w')
    for l in open(trytrip_file):
        item_list = l.split('\t')
        gene_id = item_list[0]
        for item in item_list:
            if 'GO:' in item:
                temp_goes = [strip(n) for  n in item.split(',')]
                for go_id in temp_goes:
                    res_file.write(gene_id+'\t'+go_id+'\n')
    res_file.close()

def prepare_data():
    #download the last obo
    #obo_fname = download_go_basic_obo()
    #format the last trypdb search res file
    #for goterms
    download_go_basic_obo()
    trytrip_file='data/GenesByGoTerm_Summary.txt'
    associations_file = 'data/associations.txt'
    format_file(trytrip_file=trytrip_file, res_file = associations_file)    

        

    
def make_clustering_df(in_df = pd.DataFrame(), 
                       min_number_of_clusters=4,
                       res_df_name = '',
                       method='', metric=''):
    
    link = sch.linkage(in_df, method, metric)
    
    den = sch.dendrogram(link, 
                         color_threshold=20, 
                         orientation='left')
    #plt.show()
    print'done den'
    cut_distances = np.arange(0,1000,0.01)
    done = set()
    for cut_distance in cut_distances:
        clusters = sch.fcluster(link, cut_distance, criterion='distance')
        in_df['clusters_'+str(cut_distance)]=clusters
        n_of_clusters =  in_df['clusters_'+str(cut_distance)].value_counts().shape[0]
        
        if n_of_clusters not in done:
            done.add(n_of_clusters)
            print cut_distance, n_of_clusters
        else:
            del in_df['clusters_'+str(cut_distance)]
        if n_of_clusters <=min_number_of_clusters:
            break
        
    res_df = in_df[[n for n in in_df.columns if 'clusters_' in n]]
    res_df.to_csv(res_df_name)
    print 'done'
    print res_df.shape    
    

          
def score(in_genes, index_res, cut_distance, cluster_id):
    def count(goea_results_sig):
            res = {'BP' : 0,
                   'CC' : 0,
                   'MF' : 0}
            for item in goea_results_sig:
                item_list = str(item).split()            
                go_type = item_list[1] 
                res[go_type]+=1
            return res   
    goea_results_all = goeaobj.run_study(in_genes)
    goea_results_sig = [r for r in goea_results_all if r.p_fdr_bh < 0.05]
    res = count(goea_results_sig)
    if sum(res.values())>0:
        res_file = open(res_path+'/'+'res_'+str(cut_distance)+'_'+str(cluster_id)+'.txt','w')
        res_file.write('\t'.join(['cut_distance','cluster_id','CC','MF','BP'])+'\n')
        res_file.write('\t'.join( [str(n) for n in [cut_distance,cluster_id,res['CC'],res['MF'],res['BP']]] ))
        res_file.close()   
    


            
if __name__ == '__main__':
    print 'start'
    
#==============================================================================
#     prepare_data()
#     norMax = lambda x: x / x.max()
#     in_df = pd.DataFrame.from_csv('data/in_df.txt',sep='\t')
#     #drop all zeros
#     in_df = in_df[(in_df.T != 0).any()] 
#     in_df = in_df.dropna()
#     in_df = in_df.apply(norMax,1)
#     in_df[in_df<0.01]=0        
#     make_clustering_df(in_df = in_df, 
#                        min_number_of_clusters=4,
#                        res_df_name = 'res_clustering.csv',
#                        method='ward' ,
#                        metric='euclidean'
#                           )    
#==============================================================================



    data_path = os.path.abspath("data")
    res_path = os.path.abspath("results")    
    in_go=data_path+'/go-basic.obo'
    in_assoc=data_path+'/associations.txt'
    gene_found = [strip(n.split('\t')[0]) for n in open(data_path+'/in_df.txt').read().split('\n')]
    obodag = GODag(in_go)
    geneid2gos = read_associations(in_assoc)
    goeaobj = GOEnrichmentStudy(
                gene_found, # List of mouse protein-coding genes
                geneid2gos, # geneid/GO associations
                obodag, # Ontologies
                propagate_counts = False,
                alpha = 0.05, # default significance cut-off
                methods = ['fdr_bh'])      
    
    clustering_df = pd.DataFrame.from_csv('res_clustering.csv')
    for index, col in enumerate(clustering_df.columns):
        in_data = clustering_df[col].value_counts()
        in_data = in_data[in_data>=2]
        cut_distance = col.split('_')[-1]
        print index, col, cut_distance
        for cluster in in_data.index.values:
            genes = clustering_df[col][clustering_df[col]==cluster]
            genes = genes.index.values
            #print cluster, len(genes)
            #print 'send', index, cut_distance, cluster, 'to controller'
            score(genes, index, cut_distance, cluster)    
    

    
    
    
    #now to the notebook!
     
    
    
    
    

    
    
                             
        
    

