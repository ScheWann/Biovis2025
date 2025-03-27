from __future__ import division
from scipy.stats import fisher_exact, hypergeom,binom
from random import shuffle
import operator
import random

import math
from sklearn.metrics import mean_squared_error,r2_score
import numpy as np
import pandas as pd
import scanpy as sc
import anndata
from scipy import sparse
from sympy import *
import matplotlib.pyplot as plt

from multiprocessing import Pool
from functools import partial

import datetime  #change 1
import fisher

import networkx as nx
import itertools

def get_DEG_uniq(rdata:'raw annoData',
                   adata:'annoData',
                   group_key:'cell type of annoData',
                   power=11,
                   ratio = 0.2,
                   p_threshold=0.01,
                   q_threshold=0.05):
    
    """
    get the differentially expressed genes of each cell cluster. these genes is uniq in each cell cluster.
    
    Parameters
    ----------
    rdata : raw annoData,we recommoned you use the normalized and log2 transformated raw data;
    adata : annoData, the annoData filtered by HVG is recommended;
    group_key : the cluster of cell,default,'louvain',also 'leiden' or cell_type which defined by user;
    ratio : the possibility whether gene is a  differentially expressed gene. the value is between 0 and 1;
    p_threshold : the threshold of p-value. the value is between 0 and 1;
    q_threshold : the threshold of q_value. the value is between 0 and 1;
    
    Returns
    -------
    a dataFrame of genes which is differentially expressed in each cell type and uniq in each cell type,
    which contains cell type,gene name ,ratio, p-value and q-value.  
    """
    
    #start = datetime.datetime.now()
    #print(start);
    
    
    # get the raw data frame
    power = power;
    print('power: ',power);
    print('get the raw data frame...');
    rdata_df = rdata.to_df();
    num_allCells = len(rdata_df.index);
    
    # struct the cell type sets for enrichment analysis
    print('struct the cell type sets for enrichment analysis...')
    adata_cell_type_df = pd.DataFrame(adata.obs[group_key]);
    cell_type_index = pd.Categorical(adata.obs[group_key]).categories;
    
    cell_sets = dict();
    
    for ct in cell_type_index:
        cell_sets[ct] = list(adata_cell_type_df.loc[adata_cell_type_df[group_key] ==ct,:].index);
    
    iter_genes = rdata_df.columns;
    
    print('Fisher_test_for_each_gene...')
      
    genes_ratio = dict();
    genes_pv = dict();
    genes_qv =dict();
    genes_score =dict();
    genes_means = dict();
    
    num = 0;
    for gene in iter_genes:
        num +=1; 
        gene_ratio, gene_pv, gene_qv,gene_score, gene_means = Fisher_test_for_each_gene(rdata_df,cell_sets,num_allCells,gene,power);
        if gene_ratio[gene] == []:
            continue;
        else:
            genes_ratio = dict(genes_ratio,**gene_ratio);
            genes_pv = dict(genes_pv,**gene_pv);
            genes_qv = dict(genes_qv,**gene_qv);
            genes_score = dict(genes_score,**gene_score);
            genes_means = dict(genes_means,**gene_means);
        if num%1000 ==0:
            print('whole ',num,' genes have been done.')

    print('merge differentially expressed genes...')
    genes_ratio_df = pd.DataFrame(genes_ratio);
    genes_ratio_df.index = cell_type_index;
    
    genes_pv_df = pd.DataFrame(genes_pv);
    genes_pv_df.index = cell_type_index;
    
    genes_qv_df = pd.DataFrame(genes_qv);
    genes_qv_df.index = cell_type_index;
    
    genes_score_df = pd.DataFrame(genes_score);
    genes_score_df.index = cell_type_index;
    
    genes_means_df = pd.DataFrame(genes_means);
    genes_means_df.index = cell_type_index;
    
    ra = ratio;
    pt = p_threshold;
    qt = q_threshold;

    ct_gene_r_p_q_s_list =[];
    
    for gene in genes_ratio_df.columns:
        gene = gene;
        #gene_on_ct = genes_ratio_df.loc[:,gene].idxmax(axis=0);
        gene_on_ct = genes_score_df.loc[:,gene].idxmax(axis=0);
        
        if genes_ratio_df.loc[gene_on_ct,gene]>=ra:
            
            if  genes_pv_df.loc[gene_on_ct,gene]<=pt and genes_qv_df.loc[gene_on_ct,gene]<=qt:
                ct = gene_on_ct;
                g = gene;
                r = genes_ratio_df.loc[gene_on_ct,gene];
                p = genes_pv_df.loc[gene_on_ct,gene];
                q = genes_qv_df.loc[gene_on_ct,gene];
                s = genes_score_df.loc[gene_on_ct,gene];
                m = genes_means_df.loc[gene_on_ct,gene];
                
                ct_gene_r_p_q_s_list.append([ct,g,r,p,q,s,m]);
            else:
                continue;
        else:
            continue;
    markers_s = pd.DataFrame(ct_gene_r_p_q_s_list,columns=['cell_type','gene_name','ratio','p_value','q_value','score','mean_exValue']);
    print('Done!')
    return markers_s

def get_DEG_multi(rdata:'raw annoData',
                     adata:'annoData',
                     group_key:'cell type of annoData',
                     power=11,
                     ratio=0.2,
                     p_threshold=0.01,
                     q_threshold=0.05):
    
   
    """
    Get the differentially expressed genes of each cell cluster. But these genes is not uniq in each cell cluster, 
    which means a gene could be differentially expressed in two or more cell clusters.
    
    Parameters
    ----------
    rdata : raw annoData,we recommoned you use the normalized and log2 transformated raw data;
    adata : annoData, the annoData filtered by HVG is recommended;
    power : int, default:10. a parameter for curve fitting of gene expression pattern.
    group_key : the cluster of cell,default,'louvain',also 'leiden' or cell_type which defined by user;
    power : int, default:10. a parameter for curve fitting of gene expression pattern.
    ratio : the possibility whether gene is a  differentially expressed gene. the value is between 0 and 1;
    p_threshold : the threshold of p-value. the value is between 0 and 1;
    q_threshold : the threshold of q_value. the value is between 0 and 1;
    
    Returns
    -------
    a dataFrame of cell_clusters which contain differentially expressed genes, keys are the cell types and 
    values are arrays which contain gene name ,ratio, p-value and q-value.   
    """
    
    #start = datetime.datetime.now()
    #print(start);
    
    
    
    # get the raw data frame
    power = power;
    print('get the raw data frame...')
    rdata_df = rdata.to_df();
    num_allCells = len(rdata.obs_names)
    
    # struct the cell type sets for enrichment analysis
    print('struct the cell type sets for enrichment analysis...')
    adata_cell_type_df = pd.DataFrame(adata.obs[group_key]);
    cell_type_index = pd.Categorical(adata.obs[group_key]).categories;
    
    cell_sets = dict();
    
    for ct in cell_type_index:
        cell_sets[ct] = list(adata_cell_type_df.loc[adata_cell_type_df[group_key] ==ct,:].index);
    
    iter_genes = rdata_df.columns;
    
    print('Fisher_test_for_each_gene...')
      
    genes_ratio = dict();
    genes_pv = dict();
    genes_qv = dict();
    genes_score = dict();
    genes_means = dict();
    
    num = 0;
    for gene in iter_genes:
        num +=1; 
        gene_ratio, gene_pv, gene_qv, gene_score, gene_means = Fisher_test_for_each_gene(rdata_df,cell_sets,num_allCells,gene,power);
        if gene_ratio[gene] == []:
            continue;
        else:
            genes_ratio = dict(genes_ratio,**gene_ratio);
            genes_pv = dict(genes_pv,**gene_pv);
            genes_qv = dict(genes_qv,**gene_qv);
            genes_score = dict(genes_score,**gene_score);
            genes_means = dict(genes_means,**gene_means);
            
        if num%1000 ==0:
            print('whole ',num,' genes have been done.')

    print('merge differentially expressed genes...')
    genes_ratio_df = pd.DataFrame(genes_ratio);
    genes_ratio_df.index = cell_type_index;
    genes_ratio_df = genes_ratio_df.T;
    
    genes_pv_df = pd.DataFrame(genes_pv);
    genes_pv_df.index = cell_type_index;
    genes_pv_df = genes_pv_df.T;
    
    genes_qv_df = pd.DataFrame(genes_qv);
    genes_qv_df.index = cell_type_index;
    genes_qv_df = genes_qv_df.T;
    
    genes_score_df = pd.DataFrame(genes_score);
    genes_score_df.index = cell_type_index;
    genes_score_df = genes_score_df.T;
    
    genes_means_df = pd.DataFrame(genes_means);
    genes_means_df.index = cell_type_index;
    genes_means_df = genes_means_df.T;
    
    ra = ratio;
    pt = p_threshold;
    qt = q_threshold;
    deg_genes_ratio_p_q_s_list = [];

    for ct in cell_type_index:
        up_ratio_ct_genes_list = list(genes_ratio_df.loc[genes_ratio_df[ct]>=ra,:].index);
        down_pv_ct_genes_list = list(genes_pv_df.loc[genes_pv_df[ct]<=pt,:].index);
        down_qv_ct_genes_list = list(genes_qv_df.loc[genes_qv_df[ct]<=qt,:].index);
        
        ct_deg_genes_list = list(set(up_ratio_ct_genes_list)&set(down_pv_ct_genes_list)&set(down_qv_ct_genes_list));
        
        if ct_deg_genes_list ==[]:
            continue;
            
        else:
            for gene in ct_deg_genes_list:
                ct=ct;
                g = gene;
                r = genes_ratio_df.loc[gene,ct];
                p = genes_pv_df.loc[gene,ct];
                q = genes_qv_df.loc[gene,ct];
                s = genes_score_df.loc[gene,ct];
                m = genes_means_df.loc[gene,ct];
                
                deg_genes_ratio_p_q_s_list.append([ct,g,r,p,q,s,m]);
    
    markers_m = pd.DataFrame(deg_genes_ratio_p_q_s_list,columns=['cell_type','gene_name','ratio','p_value','q_value','score','mean_exValue']);
    print('Done!')
    return markers_m

def get_genes_location_pseudotime(rdata:'annoData of rdata',
                                  adata: 'annoData of adata',
                                  group_key:'cell type of annoData',
                                  gene_matrix:'fileName of markers matrix from results by get_DEG_single or get_DEG_multiple',
                                  obsm:'obsm',
                                 power=11):
    
    """
    Calculate the pseudotemporal location and ordering of genes.
    
    Parameters
    ----------
    rdata : raw annoData,we recommoned you use the normalized and log2 transformated raw data;
    adata : annoData, the annoData filtered by HVG is recommended;
    group_key : the cluster of cell,default,'louvain',also 'leiden' or cell_type which defined by user;
    gene_matrix:'fileName of markers matrix from results by get_DEG_single or get_DEG_multiple';
    obsm:'obsm',the one of 'obsm' attribution of adata;
    
    Returns
    -------
    a dict of cell_clusters which contain differentially expressed genes, keys are the cell types and 
    values are arrays which contain gene name ,ratio, p-value and q-value.   
    """    
    start = datetime.datetime.now()
    #print(start);
    
    #pool = Pool(n_process);
    
    rdata_df = rdata.to_df();
    
    adata_obsm_df = pd.DataFrame(adata.obsm[obsm]);
    adata_obsm_df.index = adata.obs_names;
    adata_obsm_df['dpt_pseudotime'] = adata.obs.dpt_pseudotime;
    
    markers_s_LGPS = adata.uns[gene_matrix]
    markers_s_LGPS_rdata = markers_s_LGPS.loc[markers_s_LGPS['gene_name'].isin(list(rdata.var_names)),]
    #markers_s_LGPS_rdata.index = list(markers_s_LGPS_rdata['gene_name']);
    
    adata_cell_type_df = pd.DataFrame(adata.obs[group_key]);
    cell_type_index = pd.Categorical(adata.obs[group_key]).categories;
    
    cell_sets = dict();
    
    for ct in cell_type_index:
        cell_sets[ct] = list(adata_cell_type_df.loc[adata_cell_type_df[group_key] ==ct,:].index);

    iter_genes = rdata_df.columns;
    gene_pseudotime_locates =list();
    num = 0;
    for g1 in  iter_genes:
        num +=1; 
        g1_pseudotime_locates = calculate_genes_pseudotime_location(rdata_df,
                                                                    cell_sets,
                                                                    markers_s_LGPS_rdata,
                                                                    adata_obsm_df,
                                                                    g1,
                                                                   power);
        if g1_pseudotime_locates == []:
            continue;
        else:
            for i in range(0,len(g1_pseudotime_locates)):
                gene_pseudotime_locates.append(g1_pseudotime_locates[i])
        
        if num%1000 ==0:
            print('whole ',num,' genes have been done.')
    
    #func = partial(calculate_genes_pseudotime_location,rdata_df,cell_sets,markers_s_LGPS_rdata,adata_obsm_df);
    #gene_pseudotime_locate_dicts = pool.map(func,iter_genes);
    #pool.close();
    #pool.join();

    gene_pseudotime_locates_df = pd.DataFrame(gene_pseudotime_locates);
    
    end = datetime.datetime.now()
    #print(end);
    print('Running time : %s Seconds' %(end-start))
    print('Done!');
    return gene_pseudotime_locates_df;

def Fisher_test_for_each_gene(rdata_df: 'dataFrame of raw data',cell_sets:'all cell sets',num_allCells:'number of all cells',gene:'gene name',power=11):
    """
    get ratio, p value and q value of Fisher test for given gene.
    
    Parameters
    ----------
    rdata_df : data frame format of raw annoData;
    cell_sets : a dict. key is the cell tyoe; value is a list of cell name;
    num_allCells: the number of all cells which is used as the background cells;
    gene: a given gene name;
    
    return
    ----------
    gene_ratio : a dict of which key is gene and value is ratio 
    gene_pv : a dict of which key is gene and value is p value from fisher test
    gene_qv : a dict of which key is gene and claue is q value which is FDR for value by using BH method
    """
    power=power;
    gene = gene;
    gene_highly_cells,gene_mean_exvalue = get_highly_cells_for_each_gene(rdata_df,gene,power);
    test_cells = gene_highly_cells[gene];
    gene_mean = gene_mean_exvalue[gene];
    num_test_cells = len(test_cells); 
    
    gene_pv = dict();
    gene_qv = dict();
    gene_ratio = dict();
    gene_score = dict();
    gene_means = dict();
    
    all_ratio = [];
    all_pv = [];
    all_qv = [];
    all_score = [];
    all_means = [];
    
    if num_test_cells == 0:
        gene_ratio[gene] = [];
        gene_pv[gene] = [];
        gene_qv[gene] = [];
        gene_score[gene] = [];
        gene_means[gene] = [];
        
        return (gene_ratio,gene_pv, gene_qv, gene_score,gene_means);
    else:
        for cs in cell_sets.keys():
            cell_set_cs = cell_sets[cs];
            num_cell_set = len(cell_set_cs);
            a = len(set(test_cells)&set(cell_set_cs));
            b = num_test_cells-a;
            c = num_cell_set-a;
            d = num_allCells-num_cell_set;
            ra = a/float(num_test_cells);
            pv = (fisher.pvalue(a, b, c, d).right_tail)+1e-300;  # change
            #pv = fisher_exact([[a,b],[c,d]])[1];
            score = ((((-np.log10(pv))*ra)*gene_mean)*100)/num_cell_set;
            all_pv.append(pv);
            all_ratio.append(ra);
            all_score.append(score);
            all_means.append(gene_mean);
    
        all_qv = bh_qvalues(all_pv);
        all_qv = [i+1e-300 for i in all_qv];
        
        gene_pv[gene] = all_pv;
        gene_qv[gene] = all_qv;
        gene_ratio[gene] = all_ratio;
        gene_score[gene] = all_score;
        gene_means[gene] = all_means;
        
        return (gene_ratio, gene_pv, gene_qv, gene_score,gene_means);

def bh_qvalues(pv):
    if pv == []:
        return [];
    m = len(pv);
    args,pv = zip(*sorted(enumerate(pv),key=operator.itemgetter(1)));
    if pv[0] < 0 or pv[-1] >1:
        raise ValueError("p-values must be between 0 and 1");
    qvalues = m*[0];
    mincoeff = pv[-1];
    qvalues[args[-1]] = mincoeff;
    for j in range(m-2,-1,-1):
        coeff = m*pv[j]/float(j+1);
        if coeff < mincoeff:
            mincoeff = coeff;
        qvalues[args[j]] = mincoeff;
    
    return qvalues
     
# get highly expressed cells for given gene
def get_highly_cells_for_each_gene(rdata_df: 'dataFrame of raw data', gene: 'gene name',power=11):
    """
    get highly expressed cells for given gene.
    
    Parameters
    ---------
    rdata_df:data frame format of raw annoData;
    gene : a given gene name;
    """
    power=power;
    gene_list = rdata_df[gene];
    gene_filter = gene_list.loc[gene_list>0,];
    gene_sort_df_filter = gene_filter.sort_values(ascending=False);
    
    gene_highly_cells = dict();
    gene_mean_exvalue = dict();
    #gene_sort_df = rdata_df.sort_values(by=gene,ascending=False);
    #gene_sort_df_filter = gene_sort_df.loc[gene_sort_df[gene]>0,];
    # get the cells which highly express given gene
    if len(gene_sort_df_filter)<=5:
        #print('the '+gene+' has less than 10 highly expressed cells ');
        gene_highly_cells[gene] = [];
        gene_mean_exvalue[gene] = [];
        return (gene_highly_cells,gene_mean_exvalue);
    #Choose the minimum real root and use the "try...except..." to jump the "ValueError".Corrected code as follow.
    else:
        matAA = curve_fitting_2(gene_sort_df_filter,power=power);
        df_3_solution = calculate_derivation(matAA);
        df_3_solution = [complex(x) for x in df_3_solution];
        real_roots = [int(x.real) for x in df_3_solution if x.imag==0];
        real_roots.sort();
        real_roots_p = [x for x in real_roots if x>0];
        if len(real_roots_p)>0:
            highly_cells = list(gene_sort_df_filter.index[0:real_roots_p[0]]);
            gene_mean = gene_sort_df_filter.iloc[0:real_roots_p[0]].mean();
            gene_highly_cells[gene] = highly_cells;
            gene_mean_exvalue[gene] = gene_mean;
        else:
            gene_highly_cells[gene] = [];
            gene_mean_exvalue[gene] = [];
        
        return (gene_highly_cells,gene_mean_exvalue);


def calculate_genes_pseudotime_location(rdata_df: 'dataFrame of annoData',
                                        cell_sets,
                                        markers_s_LGPS_rdata,
                                        adata_obsm_df,
                                        gene: 'gene name',
                                       power=11):
    """
    Calculate the pseudotemporal location and ordering of a gaven gene.
    
    Parameters
    ---------
    rdata_df: data frame format of raw annoData;
    cell_sets: a dict. key is the cell tyoe; value is a list of cell name;
    markers_s_LGPS_rdata: markers matrix from results by get_DEG_single or get_DEG_multiple';
    adata_obsm_df: a data frame which contains the one of 'obsm' attribution of adata;
    gene : a given gene name;
    """
    gene = gene;
    gene_pseudotime_locate = list();
    gene_highly_cells = get_highly_cells_for_each_gene(rdata_df,gene,power)[0][gene];
    
    if gene_highly_cells ==[]:
        pass;
    else:
        if gene in list(markers_s_LGPS_rdata['gene_name']):
            gene_cell_types = list(markers_s_LGPS_rdata.loc[markers_s_LGPS_rdata['gene_name']==gene,'cell_type']);
            gene_ratio = list(markers_s_LGPS_rdata.loc[markers_s_LGPS_rdata['gene_name']==gene,'ratio']);
            gene_Qvalue = list(markers_s_LGPS_rdata.loc[markers_s_LGPS_rdata['gene_name']==gene,'q_value']);
            gene_mean_exValue = list(markers_s_LGPS_rdata.loc[markers_s_LGPS_rdata['gene_name']==gene,'mean_exValue']);
            for i in range(len(gene_cell_types)):
                CT = gene_cell_types[i];
                RT = gene_ratio[i];
                QV = gene_Qvalue[i];
                MEV = gene_mean_exValue[i];
                gene_pseudotime_locate_dict = dict();
                gene_cell_type_celllist = cell_sets[CT];
    
                gene_cell_type_highly_cell = list(set(gene_highly_cells) & set(gene_cell_type_celllist));
    
                gene_pseudotime_cells_df = adata_obsm_df.loc[adata_obsm_df.index.isin(gene_cell_type_highly_cell),:];
                gene_pseudotime = np.median(gene_pseudotime_cells_df['dpt_pseudotime']); 
                # calculating the x and y location for each gene
                gene_pseu_dist = abs(gene_pseudotime_cells_df['dpt_pseudotime']-gene_pseudotime);
                gene_min_pseu_dist_cell = gene_pseu_dist.idxmin();
        
                if len(gene_pseudotime_cells_df.columns) >= 4:
                    gene_x_locate = gene_pseudotime_cells_df.loc[gene_min_pseu_dist_cell,0];
                    gene_y_locate = gene_pseudotime_cells_df.loc[gene_min_pseu_dist_cell,1];
                    gene_z_locate = gene_pseudotime_cells_df.loc[gene_min_pseu_dist_cell,2];
            
                    gene_pseudotime_locate_dict['gene_name'] = gene;
                    gene_pseudotime_locate_dict['x_location'] = gene_x_locate;
                    gene_pseudotime_locate_dict['y_location'] = gene_y_locate;
                    gene_pseudotime_locate_dict['z_location'] = gene_z_locate;
                    gene_pseudotime_locate_dict['dpt_pseudotime'] = gene_pseudotime;
                    gene_pseudotime_locate_dict['cell_type'] = CT;
                    gene_pseudotime_locate_dict['ratio'] = RT;
                    gene_pseudotime_locate_dict['q_value'] = QV;
                    gene_pseudotime_locate_dict['mean_exValue'] = MEV;
        
                if len(gene_pseudotime_cells_df.columns) == 3:
                    gene_x_locate = gene_pseudotime_cells_df.loc[gene_min_pseu_dist_cell,0];
                    gene_y_locate = gene_pseudotime_cells_df.loc[gene_min_pseu_dist_cell,1];
            
                    gene_pseudotime_locate_dict['gene_name'] = gene;
                    gene_pseudotime_locate_dict['x_location'] = gene_x_locate;
                    gene_pseudotime_locate_dict['y_location'] = gene_y_locate;
                    gene_pseudotime_locate_dict['dpt_pseudotime'] = gene_pseudotime;
                    gene_pseudotime_locate_dict['cell_type'] = CT;
                    gene_pseudotime_locate_dict['ratio'] = RT;
                    gene_pseudotime_locate_dict['q_value'] = QV;
                    gene_pseudotime_locate_dict['mean_exValue'] = MEV;
                
                gene_pseudotime_locate.append(gene_pseudotime_locate_dict);
        else:
            gene_pseudotime_locate_dict = dict();
            gene_cell_type_highly_cell = list(set(gene_highly_cells));
    
            gene_pseudotime_cells_df = adata_obsm_df.loc[adata_obsm_df.index.isin(gene_cell_type_highly_cell),:];
            gene_pseudotime = np.median(gene_pseudotime_cells_df['dpt_pseudotime']); 
            # calculating the x and y location for each gene
            gene_pseu_dist = abs(gene_pseudotime_cells_df['dpt_pseudotime']-gene_pseudotime);
            gene_min_pseu_dist_cell = gene_pseu_dist.idxmin();
        
            if len(gene_pseudotime_cells_df.columns) >= 4:
                gene_x_locate = gene_pseudotime_cells_df.loc[gene_min_pseu_dist_cell,0];
                gene_y_locate = gene_pseudotime_cells_df.loc[gene_min_pseu_dist_cell,1];
                gene_z_locate = gene_pseudotime_cells_df.loc[gene_min_pseu_dist_cell,2];
            
                gene_pseudotime_locate_dict['gene_name'] = gene;
                gene_pseudotime_locate_dict['x_location'] = gene_x_locate;
                gene_pseudotime_locate_dict['y_location'] = gene_y_locate;
                gene_pseudotime_locate_dict['z_location'] = gene_z_locate;
                gene_pseudotime_locate_dict['dpt_pseudotime'] = gene_pseudotime;
                gene_pseudotime_locate_dict['cell_type'] = 'None';
                gene_pseudotime_locate_dict['ratio'] = 0;
                gene_pseudotime_locate_dict['q_value'] = 1;
                gene_pseudotime_locate_dict['mean_exValue'] = 0;
        
            if len(gene_pseudotime_cells_df.columns) == 3:
                gene_x_locate = gene_pseudotime_cells_df.loc[gene_min_pseu_dist_cell,0];
                gene_y_locate = gene_pseudotime_cells_df.loc[gene_min_pseu_dist_cell,1];
            
                gene_pseudotime_locate_dict['gene_name'] = gene;
                gene_pseudotime_locate_dict['x_location'] = gene_x_locate;
                gene_pseudotime_locate_dict['y_location'] = gene_y_locate;
                gene_pseudotime_locate_dict['dpt_pseudotime'] = gene_pseudotime;
                gene_pseudotime_locate_dict['cell_type'] = 'NS';
                gene_pseudotime_locate_dict['ratio'] = 0;
                gene_pseudotime_locate_dict['q_value'] = 1;
                gene_pseudotime_locate_dict['mean_exValue'] = 0;
                
            gene_pseudotime_locate.append(gene_pseudotime_locate_dict);
    return gene_pseudotime_locate;

def projection(A,b):
    AA = A.T.dot(A)
    w = np.linalg.inv(AA).dot(A.T).dot(b)
    #print(w)
    return w
def curve_fitting_2(sort_data: 'sorted dataFrame of annoData',power=11):
    if power:
        matAA = curve_fitting(sort_data,power=power);
        return matAA;
    else:
        m = len(sort_data);
        frrs = [];
        xa = [x for x in range(1,m+1)];
        ya = list(sort_data);
        min_frr,min_deg = 1e10, 0;
        if m<=50:
            degrees = np.arange(1,int(m/2));
            for deg in degrees:
                matAA = curve_fitting(sort_data,power=deg);
                xxa=[x for x in range(1,m+1)];
                yya=[];
                for i in range(0,len(xxa)):
                    yy=0.0;
                    for j in range(0,deg):
                        dy=1.0;
                        for k in range(0,j):
                            dy*=xxa[i];
                        dy*=matAA[j];
                        yy+=dy;
                    yya.append(yy);
                poly_rmse = np.sqrt(mean_squared_error(ya, yya));
                poly_frr = poly_rmse;
                frrs.append(poly_frr);
                if min_frr > poly_frr:
                    min_frr = poly_frr;
                    min_deg = deg;
            power = min_deg;
            matAA = curve_fitting(sort_data,power=power);
            return matAA;
        else:
            degrees = np.arange(1, 26);
            for deg in degrees:
                matAA = curve_fitting(sort_data,power=deg);
                xxa=[x for x in range(1,m+1)];
                yya=[];
                for i in range(0,len(xxa)):
                    yy=0.0;
                    for j in range(0,deg):
                        dy=1.0;
                        for k in range(0,j):
                            dy*=xxa[i];
                        dy*=matAA[j];
                        yy+=dy;
                    yya.append(yy);
                poly_rmse = np.sqrt(mean_squared_error(ya, yya));
                poly_frr = np.sqrt(poly_rmse*poly_rmse*m/(m-deg-1));
                frrs.append(poly_frr);
                if min_frr > poly_frr:
                    min_frr = poly_frr;
                    min_deg = deg;
            power = min_deg;
            matAA = curve_fitting(sort_data,power=power);
            return matAA;

# curve fitting using ordinary least squares techniques
def curve_fitting(sort_data: 'sorted dataFrame of annoData',power=11):
    xa = np.arange(1,len(sort_data)+1,dtype=float)
    #xa=[x for x in range(1,len(sort_data[gene])+1)]
    ya=list(sort_data)
    #xa = np.array(xa)
    #xa.dtype="float64"
    y1 = np.array(ya)
    # power=9; user can set different power
    power=power;#recommend power=8;

    m = []

    for i in range(power):
        a = xa**(i)
        m.append(a)
    A = np.array(m).T
    b = y1.reshape(y1.shape[0],1)

    matAA = projection(A,b)
    matAA.shape = power,

    return matAA

def calculate_derivation(matAA):
    f = np.poly1d(list(reversed(matAA)));
    df_1 = f.deriv(1);
    df_2 = f.deriv(2);
    df_3 = f.deriv(3);
    curvature_1d = (1+df_1**2)*df_3-3*df_1*df_2*df_2;
    #curvature_solution = np.roots(df_3);
    curvature_solution = np.roots(curvature_1d);
    return curvature_solution;

def draw_fitted_curve(gene_sort_df_filter,power):
    #draw the curve after fitting
    power=power;
    xa=[x for x in range(1,len(gene_sort_df_filter)+1)]
    ya=list(gene_sort_df_filter)
    matAA = curve_fitting(gene_sort_df_filter,power);
    xxa=[x for x in range(1,len(gene_sort_df_filter)+1)]
    yya=[]
    for i in range(0,len(xxa)):
        yy=0.0
        for j in range(0,power):
            dy=1.0
            for k in range(0,j):
                dy*=xxa[i]
            dy*=matAA[j]
            yy+=dy
        yya.append(yy)
    fig=plt.figure()
    ax=fig.add_subplot(111)
    ax.plot(xa,ya,color='green',linestyle='',marker='o',markersize=4)
    ax.plot(xxa,yya,color='red',linestyle='-',marker='',linewidth=2)
    return fig

def bfs_edges_modified(tree, source):
    visited, queue = [], [source]
    bfs_tree = nx.bfs_tree(tree,source=source)
    predecessors = dict(nx.bfs_predecessors(bfs_tree,source))
    edges = []
    while queue:
        vertex = queue.pop()
        if vertex not in visited:
            visited.append(vertex)
            if(vertex in predecessors.keys()):
                edges.append((predecessors[vertex],vertex))
            unvisited = set(tree[vertex]) - set(visited)
        queue.extend(unvisited)
    return edges