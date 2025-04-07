import scanpy as sc
import anndata
import os
from scipy import sparse
import numpy as np
import pandas as pd
import json
import sys
from typing import Dict, Any
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
from scipy.stats import fisher_exact, hypergeom, binom
from random import shuffle
import operator
import random
import math
from sklearn.metrics import mean_squared_error, r2_score
from sympy import *
from multiprocessing import Pool
from functools import partial
import datetime
# import fisher
import networkx as nx
import itertools

def projection(A, b):
    AA = A.T.dot(A)
    w = np.linalg.inv(AA).dot(A.T).dot(b)
    return w

def curve_fitting(sort_data, power=11):
    xa = np.arange(1, len(sort_data)+1, dtype=float)
    ya = list(sort_data)
    y1 = np.array(ya)
    power = power

    m = []
    for i in range(power):
        a = xa**(i)
        m.append(a)
    A = np.array(m).T
    b = y1.reshape(y1.shape[0], 1)

    matAA = projection(A, b)
    matAA.shape = power,
    return matAA

def calculate_derivation(matAA):
    f = np.poly1d(list(reversed(matAA)))
    df_1 = f.deriv(1)
    df_2 = f.deriv(2)
    df_3 = f.deriv(3)
    curvature_1d = (1+df_1**2)*df_3-3*df_1*df_2*df_2
    curvature_solution = np.roots(curvature_1d)
    return curvature_solution

def curve_fitting_2(sort_data, power=11):
    if power:
        matAA = curve_fitting(sort_data, power=power)
        return matAA
    else:
        m = len(sort_data)
        frrs = []
        xa = [x for x in range(1, m+1)]
        ya = list(sort_data)
        min_frr, min_deg = 1e10, 0
        if m <= 50:
            degrees = np.arange(1, int(m/2))
            for deg in degrees:
                matAA = curve_fitting(sort_data, power=deg)
                xxa = [x for x in range(1, m+1)]
                yya = []
                for i in range(0, len(xxa)):
                    yy = 0.0
                    for j in range(0, deg):
                        dy = 1.0
                        for k in range(0, j):
                            dy *= xxa[i]
                        dy *= matAA[j]
                        yy += dy
                    yya.append(yy)
                poly_rmse = np.sqrt(mean_squared_error(ya, yya))
                poly_frr = poly_rmse
                frrs.append(poly_frr)
                if min_frr > poly_frr:
                    min_frr = poly_frr
                    min_deg = deg
            power = min_deg
            matAA = curve_fitting(sort_data, power=power)
            return matAA
        else:
            degrees = np.arange(1, 26)
            for deg in degrees:
                matAA = curve_fitting(sort_data, power=deg)
                xxa = [x for x in range(1, m+1)]
                yya = []
                for i in range(0, len(xxa)):
                    yy = 0.0
                    for j in range(0, deg):
                        dy = 1.0
                        for k in range(0, j):
                            dy *= xxa[i]
                        dy *= matAA[j]
                        yy += dy
                    yya.append(yy)
                poly_rmse = np.sqrt(mean_squared_error(ya, yya))
                poly_frr = np.sqrt(poly_rmse*poly_rmse*m/(m-deg-1))
                frrs.append(poly_frr)
                if min_frr > poly_frr:
                    min_frr = poly_frr
                    min_deg = deg
            power = min_deg
            matAA = curve_fitting(sort_data, power=power)
            return matAA

def get_highly_cells_for_each_gene(rdata_df, gene, power=11):
    power = power
    gene_list = rdata_df[gene]
    gene_filter = gene_list.loc[gene_list > 0]
    gene_sort_df_filter = gene_filter.sort_values(ascending=False)
    
    gene_highly_cells = dict()
    gene_mean_exvalue = dict()
    
    if len(gene_sort_df_filter) <= 5:
        gene_highly_cells[gene] = []
        gene_mean_exvalue[gene] = []
        return (gene_highly_cells, gene_mean_exvalue)
    else:
        matAA = curve_fitting_2(gene_sort_df_filter, power=power)
        df_3_solution = calculate_derivation(matAA)
        df_3_solution = [complex(x) for x in df_3_solution]
        real_roots = [int(x.real) for x in df_3_solution if x.imag == 0]
        real_roots.sort()
        real_roots_p = [x for x in real_roots if x > 0]
        if len(real_roots_p) > 0:
            highly_cells = list(gene_sort_df_filter.index[0:real_roots_p[0]])
            gene_mean = gene_sort_df_filter.iloc[0:real_roots_p[0]].mean()
            gene_highly_cells[gene] = highly_cells
            gene_mean_exvalue[gene] = gene_mean
        else:
            gene_highly_cells[gene] = []
            gene_mean_exvalue[gene] = []
        
        return (gene_highly_cells, gene_mean_exvalue)

def bh_qvalues(pv):
    if pv == []:
        return []
    m = len(pv)
    args, pv = zip(*sorted(enumerate(pv), key=operator.itemgetter(1)))
    if pv[0] < 0 or pv[-1] > 1:
        raise ValueError("p-values must be between 0 and 1")
    qvalues = m*[0]
    mincoeff = pv[-1]
    qvalues[args[-1]] = mincoeff
    for j in range(m-2, -1, -1):
        coeff = m*pv[j]/float(j+1)
        if coeff < mincoeff:
            mincoeff = coeff
        qvalues[args[j]] = mincoeff
    
    return qvalues

def Fisher_test_for_each_gene(rdata_df, cell_sets, num_allCells, gene, power=11):
    power = power
    gene = gene
    gene_highly_cells, gene_mean_exvalue = get_highly_cells_for_each_gene(rdata_df, gene, power)
    test_cells = gene_highly_cells[gene]
    gene_mean = gene_mean_exvalue[gene]
    num_test_cells = len(test_cells)
    
    gene_pv = dict()
    gene_qv = dict()
    gene_ratio = dict()
    gene_score = dict()
    gene_means = dict()
    
    all_ratio = []
    all_pv = []
    all_qv = []
    all_score = []
    all_means = []
    
    if num_test_cells == 0:
        gene_ratio[gene] = []
        gene_pv[gene] = []
        gene_qv[gene] = []
        gene_score[gene] = []
        gene_means[gene] = []
        
        return (gene_ratio, gene_pv, gene_qv, gene_score, gene_means)
    else:
        for cs in cell_sets.keys():
            cell_set_cs = cell_sets[cs]
            num_cell_set = len(cell_set_cs)
            a = len(set(test_cells) & set(cell_set_cs))
            b = num_test_cells-a
            c = num_cell_set-a
            d = num_allCells-num_cell_set
            ra = a/float(num_test_cells)
            # pv = (fisher.pvalue(a, b, c, d).right_tail)+1e-300
            _, pv = fisher_exact([[a, b], [c, d]], alternative='greater')
            pv = pv + 1e-300
            score = ((((-np.log10(pv))*ra)*gene_mean)*100)/num_cell_set
            all_pv.append(pv)
            all_ratio.append(ra)
            all_score.append(score)
            all_means.append(gene_mean)
    
        all_qv = bh_qvalues(all_pv)
        all_qv = [i+1e-300 for i in all_qv]
        
        gene_pv[gene] = all_pv
        gene_qv[gene] = all_qv
        gene_ratio[gene] = all_ratio
        gene_score[gene] = all_score
        gene_means[gene] = all_means
        
        return (gene_ratio, gene_pv, gene_qv, gene_score, gene_means)

def get_DEG_uniq(rdata, adata, group_key, power=11, ratio=0.2, p_threshold=0.01, q_threshold=0.05):
    power = power
    print('power: ', power)
    print('get the raw data frame...')
    rdata_df = rdata.to_df()
    num_allCells = len(rdata_df.index)
    
    print('struct the cell type sets for enrichment analysis...')
    adata_cell_type_df = pd.DataFrame(adata.obs[group_key])
    cell_type_index = pd.Categorical(adata.obs[group_key]).categories
    
    cell_sets = dict()
    
    for ct in cell_type_index:
        cell_sets[ct] = list(adata_cell_type_df.loc[adata_cell_type_df[group_key] == ct, :].index)
    
    iter_genes = rdata_df.columns
    
    print('Fisher_test_for_each_gene...')
      
    genes_ratio = dict()
    genes_pv = dict()
    genes_qv = dict()
    genes_score = dict()
    genes_means = dict()
    
    num = 0
    for gene in iter_genes:
        num += 1
        gene_ratio, gene_pv, gene_qv, gene_score, gene_means = Fisher_test_for_each_gene(rdata_df, cell_sets, num_allCells, gene, power)
        if gene_ratio[gene] == []:
            continue
        else:
            genes_ratio = dict(genes_ratio, **gene_ratio)
            genes_pv = dict(genes_pv, **gene_pv)
            genes_qv = dict(genes_qv, **gene_qv)
            genes_score = dict(genes_score, **gene_score)
            genes_means = dict(genes_means, **gene_means)
        if num % 1000 == 0:
            print('whole ', num, ' genes have been done.')

    print('merge differentially expressed genes...')
    genes_ratio_df = pd.DataFrame(genes_ratio)
    genes_ratio_df.index = cell_type_index
    
    genes_pv_df = pd.DataFrame(genes_pv)
    genes_pv_df.index = cell_type_index
    
    genes_qv_df = pd.DataFrame(genes_qv)
    genes_qv_df.index = cell_type_index
    
    genes_score_df = pd.DataFrame(genes_score)
    genes_score_df.index = cell_type_index
    
    genes_means_df = pd.DataFrame(genes_means)
    genes_means_df.index = cell_type_index
    
    ra = ratio
    pt = p_threshold
    qt = q_threshold

    ct_gene_r_p_q_s_list = []
    
    for gene in genes_ratio_df.columns:
        gene = gene
        gene_on_ct = genes_score_df.loc[:, gene].idxmax(axis=0)
        
        if genes_ratio_df.loc[gene_on_ct, gene] >= ra:
            if genes_pv_df.loc[gene_on_ct, gene] <= pt and genes_qv_df.loc[gene_on_ct, gene] <= qt:
                ct = gene_on_ct
                g = gene
                r = genes_ratio_df.loc[gene_on_ct, gene]
                p = genes_pv_df.loc[gene_on_ct, gene]
                q = genes_qv_df.loc[gene_on_ct, gene]
                s = genes_score_df.loc[gene_on_ct, gene]
                m = genes_means_df.loc[gene_on_ct, gene]
                
                ct_gene_r_p_q_s_list.append([ct, g, r, p, q, s, m])
            else:
                continue
        else:
            continue
    markers_s = pd.DataFrame(ct_gene_r_p_q_s_list, columns=['cell_type', 'gene_name', 'ratio', 'p_value', 'q_value', 'score', 'mean_exValue'])
    print('Done!')
    return markers_s

def get_DEG_multi(rdata, adata, group_key, power=11, ratio=0.2, p_threshold=0.01, q_threshold=0.05):
    power = power
    print('get the raw data frame...')
    rdata_df = rdata.to_df()
    num_allCells = len(rdata.obs_names)
    
    print('struct the cell type sets for enrichment analysis...')
    adata_cell_type_df = pd.DataFrame(adata.obs[group_key])
    cell_type_index = pd.Categorical(adata.obs[group_key]).categories
    
    cell_sets = dict()
    
    for ct in cell_type_index:
        cell_sets[ct] = list(adata_cell_type_df.loc[adata_cell_type_df[group_key] == ct, :].index)
    
    iter_genes = rdata_df.columns
    
    print('Fisher_test_for_each_gene...')
      
    genes_ratio = dict()
    genes_pv = dict()
    genes_qv = dict()
    genes_score = dict()
    genes_means = dict()
    
    num = 0
    for gene in iter_genes:
        num += 1
        gene_ratio, gene_pv, gene_qv, gene_score, gene_means = Fisher_test_for_each_gene(rdata_df, cell_sets, num_allCells, gene, power)
        if gene_ratio[gene] == []:
            continue
        else:
            genes_ratio = dict(genes_ratio, **gene_ratio)
            genes_pv = dict(genes_pv, **gene_pv)
            genes_qv = dict(genes_qv, **gene_qv)
            genes_score = dict(genes_score, **gene_score)
            genes_means = dict(genes_means, **gene_means)
            
        if num % 1000 == 0:
            print('whole ', num, ' genes have been done.')

    print('merge differentially expressed genes...')
    genes_ratio_df = pd.DataFrame(genes_ratio)
    genes_ratio_df.index = cell_type_index
    genes_ratio_df = genes_ratio_df.T
    
    genes_pv_df = pd.DataFrame(genes_pv)
    genes_pv_df.index = cell_type_index
    genes_pv_df = genes_pv_df.T
    
    genes_qv_df = pd.DataFrame(genes_qv)
    genes_qv_df.index = cell_type_index
    genes_qv_df = genes_qv_df.T
    
    genes_score_df = pd.DataFrame(genes_score)
    genes_score_df.index = cell_type_index
    genes_score_df = genes_score_df.T
    
    genes_means_df = pd.DataFrame(genes_means)
    genes_means_df.index = cell_type_index
    genes_means_df = genes_means_df.T
    
    ra = ratio
    pt = p_threshold
    qt = q_threshold
    deg_genes_ratio_p_q_s_list = []

    for ct in cell_type_index:
        up_ratio_ct_genes_list = list(genes_ratio_df.loc[genes_ratio_df[ct] >= ra, :].index)
        down_pv_ct_genes_list = list(genes_pv_df.loc[genes_pv_df[ct] <= pt, :].index)
        down_qv_ct_genes_list = list(genes_qv_df.loc[genes_qv_df[ct] <= qt, :].index)
        
        ct_deg_genes_list = list(set(up_ratio_ct_genes_list) & set(down_pv_ct_genes_list) & set(down_qv_ct_genes_list))
        
        if ct_deg_genes_list == []:
            continue
            
        else:
            for gene in ct_deg_genes_list:
                ct = ct
                g = gene
                r = genes_ratio_df.loc[gene, ct]
                p = genes_pv_df.loc[gene, ct]
                q = genes_qv_df.loc[gene, ct]
                s = genes_score_df.loc[gene, ct]
                m = genes_means_df.loc[gene, ct]
                
                deg_genes_ratio_p_q_s_list.append([ct, g, r, p, q, s, m])
    
    markers_m = pd.DataFrame(deg_genes_ratio_p_q_s_list, columns=['cell_type', 'gene_name', 'ratio', 'p_value', 'q_value', 'score', 'mean_exValue'])
    print('Done!')
    return markers_m

def calculate_genes_pseudotime_location(rdata_df, cell_sets, markers_s_LGPS_rdata, adata_obsm_df, gene, power=11):
    gene = gene
    gene_pseudotime_locate = list()
    gene_highly_cells = get_highly_cells_for_each_gene(rdata_df, gene, power)[0][gene]
    
    if gene_highly_cells == []:
        pass
    else:
        if gene in list(markers_s_LGPS_rdata['gene_name']):
            gene_cell_types = list(markers_s_LGPS_rdata.loc[markers_s_LGPS_rdata['gene_name'] == gene, 'cell_type'])
            gene_ratio = list(markers_s_LGPS_rdata.loc[markers_s_LGPS_rdata['gene_name'] == gene, 'ratio'])
            gene_Qvalue = list(markers_s_LGPS_rdata.loc[markers_s_LGPS_rdata['gene_name'] == gene, 'q_value'])
            gene_mean_exValue = list(markers_s_LGPS_rdata.loc[markers_s_LGPS_rdata['gene_name'] == gene, 'mean_exValue'])
            for i in range(len(gene_cell_types)):
                CT = gene_cell_types[i]
                RT = gene_ratio[i]
                QV = gene_Qvalue[i]
                MEV = gene_mean_exValue[i]
                gene_pseudotime_locate_dict = dict()
                gene_cell_type_celllist = cell_sets[CT]
    
                gene_cell_type_highly_cell = list(set(gene_highly_cells) & set(gene_cell_type_celllist))
                
                if len(gene_cell_type_highly_cell) == 0:
                    continue
    
                gene_pseudotime_cells_df = adata_obsm_df.loc[adata_obsm_df.index.isin(gene_cell_type_highly_cell), :]
                
                if len(gene_pseudotime_cells_df) == 0:
                    continue
                    
                gene_pseudotime = np.nanmedian(gene_pseudotime_cells_df['dpt_pseudotime'])
                if np.isnan(gene_pseudotime):
                    continue
                    
                gene_pseu_dist = abs(gene_pseudotime_cells_df['dpt_pseudotime']-gene_pseudotime)
                if gene_pseu_dist.isna().all():
                    continue
                    
                gene_min_pseu_dist_cell = gene_pseu_dist.idxmin()
                if pd.isna(gene_min_pseu_dist_cell):
                    continue
        
                if len(gene_pseudotime_cells_df.columns) >= 4:
                    try:
                        gene_x_locate = gene_pseudotime_cells_df.loc[gene_min_pseu_dist_cell, 0]
                        gene_y_locate = gene_pseudotime_cells_df.loc[gene_min_pseu_dist_cell, 1]
                        gene_z_locate = gene_pseudotime_cells_df.loc[gene_min_pseu_dist_cell, 2]
                        
                        if np.isnan(gene_x_locate) or np.isnan(gene_y_locate) or np.isnan(gene_z_locate):
                            continue
                            
                        gene_pseudotime_locate_dict['gene_name'] = gene
                        gene_pseudotime_locate_dict['x_location'] = gene_x_locate
                        gene_pseudotime_locate_dict['y_location'] = gene_y_locate
                        gene_pseudotime_locate_dict['z_location'] = gene_z_locate
                        gene_pseudotime_locate_dict['dpt_pseudotime'] = gene_pseudotime
                        gene_pseudotime_locate_dict['cell_type'] = CT
                        gene_pseudotime_locate_dict['ratio'] = RT
                        gene_pseudotime_locate_dict['q_value'] = QV
                        gene_pseudotime_locate_dict['mean_exValue'] = MEV
                    except:
                        continue
        
                if len(gene_pseudotime_cells_df.columns) == 3:
                    try:
                        gene_x_locate = gene_pseudotime_cells_df.loc[gene_min_pseu_dist_cell, 0]
                        gene_y_locate = gene_pseudotime_cells_df.loc[gene_min_pseu_dist_cell, 1]
                        
                        if np.isnan(gene_x_locate) or np.isnan(gene_y_locate):
                            continue
                            
                        gene_pseudotime_locate_dict['gene_name'] = gene
                        gene_pseudotime_locate_dict['x_location'] = gene_x_locate
                        gene_pseudotime_locate_dict['y_location'] = gene_y_locate
                        gene_pseudotime_locate_dict['dpt_pseudotime'] = gene_pseudotime
                        gene_pseudotime_locate_dict['cell_type'] = CT
                        gene_pseudotime_locate_dict['ratio'] = RT
                        gene_pseudotime_locate_dict['q_value'] = QV
                        gene_pseudotime_locate_dict['mean_exValue'] = MEV
                    except:
                        continue
                
                gene_pseudotime_locate.append(gene_pseudotime_locate_dict)
        else:
            gene_pseudotime_locate_dict = dict()
            gene_cell_type_highly_cell = list(set(gene_highly_cells))
            
            if len(gene_cell_type_highly_cell) == 0:
                return gene_pseudotime_locate
    
            gene_pseudotime_cells_df = adata_obsm_df.loc[adata_obsm_df.index.isin(gene_cell_type_highly_cell), :]
            
            if len(gene_pseudotime_cells_df) == 0:
                return gene_pseudotime_locate
                
            gene_pseudotime = np.nanmedian(gene_pseudotime_cells_df['dpt_pseudotime'])
            if np.isnan(gene_pseudotime):
                return gene_pseudotime_locate
                
            gene_pseu_dist = abs(gene_pseudotime_cells_df['dpt_pseudotime']-gene_pseudotime)
            if gene_pseu_dist.isna().all():
                return gene_pseudotime_locate
                
            gene_min_pseu_dist_cell = gene_pseu_dist.idxmin()
            if pd.isna(gene_min_pseu_dist_cell):
                return gene_pseudotime_locate
        
            if len(gene_pseudotime_cells_df.columns) >= 4:
                try:
                    gene_x_locate = gene_pseudotime_cells_df.loc[gene_min_pseu_dist_cell, 0]
                    gene_y_locate = gene_pseudotime_cells_df.loc[gene_min_pseu_dist_cell, 1]
                    gene_z_locate = gene_pseudotime_cells_df.loc[gene_min_pseu_dist_cell, 2]
                    
                    if np.isnan(gene_x_locate) or np.isnan(gene_y_locate) or np.isnan(gene_z_locate):
                        return gene_pseudotime_locate
                        
                    gene_pseudotime_locate_dict['gene_name'] = gene
                    gene_pseudotime_locate_dict['x_location'] = gene_x_locate
                    gene_pseudotime_locate_dict['y_location'] = gene_y_locate
                    gene_pseudotime_locate_dict['z_location'] = gene_z_locate
                    gene_pseudotime_locate_dict['dpt_pseudotime'] = gene_pseudotime
                    gene_pseudotime_locate_dict['cell_type'] = 'None'
                    gene_pseudotime_locate_dict['ratio'] = 0
                    gene_pseudotime_locate_dict['q_value'] = 1
                    gene_pseudotime_locate_dict['mean_exValue'] = 0
                except:
                    return gene_pseudotime_locate
        
            if len(gene_pseudotime_cells_df.columns) == 3:
                try:
                    gene_x_locate = gene_pseudotime_cells_df.loc[gene_min_pseu_dist_cell, 0]
                    gene_y_locate = gene_pseudotime_cells_df.loc[gene_min_pseu_dist_cell, 1]
                    
                    if np.isnan(gene_x_locate) or np.isnan(gene_y_locate):
                        return gene_pseudotime_locate
                        
                    gene_pseudotime_locate_dict['gene_name'] = gene
                    gene_pseudotime_locate_dict['x_location'] = gene_x_locate
                    gene_pseudotime_locate_dict['y_location'] = gene_y_locate
                    gene_pseudotime_locate_dict['dpt_pseudotime'] = gene_pseudotime
                    gene_pseudotime_locate_dict['cell_type'] = 'NS'
                    gene_pseudotime_locate_dict['ratio'] = 0
                    gene_pseudotime_locate_dict['q_value'] = 1
                    gene_pseudotime_locate_dict['mean_exValue'] = 0
                except:
                    return gene_pseudotime_locate
                
            gene_pseudotime_locate.append(gene_pseudotime_locate_dict)
    return gene_pseudotime_locate

def get_genes_location_pseudotime(rdata, adata, group_key, gene_matrix, obsm, power=11):
    start = datetime.datetime.now()
    
    rdata_df = rdata.to_df()
    
    adata_obsm_df = pd.DataFrame(adata.obsm[obsm])
    adata_obsm_df.index = adata.obs_names
    adata_obsm_df['dpt_pseudotime'] = adata.obs.dpt_pseudotime
    
    markers_s_LGPS = adata.uns[gene_matrix]
    markers_s_LGPS_rdata = markers_s_LGPS.loc[markers_s_LGPS['gene_name'].isin(list(rdata.var_names)), ]
    
    adata_cell_type_df = pd.DataFrame(adata.obs[group_key])
    cell_type_index = pd.Categorical(adata.obs[group_key]).categories
    
    cell_sets = dict()
    
    for ct in cell_type_index:
        cell_sets[ct] = list(adata_cell_type_df.loc[adata_cell_type_df[group_key] == ct, :].index)

    iter_genes = rdata_df.columns
    gene_pseudotime_locates = list()
    num = 0
    for g1 in iter_genes:
        num += 1
        g1_pseudotime_locates = calculate_genes_pseudotime_location(rdata_df,
                                                                    cell_sets,
                                                                    markers_s_LGPS_rdata,
                                                                    adata_obsm_df,
                                                                    g1,
                                                                   power)
        if g1_pseudotime_locates == []:
            continue
        else:
            for i in range(0, len(g1_pseudotime_locates)):
                gene_pseudotime_locates.append(g1_pseudotime_locates[i])
        
        if num % 1000 == 0:
            print('whole ', num, ' genes have been done.')
    
    gene_pseudotime_locates_df = pd.DataFrame(gene_pseudotime_locates)
    
    end = datetime.datetime.now()
    print('Running time : %s Seconds' % (end-start))
    print('Done!')
    return gene_pseudotime_locates_df

def run_deaplog_analysis(rdata, adata, sample_percent=None):
    """
    Run DEAPLOG analysis on the data and return results in format suitable for frontend
    """
    try:
        print(f"Initial data shape: {adata.shape}")
        print(f"Sample percent received: {sample_percent}")
        
        # Sample cells if specified
        if sample_percent is not None:
            sample_percent = float(sample_percent)
            if sample_percent < 1.0:
                n_cells = int(adata.n_obs * sample_percent)
                print(f"Sampling {n_cells} cells from {adata.n_obs} total cells")
                # Random sampling
                indices = np.random.choice(adata.n_obs, size=n_cells, replace=False)
                adata = adata[indices].copy()
                rdata = rdata[indices].copy()
                print(f"After sampling, data shape: {adata.shape}")
            else:
                print("Using full dataset (sample_percent >= 1.0)")
        
        # Convert to dense array for preprocessing
        print("Converting to dense array...")
        if sparse.issparse(adata.X):
            adata.X = adata.X.toarray()
        if sparse.issparse(rdata.X):
            rdata.X = rdata.X.toarray()
            
        # Convert to float32 to reduce memory usage
        adata.X = adata.X.astype('float32')
        rdata.X = rdata.X.astype('float32')
        
        # Basic preprocessing
        print("Performing basic preprocessing...")
        # Replace inf/-inf with max/min finite values
        adata.X = np.nan_to_num(adata.X, nan=0)
        rdata.X = np.nan_to_num(rdata.X, nan=0)
        
        # Clip extreme values
        print("Clipping extreme values...")
        upper = np.percentile(adata.X[adata.X > 0], 99)
        adata.X = np.clip(adata.X, 0, upper)
        rdata.X = np.clip(rdata.X, 0, upper)
        
        # Log transform
        print("Log transforming data...")
        sc.pp.log1p(adata)
        sc.pp.log1p(rdata)
        
        # Calculate highly variable genes manually
        print("Calculating highly variable genes...")
        means = np.mean(adata.X, axis=0)
        vars = np.var(adata.X, axis=0)
        dispersion = vars / means
        
        # Select top 2000 genes by dispersion
        n_top = 2000
        dispersion[np.isnan(dispersion)] = 0  # Handle any remaining NaNs
        top_genes_idx = np.argsort(dispersion)[-n_top:]
        
        # Create highly_variable column in var
        adata.var['highly_variable'] = False
        adata.var.iloc[top_genes_idx, adata.var.columns.get_loc('highly_variable')] = True
        
        rdata.var['highly_variable'] = adata.var['highly_variable'].copy()
        
        print(f"Number of highly variable genes: {sum(adata.var['highly_variable'])}")
        
        # Run PCA with reduced number of components
        print("Running PCA...")
        sc.pp.pca(adata, n_comps=50)
        
        # Compute neighborhood graph
        print("Computing neighborhood graph...")
        sc.pp.neighbors(adata, n_neighbors=10, n_pcs=50)
        
        # Run UMAP
        print("Running UMAP...")
        sc.tl.umap(adata)
        
        # Run Leiden clustering
        print("Running Leiden clustering...")
        sc.tl.leiden(adata, resolution=0.5)
        
        # Calculate pseudotime using diffusion pseudotime
        print("Calculating diffusion pseudotime...")
        sc.tl.diffmap(adata)
        adata.uns['iroot'] = np.argmin(adata.obsm['X_diffmap'][:, 0])
        sc.tl.dpt(adata)
        
        # Get marker genes
        print("Getting marker genes...")
        markers_uniq = get_DEG_uniq(rdata, adata, group_key='leiden', power=11, ratio=0.2, p_threshold=0.01, q_threshold=0.05)
        
        # Extract UMAP coordinates and metadata
        print("Extracting visualization data...")
        umap_coords = adata.obsm['X_umap']
        cell_types = adata.obs['cell_type'].tolist() if 'cell_type' in adata.obs else ['Unknown'] * len(adata)
        clusters = adata.obs['leiden'].tolist()
        pseudotime = adata.obs['dpt_pseudotime'].tolist()
        
        # Create response data structure
        response_data = {
            "umap_coordinates": {
                "x": umap_coords[:, 0].tolist(),
                "y": umap_coords[:, 1].tolist()
            },
            "cell_types": cell_types,
            "cell_clusters": clusters,
            "pseudotime": pseudotime,
            "metadata": {
                "unique_cell_types": sorted(list(set(cell_types))),
                "unique_clusters": sorted(list(set(clusters))),
                "total_cells": len(adata),
                "sampled_cells": len(adata) if sample_percent is None else int(len(adata) * sample_percent),
                "root_cell_index": int(adata.uns['iroot'])
            }
        }
        
        # Save results
        output_dir = 'output'
        os.makedirs(output_dir, exist_ok=True)
        output_file = os.path.join(output_dir, 'umap_data.json')
        with open(output_file, 'w') as f:
            json.dump(response_data, f)
        
        print("Analysis completed successfully")
        print(f"Response data shape: {len(response_data['umap_coordinates']['x'])} cells")
        return response_data
        
    except Exception as e:
        print(f"Error in run_deaplog_analysis: {str(e)}")
        raise

def main():
    parser = argparse.ArgumentParser(description='Run DEAPLOG analysis')
    parser.add_argument('--sample_percent', type=float, default=0.01,
                      help='Percentage of cells to sample for analysis (default: 0.01 for 1% of data)')
    parser.add_argument('--step', type=int, default=0,
                      help='Analysis step to run')
    parser.add_argument('--data_path', type=str, required=True,
                      help='Path to the AnnData file')
    
    args = parser.parse_args()
    
    print(f"Loading data from {args.data_path}")
    try:
        # Load data
        adata = sc.read_h5ad(args.data_path)
        print(f"Data loaded successfully. Shape: {adata.shape}")
        rdata = adata.copy()  # For normalized counts
        print("Data copied successfully")
        
        # Run analysis
        results = run_deaplog_analysis(rdata, adata, args.sample_percent)
        
        # Print results as JSON
        print(json.dumps(results))
    except Exception as e:
        print(f"Error loading or processing data: {str(e)}")
        raise

if __name__ == "__main__":
    main()
