#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 11:28:10 2019

@author: Nick Wawee w/ adaptation by Warren Sink

This script will load correlation columns from a given directory and perform kmeans clustering
"""

import seaborn as sns
from sklearn.cluster import KMeans
from sklearn.cluster import AgglomerativeClustering
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

adpd=pd.read_csv("/Users/Warren.Sink/Desktop/OneDrive-Van-Andel-Institute/schoeller_all_dist_pheno_delta.csv",header=0)
adpfc=pd.read_csv("/Users/Warren.Sink/Desktop/OneDrive-Van-Andel-Institute/schoeller_all_dist_pheno_fc.csv",header=0)

#delta cor
adpd_val = adpd.iloc[:,3:353] 
adpd_corr=adpd_val.corr(method='pearson')

#fc cor
adpfc_val = adpfc.iloc[:,3:353]
adpfc_corr=adpfc_val.corr(method='pearson')

#1st heatmap
heatmap=sns.clustermap(adpd_corr.abs(), metric="euclidean", method="ward")
heatmap.savefig("/Users/Warren.Sink/Desktop/OneDrive-Van-Andel-Institute/heatmap_adpd_py.png",dpi=600)

heatmap=sns.clustermap(adpfc_corr.abs(), metric="euclidean", method="ward")
heatmap.savefig("/Users/Warren.Sink/Desktop/OneDrive-Van-Andel-Institute/heatmap_adpfc_py.png",dpi=600)

#use if hierarchial does not work, finding elbow
numclust=20
ssd=[]
for i in range(1,numclust):
    kmeans=KMeans(n_clusters=i, random_state=123)
    ssd.append(kmeans.fit(adpd_corr).inertia_)

plt.scatter(range(1,numclust), ssd)
plt.xlabel('k')
plt.ylabel('Total Within Variation Sum of Squares')
plt.xticks(range(1,numclust,2))
plt.savefig("elbow_adpd.png")

#fitting kmeans
numclust=3
kmeans=KMeans(n_clusters=numclust, random_state=123)
res=kmeans.fit(adpd_corr)
centriods=res.cluster_centers_
ssd=res.inertia_

#hierarchial clustering on centriods
cluster = AgglomerativeClustering(n_clusters=numclust, affinity='euclidean', linkage='ward')
clusters=cluster.fit_predict(centriods.transpose())
clustersdf=pd.DataFrame(clusters)
clustersdf.index=adpd_corr.index

#sorting rows
alldf=pd.concat([clustersdf,adpd_corr], axis=1)
alldf=alldf.sort_values(by=[0])
alldf=alldf.drop([0], axis=1)

#sorting columns
alldf=pd.concat([clustersdf.transpose(),alldf])
alldf=alldf.sort_values(by=[0], axis=1)
alldf=alldf.drop([0], axis=0)

#plotting heatmap
heatmap=sns.clustermap(alldf.abs())
heatmap.savefig("/Users/Warren.Sink/Desktop/OneDrive-Van-Andel-Institute/heatmap_adpd_elbowcorrected_5.png",dpi=600)
#clustersdf=clustersdf.sort_values(by=[0])
