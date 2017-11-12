# JRmGRN
Joint reconstruction of multiple gene regulatory networks with common hub genes

Abstract

Motivation: Joint reconstruction of multiple gene regulatory networks (GRNs) using gene expression data from multiple tissues/conditions is very important for understanding common and tissue- or condition-specific regulation. However, there are currently no computational models and methods available for directly constructing such multiple GRNs that not only share some common hub genes but also possess tissue- or condition-specific regulatory edges.  

Results: In this paper, we proposed a new graphic Gaussian model for joint reconstruction of multiple gene regulatory networks (JRmGRN), which highlighted hub genes, using gene expression data from several tissues/conditions. Under the framework of Gaussian graphical model, JRmGRN method constructs the GRNs through maximizing a penalized log likelihood function. We formulated it as a convex optimization problem, and then solved it with an alternating direction method of multipliers (ADMM) algorithm. The performance of JRmGRN was first evaluated with synthetic data and the results showed that JRmGRN outperformed several other methods for reconstruction of GRNs. We also applied our method to real Arabidopsis thaliana RNA-seq data from two light regime conditions, and both common hub genes and some conditions-specific hub genes were identified.  
