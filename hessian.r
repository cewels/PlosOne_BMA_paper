#find the hessian matrix

library(MASS)
library(corpcor)
z2<-read.table("paras_LCA.txt", header=T)
para_l<-z2[, -c(1:4)] 

min_inv_hes_L<-cov(para_l)

det_l<-determinant(min_inv_hes_L)

z1<-read.table("paras_GOM_gamma.txt", header=T)
para_g<-z1[, -c(1:4)]

min_inv_hes_g<-cov(para_g)
det_g<-determinant(min_inv_hes_g)
