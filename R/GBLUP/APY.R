#simu
Gc=GenerateGenoDprimeRcpp(frq, Dp, N)
Gn=GenerateGenoDprimeRcpp(frq, Dp, N)
plot(colMeans(Gc)/2, frq)

sGc=apply(Gc, 2, scale)
sGn=apply(Gn, 2, scale)

grm_c=sGc%*%t(sGc)/M
GRM_C=CorMatrixRcpp(sGc)
GRM_N=CorMatrixRcpp(sGn)
GRM_CN=CorMatrix2Rcpp(sGc, sGn)


GRM_C=matrix(c(0.81, 0, 0, 0.81), 2, 2, byrow = T)
GRM_N=matrix(0, 3, 3)
diag(GRM_N)=c(.01, 1.61, 1.61)
GRM_N[2,3]=GRM_N[3,2]=-1.6
GRM_CN=matrix(0, 2, 3)
GRM_CN[,2]=0.8
GRM_CN[,3]=-0.8

IC=solve(GRM_C)
Pcn=IC %*% GRM_CN
Pnc=t(GRM_CN) %*% IC
Pcn_Gcn=t(Pcn) %*% GRM_CN
Mnn=diag(diag(GRM_N)-diag(Pcn_Gcn), nrow=nrow(GRM_N), ncol=nrow(GRM_N))
IMnn=solve(Mnn)

v1_1=rbind(IC, matrix(0, nrow=nrow(GRM_N), ncol=ncol(GRM_C)))
v1_2=rbind(-1*Pcn%*%IMnn, IMnn)
V1=cbind(v1_1, v1_2)

v2_1=rbind(diag(1, nrow=nrow(GRM_C), ncol=nrow(GRM_C)), -1*Pnc)
v2_2=rbind(matrix(0, nrow=nrow(GRM_C), ncol=nrow(GRM_N)), diag(1, nrow(GRM_N), ncol(GRM_N)))
V2=cbind(v2_1, v2_2)

IV=V1 %*% V2
IVI=solve(IV)
