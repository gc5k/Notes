n=100
Chr=20
Me=(1:Chr)*10
Mat=matrix(0, Chr+1, Chr+1)
D_mat=n^2/Me
diag(Mat)=c(D_mat,0)
Mat=Mat+n

Mat_I=solve(Mat)

Mat_IE=matrix(0, Chr+1, Chr+1)
Mat_IE[1:Chr, ncol(Mat_IE)]=-1/D_mat
diag(Mat_IE)=c(1/D_mat, 1/n+sum(Me)/n^2)
Mat_IE[nrow(Mat_IE), 1:(ncol(Mat_IE)-1)]=-1/D_mat

diag(Mat_I)
diag(Mat_IE)

Mat_I[1:Chr, ncol(Mat_I)]
Mat_IE[1:Chr, ncol(Mat_IE)]
Mat_I[nrow(Mat_I), ]
Mat_IE[nrow(Mat_IE), ]
