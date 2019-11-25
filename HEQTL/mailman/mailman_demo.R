############ 出于演示目的模拟用户构造算法输入 ############

sigma = c(0, 1, 2)  # 矩阵A的元素只能从sigma里面挑，并且sigma里的元素必须为整数

m = 4  # 矩阵A的行数

n = length(sigma)^m  # 矩阵A的列数；为了便于演示，假设n正好等于|sigma|^m

# 随机生成一个m行n列的矩阵A，元素全都来自sigma。这就是算法的输入矩阵A。
A = matrix(sample(sigma, m * n, replace = TRUE), nrow = m, ncol = n)

# 随机生成两个长度为n的实数向量。
x1 = runif(n, -5, 5)  # -5, 5是随便取的，可以是任意值
x2 = runif(n, -5, 5)

########### 进入算法预处理阶段，将A分解为U*P ###############

buildU <- function(m, sigma) {
	S = length(sigma)
	if (m <= 1) {
		return (t(matrix(sigma)))
	} else {
		Usub = buildU(m - 1, sigma)
		up = t(matrix(rep(sigma, each = ncol(Usub))))
		down = do.call(cbind, replicate(S, Usub, simplify = FALSE))
		return (rbind(up, down))
	}
}

U = buildU(m, sigma)

# P只是演示用的，便于理解，实际计算的时候用的是Pcompact
P = matrix(0, nrow = n, ncol = n)
Pcompact = array(list(list()), n)

for (i in 1:n) {
  for (j in 1:n) {
    if (all(U[,i]==A[,j])) {
      P[i,j] = 1
      Pcompact[[i]]=append(Pcompact[[i]],j)
    }
  }
}

############ 计算 A*x1 和 A*x2 #############

mailman_product <- function(x) {
  # Px = P * x
  Px = c()
  for (i in 1:n) {
    Px = c(Px, sum(x[unlist(Pcompact[[i]])]))
  }
  
  Ax = matrix(0, nrow = m, ncol = 1)
  z = matrix(Px, nrow = n, ncol = 1)
  for (i in 1:m) {
    n1 = length(z) / length(sigma)
    z1 = matrix(0, nrow = n1, ncol = 1)
    for (j in 1:length(sigma)) {
      z2 = z[((j-1)*n1+1):(j*n1)]
      Ax[i] = Ax[i] + sigma[j] * sum(z2)
      z1 = z1 + z2
    }
    z = z1
  }
  return (Ax)
}

Ax1 = mailman_product(x1)
print(paste("A %*% x1 == Ax1: ", all.equal(A %*% x1, Ax1)))

Ax2 = mailman_product(x2)
print(paste("A %*% x2 == Ax2: ", all.equal(A %*% x2, Ax2)))