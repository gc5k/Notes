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

buildP <- function(u, a, n) { #zzx's original code
  pcompact = array(list(list()), n)
  for (i in 1:n) {
    for (j in 1:n) {
      if (all(u[,i]==a[,j])) {
        pcompact[[i]]=append(pcompact[[i]],j)
      }
    }
  }
  return (pcompact)
}

buildPF <- function(a, bs) { #modified zzx's code by without using U
  pcompactF = array(list(list()), ncol(a))
  for(i in 1:ncol(a)) {
    ui=sum(bs*a[,i])+1
    pcompactF[[ui]]=append(pcompactF[[ui]], i)
  }
  return(pcompactF)
}

buildPF2 <- function(a, bs, n) { #further modify by directly adding the last 0-strings
  pcompactF = array(list(list()), ncol(a))
  for(i in 1:n) {
    ui=sum(bs*a[,i])+1
    pcompactF[[ui]]=append(pcompactF[[ui]], i)
  }
  if(ncol(a) > n) {
    pcompactF[[1]]=append(pcompactF[[1]], seq(n+1, ncol(a)))
  }
  return(pcompactF)
}

mailman_product2 <- function(m, n, pcompact, sigma, x) {
  # Px = P * x
  Px = c()
  for (i in 1:n) {
    Px = c(Px, sum(x[unlist(pcompact[[i]])]))
  }
  
  Ax = matrix(0, nrow = m, ncol = 1)
  z = matrix(Px, nrow = n, ncol = 1)
  for (i in 1:m) {
    n1 = length(z) / length(sigma)
    z1 = matrix(0, nrow = n1, ncol = 1)
    for (j in 1:length(sigma)) {
      z2 = z[((j-1)*n1+1):(j*n1)]
      z1 = z1 + z2
      if (sigma[j]!=0) {
        Ax[i] = Ax[i] + sigma[j] * sum(z2)
      }
    }
    z = z1
  }
  return (Ax)
}

mailman_product <- function(m, n, pcompact, sigma, x) {
  # Px = P * x
  Px = c()
  for (i in 1:n) {
    Px = c(Px, sum(x[unlist(pcompact[[i]])]))
  }

  Ax = matrix(0, nrow = m, ncol = 1)
  z = matrix(Px, nrow = n, ncol = 1)
  for (i in 1:m) {
    n1 = length(z) / length(sigma)
    z1 = matrix(0, nrow = n1, ncol = 1)
    for (j in 1:length(sigma)) {
      z2 = z[((j-1)*n1+1):(j*n1)]
      z1 = z1 + z2
      Ax[i] = Ax[i] + sigma[j] * sum(z2)
    }
    z = z1
  }
  return (Ax)
}
