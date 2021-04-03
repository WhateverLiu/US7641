


x = runif(20000, min = -5, max = 5)
y = runif(20000, min = -5, max = 5)
tmp = 
x ^ 2 + y ^ 2 <= 2 ^ 2 & x ^ 2 + y ^ 2 >= 1.8 ^ 2 |
x ^ 2 + y ^ 2 <= 3 ^ 2 & x ^ 2 + y ^ 2 >= 2.7 ^ 2 |
x ^ 2 + y ^ 2 <= 4 ^ 2 & x ^ 2 + y ^ 2 >= 3.7 ^ 2 |
(x - 0.7) ^ 2 + (y - 0.7) ^ 2 <= 0.3 ^ 2 |
(x + 0.7) ^ 2 + (y - 0.7) ^ 2 <= 0.3 ^ 2 |
(x - 0) ^ 2 + (y + 1) ^ 2 <= 0.3 ^ 2 |
(x - 0.1) ^ 2 + (y + 1) ^ 2 <= 0.3 ^ 2 |
(x + 0.1) ^ 2 + (y + 1) ^ 2 <= 0.3 ^ 2 |
(x - 0.2) ^ 2 + (y + 1) ^ 2 <= 0.3 ^ 2 |
(x + 0.2) ^ 2 + (y + 1) ^ 2 <= 0.3 ^ 2 |
(x - 0.3) ^ 2 + (y + 1) ^ 2 <= 0.3 ^ 2 |
(x + 0.3) ^ 2 + (y + 1) ^ 2 <= 0.3 ^ 2
x = x[tmp]
y = y[tmp]
plot(x, y, cex = 0.2, pch = 16)
circleface = data.frame(x = x, y = y)
save(circleface, file = "data/circleface.Rdata")








if(F)
{
  
  
  N = 300
  dat = list()
  r = 5
  theta = seq(0, 2 * pi, len = N)
  dat[[1]] = r * cbind(cos(theta), sin(theta))
  
  
  r = 4.5
  theta = seq(0, 2 * pi, len = N)
  dat[[2]] = r * cbind(cos(theta), sin(theta))
  
  
  r = 4
  theta = seq(0, 2 * pi, len = N)
  dat[[3]] = r * cbind(cos(theta), sin(theta))
  
  
  
  # x = 3 * cos(pi / 6)
  # y = 3 * sin(pi / 6)
  x = 2 * cos(pi / 6)
  y = 1 * sin(pi / 6)
  r = 0.75
  tmp = cbind(runif(10000, min = x - r, max = x + r), runif(10000, min = y - r, max = y + r))
  tmp = tmp[(tmp[,1] - x) ^ 2 + (tmp[,2] - y) ^ 2 <= r ^ 2, ]
  tmp = tmp[sample(nrow(tmp), N), ]
  dat[[4]] = tmp
  
  
  x = -2.5 * cos(pi / 6)
  y = 2.5 * sin(pi / 6)
  r = 0.5
  tmp = cbind(runif(10000, min = x - r, max = x + r), runif(10000, min = y - r, max = y + r))
  tmp = tmp[(tmp[,1] - x) ^ 2 + (tmp[,2] - y) ^ 2 <= r ^ 2, ]
  tmp = tmp[sample(nrow(tmp), N), ]
  dat[[5]] = tmp
  
  
  tmp = cbind(runif(N, min = -1.5, max = 1.5), runif(N, min = -2.5, max = -2.25))
  dat[[6]] = tmp
  
  
  dat = as.data.frame(lapply(1:ncol(dat[[1]]), function(i) unlist(lapply(dat, function(x) x[, i]))))
  colnames(dat) = c("x", "y")
  plot(dat)
  # dat = dat[sample(nrow(dat)), ]
  
  
  # Test if diffusion maps work perfectly.
  if(T)
  {
    
    
    X = dat
    tX = t(X)
    
    
    dmat = as.matrix(dist(X))
    # eps = max(dmat) ^ 2; kernelDmat = exp(-dmat ^ 2 / eps)
    dmat[order(dmat)] = sort(rnorm(length(dmat))); dmat = dmat - min(dmat); dmat[upper.tri(dmat)] = dmat[lower.tri(dmat)]; kernelDmat = dmat
    
    
    s = sqrt(colSums(kernelDmat))
    s = s %*% t(s)
    kernelDmat = kernelDmat / s
    
    
    P = kernelDmat / colSums(kernelDmat) # Symmetric
    M = P
    
    
    Rcpp::sourceCpp("src/speedup.cpp", verbose = T)
    rst = list()
    for(i in 1:1000)
    {
      lastM = M
      M = matmul(M, M, maxCore = 15) 
      M = M / rowSums(M)
      meanErr = mean(abs(M - lastM))
      cat(i, ", err = ", meanErr, "    ")
      rst[[length(rst) + 1]] = M
      if(meanErr < 1e-17) break 
    }
    stationaryP = rst[[length(rst)]]
    plot(stationaryP[1,])
    
    
    
    M1 = P; M2 = rst[[1]]; M4 = rst[[2]]
    
    
    M1svd = svd(t(M1 / stationaryP ^ 0.5))
    plot(M1svd$u[, 1:2])
    rgl::plot3d(M1svd$u[, 1:3])
    
    
    
  }
  
}















