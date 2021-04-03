

load("data/circleface.Rdata")
X = circleface
tX = t(X)


# For every k (number of clusters), perform K-means 30 times. Collect all results.
if(T)
{
  set.seed(123)
  rst = list()
  Ks = 3:100
  for(k in Ks)
  {
    cat(k, "")
    rst[[length(rst) + 1]] = list()
    for(i in 1:30)
    {
      rst[[length(rst)]][[i]] = GMKMcharlie::KM(X = tX, centroid = tX[, sample(1:nrow(X), k)], verbose = F, maxCore = 7)
    }
  } 
}


# For every k, select the clustering having the least sum of squares.
if(T)
{
  
  
  rstSelected = lapply(rst, function(x)
  {
    tmp = which.min(unlist(lapply(x, function(u)
    {
      sum(unlist(lapply(u, function(v) sum(v$member2centroidDistance ^ 2))))
    })))
    x[[tmp]]
  })
  
  
  # Compute BIC.
  clusterings = lapply(rstSelected, function(x)
  {
    lapply(x, function(y) X[y$clusterMember, ])
  })
  
  
  clusteringsBIC = unlist(lapply(clusterings, function(x) kmeansBIC(x)))
  
  
  pdf("figure/kmeansBICcircleface.pdf", width = 8, height = 8 * 0.618)
  par(mar = c(4.1, 5, 0, 0), family = "serif")
  plot(x = Ks, y = clusteringsBIC, lwd = 2, col = "darkblue", bty = "L", xlab = "K", ylab = "BIC", cex.lab = 2, cex.axis = 1.5, type = "l")
  tmp = Ks[which.min(clusteringsBIC)]
  lines(x = tmp, y = 1e100, type = "h", lty = 2)
  legend("right", bty = "n", legend = paste0("Optimal K = ", tmp), cex = 2)
  legend("top", bty = "n", legend = "Circleface data", cex = 2)
  dev.off()
  
  
  bestClustering = clusterings[[which.min(clusteringsBIC)]]
  palette = randomcoloR::distinctColorPalette(length(bestClustering))
  tmp = as.data.frame(lapply(1:ncol(bestClustering[[1]]), function(i) unlist(lapply(bestClustering, function(x) x[[i]]))))
  cols = unlist(mapply(function(x, y) rep(x, nrow(y)), palette, bestClustering, SIMPLIFY = T))
  pdf("figure/K52circlefaceClustering.pdf", width = 8, height = 8 * 0.618)
  par(mar = c(4.1, 5, 0, 0), family = "serif")
  plot(tmp, col = cols, bty = "L", cex.lab = 2, cex.axis = 1.5, xlab = "x", ylab = "y")
  dev.off()
  
  
  
}




# For every g (number of gaussian kernels), perform GMM 30 times. Collect all results.
if(T)
{
  
  
  set.seed(123)
  rst = list()
  Ks = 3:100
  for(k in Ks)
  {
    cat(k, "")
    rst[[length(rst) + 1]] = list()
    for(i in 1:30)
    {
      mu = tX[, sample(1:nrow(X), k)]
      sig = as.numeric(diag(diag(cov(X)) / ncol(X), ncol = ncol(X)))
      sig = matrix(rep(sig, k), ncol = k)
      rst[[length(rst)]][[i]] = GMKMcharlie::GM(X = tX, alpha = rep(1 / k, k), mu = mu, sigma = sig, verbose = F, maxCore = 7, eigenRatioLim = 1000, maxIter = 10000, loglikehoodConverge = T, loglikehoodConvergeBlock = 20, convergenceEPS = 0.01)
    }
  }
  save(rst, file = "data/circlefaceGMrst.Rdata")
}


# For every g, select the GMM having the least negative log-likihood.
source("R/rfuns.r")
if(T)
{
  
  
  rstSelected = lapply(rst, function(x)
  {
    tmp = which.min(unlist(lapply(x, function(u)
    {
      -sum(log(u$fitted))
    })))
    x[[tmp]]
  })
  
  
  clusterings = rstSelected
  clusteringsBIC = unlist(lapply(clusterings, function(x) gmBIC(nrow(X), ncol(X), length(x$alpha), x$fitted) ))
  
  
  pdf("figure/gmBICcircleface.pdf", width = 8, height = 8 * 0.618)
  par(mar = c(4.1, 5, 0, 0), family = "serif")
  plot(x = Ks, y = clusteringsBIC, lwd = 2, col = "darkblue", bty = "L", xlab = "K", ylab = "BIC", cex.lab = 2, cex.axis = 1.5, type = "l")
  tmp = Ks[which.min(clusteringsBIC)]
  lines(x = tmp, y = 1e100, type = "h", lty = 2)
  legend("right", bty = "n", legend = paste0("Optimal K = ", tmp), cex = 2)
  legend("top", bty = "n", legend = "Circleface data", cex = 2)
  dev.off()
  
  
  bestClustering = clusterings[[which.min(clusteringsBIC)]]
  bestClustering = lapply(bestClustering$clusterMember, function(x) X[x, ]    )
  palette = randomcoloR::distinctColorPalette(length(bestClustering))
  tmp = as.data.frame(lapply(1:ncol(bestClustering[[1]]), function(i) unlist(lapply(bestClustering, function(x) x[[i]]))))
  cols = unlist(mapply(function(x, y) rep(x, nrow(y)), palette, bestClustering, SIMPLIFY = T))
  pdf("figure/g48circlefaceClustering.pdf", width = 8, height = 8 * 0.618)
  par(mar = c(4.1, 5, 0, 0), family = "serif")
  plot(tmp, col = cols, bty = "L", cex.lab = 2, cex.axis = 1.5, xlab = "x", ylab = "y")
  dev.off()
  
  
}




# Extend X to high dimensional space.
if(T)
{
  load("data/circleface.Rdata")
  tX = t(circleface)
  distmat = as.matrix(dist(circleface))
  dimnames(distmat) = NULL
  NN = 50L
  tX = apply(distmat, 2, function(x) as.numeric(tX[, order(x)[1:NN]]))
  X = t(tX)
}




# PCA
if(F)
{
  
  
  Xcentered = apply(X, 2, function(x) x - mean(x))
  tmp = svd(Xcentered)
  pdf("figure/highDcirclefacePCA.pdf", width = 8, height = 8 * 0.618)
  par(mar = c(4.1, 5, 0, 0), family = "serif")
  d = rev(cumsum(rev(tmp$d ^ 2)))
  d = d / d[1]
  plot(d, xlab = "Eigen index", ylab = "Cumulative variance", cex.lab = 2, cex.axis = 1.5, bty = "L")
  dev.off()
  
  
  pdf("figure/highDcirclefacePCAplotInEigenSpace.pdf", width = 8, height = 8 * 0.618)
  par(mar = c(4.1, 5, 0, 0), family = "serif")
  plot(tmp$u[, 1:2] %*% diag(tmp$d[1:2]), xlab = "x", ylab = "y", cex.lab = 2, cex.axis = 1.5, bty = "L")
  dev.off()
  
  
  
}




# ICA
if(T)
{
  
  
  # fastICA(X, n.comp, alg.typ = c("parallel","deflation"), fun = c("logcosh","exp"), alpha = 1.0, method = c("R","C"), row.norm = FALSE, maxit = 200, tol = 1e-04, verbose = FALSE, w.init = NULL)
  
  
  # tmp2 = fastICA::fastICA(X, n.comp = 100, alg.typ = c("parallel"), fun = c("logcosh"), alpha = 1.0, method = c("C"), row.norm = FALSE, maxit = 200, tol = 1e-04, verbose = FALSE, w.init = NULL)
  
  
  set.seed(42)
  tmp = fastICA::fastICA(X, n.comp = 2, alg.typ = c("deflation"), fun = c("logcosh"), alpha = 1.0, method = c("C"), row.norm = FALSE, maxit = 200, tol = 1e-04, verbose = FALSE, w.init = NULL)
  pdf("figure/highDcirclefaceICA2components.pdf", width = 8, height = 8 * 0.618)
  par(mar = c(4.1, 5, 0, 0), family = "serif")
  tmpkurto = round(kurto(tmp$S), 2)
  plot(tmp$S, xlab = paste0("x, kurtosis = ", tmpkurto[1]), ylab = paste0("y, kurtosis = ", tmpkurto[2]), cex.lab = 2, cex.axis = 1.5, bty = "L")
  dev.off()
  
  
  set.seed(42)
  tmp = fastICA::fastICA(X, n.comp = 4, alg.typ = c("deflation"), fun = c("logcosh"), alpha = 1.0, method = c("C"), row.norm = FALSE, maxit = 200, tol = 1e-04, verbose = FALSE, w.init = NULL)
  png("figure/highDcirclefaceICA4components.png", width = 8, height = 8 * 0.618, res = 120, unit = "in")
  par(mar = c(1, 1, 0, 0), family = "serif")
  plot(as.data.frame(tmp$S), cex.lab = 2, cex.axis = 1.5, bty = "L")
  dev.off()
  
  
  set.seed(42)
  icaRst = list()
  for(i in 2:ncol(X))
  {
    cat(i, "")
    icaRst[[length(icaRst) + 1]] = fastICA::fastICA(X, n.comp = i, alg.typ = c("deflation"), fun = c("logcosh"), alpha = 1.0, method = c("C"), row.norm = FALSE, maxit = 200, tol = 1e-04, verbose = FALSE, w.init = NULL)
  }
  for(i in 2:length(icaRst)) icaRst[[i]]$X = NULL
  # save(icaRst, file = "data/icaRst.Rdata")
  kurtoAll = lapply(icaRst, function(x) kurto(x$S))
  abskurtoAll = as.numeric(lapply(kurtoAll, function(x) mean(abs(x))))
  
  
  pdf("figure/highDcirclefaceICAmeanExcessKurto100.pdf", width = 8, height = 8 * 0.618)
  par(mar = c(4.1, 5, 0, 0), family = "serif")
  plot(abskurtoAll, xlab = "Number of ICA components", ylab = "Mean excess kurtosis", cex.lab = 2, cex.axis = 1.5, bty = "L", type = "l", col = "darkblue")
  dev.off()
  
  
  i = 99L
  tmp = icaRst[[99]]$S[, order(-abs(kurtoAll[[99]]))[1:3]]
  # png("figure/highDcirclefaceICA100componentsLeastGaussian.png", width = 8, height = 8 * 0.618, res = 120, unit = "in")
  breaks = seq(min(as.numeric(tmp)), max(as.numeric(tmp)), len = 100)
  pdf("figure/highDcirclefaceICA100componentsLeastGaussian.pdf", width = 8, height = 8 * 0.618)
  par(mar = c(4.1, 5, 0, 0), family = "serif", mfrow = c(3, 1))
  hist(tmp[,1], breaks = breaks, main = "", xlab = "", ylab = "Frequency", col = scales::alpha("blue", 0.5),cex.lab = 2, cex.axis =1.5, border = NA)
  legend("top", legend = paste0("Least Gaussian, excess kurtosis = ", round(kurto(tmp[,1]), 2)), bty = "n", cex = 2)
  hist(tmp[,2], breaks = breaks, main = "", xlab = "", ylab = "", col = scales::alpha("blue", 0.5),cex.lab = 2, cex.axis =1.5, border = NA)
  legend("top", legend = paste0("Second least Gaussian, excess kurtosis = ", round(kurto(tmp[,2]), 2)), bty = "n", cex = 2)
  hist(tmp[,3], breaks = breaks, main = "", xlab = "x", ylab = "", col = scales::alpha("blue", 0.5),cex.lab = 2, cex.axis =1.5, border = NA)
  legend("top", legend = paste0("Third least Gaussian, excess kurtosis = ", round(kurto(tmp[,3]), 2)), bty = "n", cex = 2)
  dev.off()
  
}




# Random projection.
# (1 - eps)||u - v||^2 < ||p(u) - p(v)||^2 < (1 + eps)||u - v||^2
findNcomp = function(eps, Ndata) { 4 * log(Ndata) / (eps ^ 2 / 2 - eps ^ 3 / 3) }
findNcomp(1, nrow(X)); findNcomp(0.5, nrow(X)); findNcomp(0.2, nrow(X)); findNcomp(0.1, nrow(X))
# 193; 386; 1857; 6897;
if(T)
{
  
  
  set.seed(42)
  
  
  XorignalPairwiseD = dist(X)
  Ncomponent = c(2L, 10L, 50L, 100L)
  
  
  # Gaussian random projection
  if(T)
  {
    XlowDs = lapply(Ncomponent, function(x) 
    {
      lapply(1:20, function(i)
      {
        R = matrix(rnorm(x * ncol(X)) / x ^ 0.5, ncol = x)
        X %*% R
      })
    })
  }
  
  
  reconsErr = lapply(XlowDs, function(x) 
  {
    cat(".")
    tmp = as.data.frame(lapply(x, function(u) as.numeric(dist(u) / XorignalPairwiseD)))
    rowMeans(tmp)
  })
  
  
  pdf("figure/highDcirclefaceRPreconErrHist.pdf", width = 8, height = 8 * 0.618)
  par(mar = c(4.1, 5, 0, 0), family = "serif", mfrow = c(2, 2))
  tmp = lapply(reconsErr, function(x) x ^ 2)
  hist(tmp[[1]], breaks = 100, main = "", xlim = c(0, 2), cex.lab = 2, cex.axis = 1.5, col = scales::alpha("darkblue", 0.3), border = NA, xlab = "", ylab = "Frequency")
  legend("topright", bty = "n", cex = 2, legend = "N(component) = 2")
  hist(tmp[[2]], breaks = 100, main = "", xlim = c(0, 2), cex.lab = 2, cex.axis = 1.5, col = scales::alpha("darkblue", 0.3), border = NA, xlab = "", ylab = "")
  legend("topright", bty = "n", cex = 2, legend = "10")
  hist(tmp[[3]], breaks = 100, main = "", xlim = c(0, 2), cex.lab = 2, cex.axis = 1.5, col = scales::alpha("darkblue", 0.3), border = NA, xlab = "Projected / original", ylab = "")
  legend("topright", bty = "n", cex = 2, legend = "50")
  hist(tmp[[4]], breaks = 100, main = "", xlim = c(0, 2), cex.lab = 2, cex.axis = 1.5, col = scales::alpha("darkblue", 0.3), border = NA, xlab = "", ylab = "")
  legend("topright", bty = "n", cex = 2, legend = "100")
  dev.off()
  
  
  # Average projections.
  if(T)
  {
    set.seed(42)
    XlowDs = lapply(Ncomponent, function(x) 
    {
      tmp = lapply(1:30, function(i)
      {
        R = matrix(rnorm(x * ncol(X)) / x ^ 0.5, ncol = x)
        X %*% R
      })
      rst = tmp[[1]]
      for(i in 2:30) rst = rst + tmp[[i]]
      rst / 30 ^ 0.5
    })
  }
  
  
  reconsErrAvgProj = lapply(XlowDs, function(x) 
  {
    cat(".")
    as.numeric(dist(x) / XorignalPairwiseD)
  })
  
  
  pdf("figure/highDcirclefaceAvgRPReconErrHist.pdf", width = 8, height = 8 * 0.618)
  par(mar = c(4.1, 5, 0, 0), family = "serif", mfrow = c(2, 2))
  tmp = lapply(reconsErrAvgProj, function(x) x)
  hist(tmp[[1]], breaks = 100, main = "", xlim = c(0, 2), cex.lab = 2, cex.axis = 1.5, col = scales::alpha("darkblue", 0.3), border = NA, xlab = "", ylab = "Frequency")
  legend("topright", bty = "n", cex = 2, legend = "N(component) = 2")
  hist(tmp[[2]], breaks = 100, main = "", xlim = c(0, 2), cex.lab = 2, cex.axis = 1.5, col = scales::alpha("darkblue", 0.3), border = NA, xlab = "", ylab = "")
  legend("topright", bty = "n", cex = 2, legend = "10")
  hist(tmp[[3]], breaks = 100, main = "", xlim = c(0, 2), cex.lab = 2, cex.axis = 1.5, col = scales::alpha("darkblue", 0.3), border = NA, xlab = "Projected / original", ylab = "")
  legend("topright", bty = "n", cex = 2, legend = "50")
  hist(tmp[[4]], breaks = 100, main = "", xlim = c(0, 2), cex.lab = 2, cex.axis = 1.5, col = scales::alpha("darkblue", 0.3), border = NA, xlab = "", ylab = "")
  legend("topright", bty = "n", cex = 2, legend = "100")
  dev.off()
  
  
  
}




# Diffusion map. test.
if(F)
{
  
  
  load("data/circleface.Rdata")
  X = circleface
  tX = t(X)
  
  
  dmat = as.matrix(dist(X))
  eps = max(dmat) ^ 2
  kernelDmat = exp(-dmat ^ 2 / eps)
  P = kernelDmat / colSums(kernelDmat) # Symmetric
  M = P
  
  
  Rcpp::sourceCpp("src/speedup.cpp", verbose = T)
  rst = list()
  for(i in 1:30)
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
  
  
  M1 = P; M2 = rst[[1]]; M4 = rst[[2]]; M8 = rst[[3]]; M16 = rst[[4]]; M32 = rst[[5]]; M64 = rst[[6]]; M128 = rst[[7]]; M256 = rst[[8]]
  M18 = matmul(M16, M2, maxCore = 15); M18 = M18 / rowSums(M18)
  # tmp = M2 / stationaryP ^ 0.5
  # tmp = M32 / stationaryP ^ 0.5
  tmp = M18 / stationaryP ^ 0.5
  if(T) # SVD approach
  {
    # tmp = apply(tmp, 2, function(x) (x - mean(x)) / sd(x))
    Msvd = svd(t(tmp))
  }
  if(F)
  {
    tmp = eigen(as.matrix(dist(tmp))) 
  }
  
  
  pdf("figure/circleFace", width = 8, height = 8 * 0.618)
  par(mar = c(4.1, 5, 0, 0), family = "serif")
  d = Msvd$d
  plot(d, xlab = "Eigen index", ylab = "Eigen value", cex.lab = 2, cex.axis = 1.5, bty = "L")
  dev.off()
  
  
  plot(Msvd$u[, 1])
  Xnew = Msvd$u[,1:2] %*% diag(Msvd$d[1:2])
  
  
  
  
  
  
  
  Xnew = Msvd$u[,1:10] %*% diag(Msvd$d[1:10])
  rgl::points3d(x = Xnew[,1], y = Xnew[,2], z = Xnew[,3])
  rgl::points3d(x = Xnew[,2], y = Xnew[,3], z = Xnew[,4])
  rgl::points3d(x = Xnew[,3], y = Xnew[,4], z = Xnew[,5])
  rgl::points3d(x = Xnew[,4], y = Xnew[,5], z = Xnew[,6])
  
  
  tmp = rev(cumsum(rev(Msvd$d) ^ 2))
  tmp = tmp / tmp[1]
  plot(tmp[1:50])
  
  
  
  
  Npc = sum(tmp >= 1 - 0.99)
  Xnew = tmp2$u[,1:Npc] %*% diag(tmp2$d[1:Npc])
  tXnew = t(Xnew)
  
  
  set.seed(42)
  k = 6L
  rstCluster = list()
  for(i in 1:30)
  {
    rstCluster[[i]] = GMKMcharlie::KM(X = tXnew, centroid = tXnew[, sample(1:nrow(Xnew), k)], verbose = F, maxCore = 7)
  }
  tmp = unlist(lapply(rstCluster, function(x) sum(unlist(lapply(x, function(y) sum(y$member2centroidDistance ^ 2))))))
  rstClusterOpt = lapply(rstCluster[[which.min(tmp)]], function(x) x$clusterMember)
  
  
  palette = randomcoloR::distinctColorPalette(k)
  # tmp = as.data.frame(lapply(1:ncol(bestClustering[[1]]), function(i) unlist(lapply(bestClustering, function(x) x[[i]]))))
  cols = unlist(mapply(function(x, y) rep(x, length(y)), palette, rstClusterOpt, SIMPLIFY = T))
  pdf("figure/K6circlefaceClusteringDiffusionmap.pdf", width = 8, height = 8 * 0.618)
  par(mar = c(4.1, 5, 0, 0), family = "serif")
  plot(X, col = cols, bty = "L", cex.lab = 2, cex.axis = 1.5, xlab = "x", ylab = "y")
  dev.off()
  
  
  
  
  
}




# Diffusion map, production.
if(F)
{
  
  
  load("data/circleface.Rdata")
  X = circleface
  tX = t(X)
  
  
  dmat = as.matrix(dist(X))
  eps = max(dmat) ^ 2
  kernelDmat = exp(-dmat ^ 2 / eps)
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
  M1 = P
  M1 = apply(M1, 2, function(x) x - mean(x))
  M1svd = svd(t(M1 / stationaryP ^ 0.5))
  
  
  pdf("figure/circleFaceDiffusionMapStationary.pdf", width = 8, height = 8 * 0.618)
  par(mar = c(4.1, 5, 0, 0), family = "serif")
  plot(stationaryP[1,], xlab = "Data point index", ylab = "Value in 1D", cex.lab = 2, cex.axis = 1.5, bty = "L")
  dev.off()
  Xnew1d = stationaryP[1,]
  
  
  # Largest 5 eigen values: 55.89, 8.52, 8.21, 0.81, 0.77
  pdf("figure/circleFaceDiffusionMapIn2D0step.pdf", width = 8, height = 8 * 0.618)
  par(mar = c(4.1, 5, 0, 0), family = "serif")
  plot(M1svd$u[, 1:2], xlab = "x", ylab = "y", cex.lab = 2, cex.axis = 1.5, bty = "L")
  dev.off()
  Xnew2d = M1svd$u[, 1:2] %*% diag(M1svd$d[1:2])
  
  
  pdf("figure/circleFaceDiffusionMapIn3D0step.pdf", width = 8, height = 8 * 0.618)
  par(mar = c(0, 0, 0, 0), family = "serif")
  plot3D::points3D(x = M1svd$u[, 1], y = M1svd$u[, 2], z = M1svd$u[, 3], phi = 30, theta = -40)
  dev.off()
  Xnew3d = M1svd$u[, 1:3] %*% diag(M1svd$d[1:3])
  
  
  
  
  source("R/rfuns.r")
  # Cluster diffused 2D data
  if(T)
  {
    
    
    load("data/circleface.Rdata")
    X = circleface
    tX = t(X)
    originalX = X; originaltX = tX
    
    
    X = Xnew2d; X = apply(X, 2, function(x) (x - mean(x)) / sd(x))
    tX = t(X)
    set.seed(123)
    k = 5L
    rst = list()
    for(i in 1:100)
    {
      rst[[length(rst) + 1]] = GMKMcharlie::KM(X = tX, centroid = tX[, sample(1:nrow(X), k)], verbose = F, maxCore = 7)
    }
    bestClustering = rst[[which.min(unlist(lapply(rst, function(u) sum(unlist(lapply(u, function(v) sum(v$member2centroidDistance ^ 2)))))))]]
    bestClustering = lapply(bestClustering, function(x) as.data.frame(X[x$clusterMember, ]))
    
    
    palette = randomcoloR::distinctColorPalette(length(bestClustering))
    tmp = as.data.frame(lapply(1:ncol(bestClustering[[1]]), function(i) unlist(lapply(bestClustering, function(x) x[[i]]))))
    cols = unlist(mapply(function(x, y) rep(x, nrow(y)), palette, bestClustering, SIMPLIFY = T))
    pdf("figure/diffusionmap2dcirclefaceClusteringInDiffusionSpace.pdf", width = 8, height = 8 * 0.618)
    par(mar = c(4.1, 5, 0, 0), family = "serif")
    plot(tmp, col = cols, bty = "L", cex.lab = 2, cex.axis = 1.5, xlab = "x", ylab = "y")
    dev.off()
    
    
    
    
    X = Xnew1d
    set.seed(123)
    k = 4L
    rst = list()
    for(i in 1:30)
    {
      rst[[i]] = kmeans(X, centers = X[sample(length(X), k)])
    }
    rst = rst[[which.min(unlist(lapply(rst, function(x) x$tot.withinss)))]]
    rst = data.frame(ind = 1:length(X), cid = rst$cluster)
    
    
    palette = randomcoloR::distinctColorPalette(length(unique(rst$cid)))
    cols = palette[rst[order(X), ]$cid]
    pdf("figure/diffusionmap1dcirclefaceClusteringInDiffusionSpace.pdf", width = 8, height = 8 * 0.618)
    par(mar = c(4.1, 5, 0, 0), family = "serif")
    plot(sort(X), col = cols, bty = "L", cex.lab = 2, cex.axis = 1.5, xlab = "Sorted data point index", ylab = "y")
    dev.off()
    
    
    pdf("figure/diffusionmap1dcirclefaceClusteringInOriginalSpace.pdf", width = 8, height = 8 * 0.618)
    par(mar = c(4.1, 5, 0, 0), family = "serif")
    cols = palette[rst$cid]
    plot(originalX, col = cols, bty = "L", cex.lab = 2, cex.axis = 1.5, xlab = "x", ylab = "y")
    dev.off()
   
     
  } 
}
  
  
  
  
  
  
  
  
  
  
  
  
  
  