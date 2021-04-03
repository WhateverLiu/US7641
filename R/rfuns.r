

# cl is a list of data points.
kmeansBIC = function(cl)
{
  X = as.matrix(as.data.frame(lapply(1:ncol(cl[[1]]), function(i) unlist(lapply(cl, function(x) x[[i]])))))
  tX = t(X)
  P = t(as.data.frame(lapply(cl, function(x)
  {
    mu = colMeans(x); sig = cov(x)
    tX_mu = tX - mu
    k = length(mu)
    tmp = diag(sig)
    diag(sig) = tmp + 1e-6
    rst = mvtnorm::dmvnorm(x = X, mean = mu, sigma = sig, log = T)
    rst
  })))
  colMax = apply(P, 2, function(x) max(x))
  colS = colSums(apply(P, 2, function(x) exp(x - max(x))))
  loglike = sum(log(colS) + colMax) / length(cl)
  length(cl) * ncol(cl[[1]]) * log(nrow(X)) - 2 * loglike
}


gmBIC = function(Ndata, d, Ncluster, P)
{
  (d + d * (d + 1) / 2) * Ncluster * log(Ndata) - 2 * sum(log(P))
}


kurto = function(X)
{
  if(is.vector(X)) return(mean( ( (X - mean(X)) / sd(X) ) ^ 4 ) - 3)
  apply(X, 2, function(x) mean( ( (x - mean(x)) / sd(x) ) ^ 4 ) - 3 )
}


epsilonCompute <- function(D, p = 0.01)
{
  n = dim(D)[1]
  k = ceiling(p * n)
  k = ifelse(k < 2, 2, k) # use k of at least 2
  D.sort = apply(D, 1, sort)
  dist.knn = D.sort[k + 1, ] # find dists. to kth nearest neighbor
  epsilon = 2 * median(dist.knn) ^ 2
  return(epsilon)
}


# gauKernel <- function(x, sig = )
























