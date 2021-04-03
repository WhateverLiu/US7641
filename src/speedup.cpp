// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::plugins(cpp11)]]
// # include <Rcpp.h>
# include <RcppArmadillo.h>
# include <RcppParallel.h>
# include <atomic>
using namespace Rcpp;
// # define vec std::vector


struct dynamicTasking
{
  std::size_t NofCore;
  std::size_t NofAtom;
  std::atomic<std::size_t> counter;
  
  
  void reset(std::size_t NofCPU, std::size_t NofTask)
  {
    NofCore = NofCPU;
    if(NofCore > NofTask) NofCore = NofTask;
    NofAtom = NofTask;
    counter = 0;
  }
  
  
  dynamicTasking(std::size_t NofCPU, std::size_t NofTask)
  {
    reset(NofCPU, NofTask);
  }
  
  
  bool nextTaskID(std::size_t &taskID, std::size_t increment = 1)
  {
    taskID = counter.fetch_add(increment);
    return taskID < NofAtom;
  }
};


template<typename ing, typename num>
inline void vecxscalerAdd(num *x, ing size, num scalar, num *rst)
{
  for(ing i = 0; i < size; ++i)
    rst[i] += x[i] * scalar;
}
// y is a column vector and of size ncol.
template<typename ing, typename num>
inline void matxcol(num *X, ing nrow, ing ncol, num *y, num *rst)
{
  std::fill(rst, rst + nrow, 0);
  for(std::size_t i = 0, iend = ncol; i < iend; ++i)
    vecxscalerAdd(X + i * nrow, nrow, y[i], rst);
}


template<typename ing, typename num>
struct paraMatMul: public RcppParallel::Worker
{
  ing grainSize;
  ing N, P; // N x P, P x M
  num *X, *Y, *rst;
  dynamicTasking *dT;
  void operator() (std::size_t st, std::size_t end)
  {
    for(;;)
    {
      std::size_t objI = 0;
      if(!dT->nextTaskID(objI, grainSize)) break; // Number of tasks equals M
      for(std::size_t i = objI, iend = std::min(i + grainSize, dT->NofAtom);
          i < iend; ++i)
      {
        matxcol(X, N, P, Y + P * i, rst + N * i);
      }
    }
  }
  paraMatMul(num *X, num *Y, num *rst,
             ing N, ing P, ing M, ing maxCore): // (N x P) x (P x M)
    N(N), P(P), X(X), Y(Y), rst(rst)
  {
    dynamicTasking dt(maxCore, M); dT = &dt;
    grainSize = M / (maxCore * maxCore);
    parallelFor(0, maxCore, *this);
  }
};


// [[Rcpp::export]]
NumericMatrix matmul(NumericMatrix X, NumericMatrix Y, int maxCore = 16) 
{
  NumericMatrix rst(X.nrow(), Y.ncol());
  paraMatMul<int, double> (
      &X[0], &Y[0], &rst[0], X.nrow(), X.ncol(), Y.ncol(), maxCore);
  return rst;
}


// [[Rcpp::export]]
arma::mat matmulArmadillo(arma::mat &X, arma::mat &Y)
{
  return X * Y;
}


// [[Rcpp::export]]
List svdArmadillo(arma::mat &X)
{
  arma::mat U, V;
  arma::vec s;
  arma::svd(U, s, V, X);
  return List::create(Named("d") = s, Named("u") = U, Named("v") = V);
}




















