//disable armadillo matrix/vector boundary checking
#define ARMA_NO_DEBUG

#ifdef _OPENMP
#include <omp.h>
#endif

//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

//[[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>

#include <Rcpp.h>
#include <R.h>

#define TINY_NUM 1e-16
#define NNLM_REL_TOL 1e-8
#define NNMF_REL_TOL 1e-6
#define MAX_ITER 500
#define NNMF_INNER_MAX_ITER 10
#define N_THREADS 1
#define TRACE_STEP 10
#define SHOW_WARNING true
#define DEFAULT_METHOD 1

//using namespace Rcpp;
using namespace arma;



Rcpp::List c_semisupnnmf(const arma::mat & A, const unsigned int k, arma::mat W, arma::mat H, arma::umat Wm, arma::umat Hm,
              const unsigned int max_iter, const double rel_tol,
              const int n_threads, const int verbose, const bool show_warning, const unsigned int inner_max_iter,
              const double inner_rel_tol,  unsigned int trace);


int update_with_marker(mat & H, const mat & Wt, const mat & A, const umat & mask, const bool useMask,
                       unsigned int max_iter, double rel_tol, int n_threads);
int  scd_kl_marker_update(subview_col<double> Hj, const mat & Wt, const vec & Aj, const vec & sumW,
                           const subview_col<uword> mask,const bool useMask,
                           const unsigned int & max_iter, const double & rel_tol);
