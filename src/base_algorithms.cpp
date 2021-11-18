#include "nnlm.h"


int update_with_marker(mat & H, const mat & Wt, const mat & A, const umat & mask, const bool useMask,
                       unsigned int max_iter, double rel_tol, int n_threads)
{

  unsigned int m = A.n_cols;
  int total_raw_iter = 0;
  if (n_threads < 0) n_threads = 0;
  bool is_masked = !mask.empty();

  mat WtW;
  vec mu, sumW;
  sumW = sum(Wt, 1);

#pragma omp parallel for num_threads(n_threads) schedule(dynamic) private(mu)
  for (unsigned int j = 0; j < m; j++) // by columns of H
  {
    // break if all entries of col_j are masked
    if (is_masked && arma::all(mask.col(j)) && useMask)
      continue;


    int iter = 0;

    iter = scd_kl_marker_update(H.col(j), Wt, A.col(j), sumW, mask.col(j),useMask, max_iter, rel_tol);

#pragma omp critical
    total_raw_iter += iter;
  }
  return total_raw_iter;
}




int scd_kl_marker_update(subview_col<double> Hj, const mat & Wt, const vec & Aj, const vec & sumW, const subview_col<uword> mask,
                         const bool useMask,  const unsigned int & max_iter, const double & rel_tol)
{
  // Problem:  Aj = W * Hj
  // Method: Sequentially minimize KL distance using quadratic approximation
  // Wt = W^T
  // sumW = column sum of W
  // mask: skip updating
  // beta: a vector of 3, for L2, angle, L1 regularization

  // Rprintf("elem %10d ; col  %10d ; row  %10d \n", Hj.n_elem,Hj.n_cols,Hj.n_rows );



  double sumHj = sum(Hj);
  vec Ajt = Wt.t()*Hj;
  vec mu;
  double a; // 2nd-order-derivative
  double b; // 1st-order-derivative
  double tmp, etmp;
  double rel_err = 1 + rel_tol;
  bool is_masked = mask.n_elem > 0;
  bool doUpdate = true;

  // if(useMask){
  //   Rprintf(" Updating W genes, sum mask %d ; is_msked %d\n", sum(mask),is_masked);
  // }
  //   Rprintf(" Updating H samples %d\n", sum(mask));
  // }

  //Rprintf("MASK : elem %3d ; col  %3d ; row  %3d ; is_masked %d ;sumAsked %d \n",
  //      mask.n_elem,mask.n_cols,mask.n_rows,is_masked, sum(mask) );

  unsigned int t = 0;
  for (; t < max_iter && rel_err > rel_tol; t++)
  {
    rel_err = 0;
    for (unsigned int k = 0; k < Wt.n_rows; k++)
    {

      mu = Wt.row(k).t()/(Ajt + TINY_NUM);
      a = dot(Aj, square(mu));
      b = dot(Aj, mu) - sumW(k); // 0.5*ax^2 - bx
      //a += beta(0);
      //b += a*Hj(k) - beta(2) - beta(1)*(sumHj - Hj(k));
      b += a*Hj(k) ;
      tmp = b/(a+TINY_NUM);
      if (tmp < 0) tmp = 0;
      doUpdate= tmp != Hj(k);



      if(useMask)
      {
        // Rprintf("k%d, wtnrows %d ",k,Wt.n_rows);
        // if(k==(Wt.n_rows)-1)Rprintf("needUpdate %d",doUpdate);
        //Rprintf("UsingMask, doUpdate %d, Updating W genes %d, is masekd %d\n",doUpdate, sum(mask),mask(k));
        //Rprintf("UsingMask, Updating W genes %d\n", sum(mask));
        if (is_masked && mask(k) > 0 && tmp < max(Hj) )
        {
          doUpdate=false;
        } else if (sum(mask)>0 && mask(k) == 0 && tmp > max(Hj) ){
          doUpdate=false;
        }
        // if(k==(Wt.n_rows)-1)Rprintf(" doingUpdate %d, has masked %d, is masekd %d, k %d\n", doUpdate, sum(mask),mask(k),k);
      }

      if (doUpdate)
      {
        Ajt += (tmp - Hj(k)) * Wt.row(k).t();
        etmp = 2*std::abs(Hj(k)-tmp) / (tmp+Hj(k) + TINY_NUM);
        if (etmp > rel_err)
          rel_err = etmp;
        sumHj += tmp - Hj(k);
        Hj(k) = tmp;
      }
    }
  }
  return int(t);
}
