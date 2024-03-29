/* Copyright(c)  Ryuichiro Nakato <rnakato@iqb.u-tokyo.ac.jp>
 * All rights reserved.
 */
#include <iostream>
#include <boost/format.hpp>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_sf.h>
#include "statistics.hpp"
#include "../submodules/SSP/common/inline.hpp"

double getStdNormDistP(const double x)
{
  return gsl_cdf_ugaussian_P(x);
}

double getStdNormDistQ(const double x)
{
  return gsl_cdf_ugaussian_Q(x);
}

double getNormDist(const double x, double m, double myu)
{
  double z((x - m)/myu);
  double d(getStdNormDistP(z));
  double p(std::min(d,1-d) *2);  // two-sided test
  return p;
}

double _getPoisson(const int32_t i, const double m)
{
  return gsl_ran_poisson_pdf(i, m);
}

double _getNegativeBinomial(const int32_t k, const double p, const double n)
{
  return gsl_ran_negative_binomial_pdf(k, p, n);
}


double getBinomial(const int32_t k, const double p, const int32_t n)
{
  return gsl_ran_binomial_pdf(k, p, n);
}

/*double _getZIP(const int32_t k, const double p, const double p0)
{
  double r(0);
  if (!k) r = p0 + (1 - p0) * gsl_ran_poisson_pdf(k, p);
  else    r =      (1 - p0) * gsl_ran_poisson_pdf(k, p);

  return r;
  }*/

double _getZINB(const int32_t k, const double p, const double n, const double p0)
{
  double r(0);
  if (!k) r = p0 + (1 - p0) * gsl_ran_negative_binomial_pdf(0, p, n);
  else    r =      (1 - p0) * gsl_ran_negative_binomial_pdf(k, p, n);

  return r;
}

gsl_multimin_fminimizer *gsl_multimin_fminimizer_new(const size_t ndim)
{
  const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex;  /* ネルダーとミードのシンプレックス法 */
  gsl_multimin_fminimizer            *s = gsl_multimin_fminimizer_alloc(T, ndim);
  return s;
}

gsl_vector *gsl_vector_new(const int32_t ndim, const double init)
{
  gsl_vector *p = gsl_vector_alloc(ndim);
  gsl_vector_set_all(p, init);
  return p;
}

double f_zinb_const(const gsl_vector *v, void *params)
{
  double *par = (double *)params;
  double p  = gsl_vector_get(v, 0);
  double n  = gsl_vector_get(v, 1);
  double p0 = gsl_vector_get(v, 2);
  if(p <= 0) p = 0.01;
  if(p >= 1) p = 0.99;
  if(p0 < 0) p0 = 0;
  if(p0 > 1) p0 = 1.0;

  double r(0), fxy(0);
  int32_t thre = par[0];
  for(int32_t i=0; i<thre; ++i) {
    if(!i) r = p0 + (1 - p0) * gsl_ran_negative_binomial_pdf(0, p, n);
    else   r =      (1 - p0) * gsl_ran_negative_binomial_pdf(i, p, n);
    fxy += (par[i+1] - r)*(par[i+1] - r);
    //std::cout << boost::format("i=%1% par[i+1]=%2% r=%3% xx=%4%\n") % i % par[i+1] % r % ((par[i+1] - r)*(par[i+1] - r));
  }
  //  std::cout << boost::format("fxy=%1% p=%2% n=%3% p0=%4% thre=%5%\n") % fxy % p % n % p0 % thre;
  return fxy;
}

#ifdef DEBUG
void func_iteration(gsl_multimin_fminimizer *s, size_t ndim)
#else
void func_iteration(gsl_multimin_fminimizer *s)
#endif
{
  size_t iter(0);
  int32_t status;
  double size;
  do {
    ++iter;
    status = gsl_multimin_fminimizer_iterate(s);  /* 繰り返し計算を1回行う。予期しない問題が発生した場合はエラーコードを返す。 */
    if(status) break;
    size = gsl_multimin_fminimizer_size(s);       /* sのその時点での最小点の最良推定値を返す */
    status = gsl_multimin_test_size(size, 1e-3);  /* sizeが閾値(1e-3)より小さければGSL_SUCCESS を、そうでなければGSL_CONTINUEを返す。 */

#ifdef DEBUG
    if(status == GSL_SUCCESS) std::cout << "converged to minimum at " << iter << std::endl;
    if(ndim==2) std::cout << boost::format("%1% p=%2% p0=%3% f() = %4% size = %5%\n") % iter % gsl_vector_get(s->x, 0) % gsl_vector_get(s->x, 1) % gsl_multimin_fminimizer_minimum(s) % size;
    else std::cout << boost::format("%1% p=%2% n=%3% p0=%4% f() = %5% size = %6%\n") % iter % gsl_vector_get(s->x, 0) % gsl_vector_get(s->x, 1) % gsl_vector_get(s->x, 2) % gsl_multimin_fminimizer_minimum(s) % size;
#endif
  } while (status == GSL_CONTINUE && iter < 1000);

  return;
}

void iterateZINB(void *par, double nb_p_pre, double nb_n_pre, double &nb_p, double &nb_n, double &nb_p0)
{
  size_t ndim(3);
  gsl_multimin_fminimizer *s = gsl_multimin_fminimizer_new(ndim);
  gsl_vector *x = gsl_vector_alloc(ndim);
  gsl_vector_set(x, 0, nb_p_pre);
  gsl_vector_set(x, 1, nb_n_pre);
  gsl_vector_set(x, 2, 0.2); // p0
  gsl_vector *ss = gsl_vector_new(ndim, 0.1); // step size

  gsl_multimin_function minex_func;
  minex_func.n = ndim;
  minex_func.f = &f_zinb_const;
  minex_func.params = par;
  gsl_multimin_fminimizer_set(s, &minex_func, x, ss);

#ifdef DEBUG
  func_iteration(s, ndim);
#else
  func_iteration(s);
#endif

  nb_p  = gsl_vector_get(s->x, 0);
  if(nb_p <= 0) nb_p = 0.01;
  if(nb_p >= 1) nb_p = 0.99;
  nb_n  = gsl_vector_get(s->x, 1);
  nb_p0 = gsl_vector_get(s->x, 2);
  if(nb_p0 < 0) nb_p0 = 0.0;
  if(nb_p0 > 1) nb_p0 = 1.0;

  gsl_vector_free(x);
  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free(s);
  return;
}

double f_poisson(const gsl_vector *v, void *params)
{
  double *par = (double *)params;
  double p = gsl_vector_get(v, 0);
  double p0 = gsl_vector_get(v, 1);
  if(p <= 0) p = 0.01;
  if(p >= 1) p = 0.99;

  double fxy(0);
  int32_t thre = par[0];
  double r;
  for(int32_t i=0; i<thre; ++i) {
    if(!i) r = p0 + (1 - p0) * _getPoisson(i, p);
    else   r =      (1 - p0) * _getPoisson(i, p);
    fxy += (par[i+1] - r)*(par[i+1] - r);
  }
  return fxy;
}

void iteratePoisson(void *par, double ave_pre, double &ave, double &p0)
{
  size_t ndim(2);

  gsl_multimin_fminimizer *s = gsl_multimin_fminimizer_new(ndim);
  gsl_vector *x = gsl_vector_alloc(ndim);
  gsl_vector_set(x, 0, ave_pre);
  gsl_vector_set(x, 1, 0.02); // p0
  gsl_vector *ss = gsl_vector_new(ndim, 0.1); // step size

  gsl_multimin_function minex_func;
  minex_func.n = ndim;
  minex_func.f = &f_poisson;
  minex_func.params = par;
  gsl_multimin_fminimizer_set(s, &minex_func, x, ss);

#ifdef DEBUG
  func_iteration(s, ndim);
#else
  func_iteration(s);
#endif

  ave = gsl_vector_get(s->x, 0);
  if(ave <= 0) ave = 0.01;
  if(ave >= 1) ave = 0.99;
  p0 = gsl_vector_get(s->x, 1);
  if(p0 < 0) p0 = 0.0;
  if(p0 > 1) p0 = 1.0;

  gsl_vector_free(x);
  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free(s);
  return;
}

double getlogpZINB(double k, double p, double n)
{
  double r, pval;
  if(!k) pval = 0;
  else pval = gsl_cdf_negative_binomial_Q(k, p, n);
  if(!pval) r = 0; else r = -log10(pval);
  return r;
}
