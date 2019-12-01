#include "HPR.h"
#include <NTL/mat_ZZ.h>
#include <NTL/vec_ZZ.h>
#include <NTL/LLL.h>
#include <assert.h>

/* precomputations */
std::array<value_t, h1> b1i_over_B1_mod_b1i;
std::array<std::array<value_t, h1>, h2> B1_over_b1j_mod_b2i;
std::array<value_t, h2> B1_inv_mod_b2i;
std::array<value_t, h2> b2i_over_B2;
std::array<value_t, h1> B1_over_b1i_sk;
value_t B1_inv_sk;
std::array<value_t, h2> b2i_inv_sk;
value_t B2_inv_sk;
std::array<std::array<value_t, h2>, h1> B2_over_b2i_mod_b1i;
std::array<value_t, h1> B2_mod_b1i;
std::array<value_t, h2> B1_mod_b2i;

/* channel arithmetic */
std::array<value_t, h1> c1;
std::array<value_t, h2> c2;
value_t mask;
value_t mask_sk;
  
value_t reduce (greater_value_t x, value_t c, value_t p);
value_t addmod (greater_value_t a, value_t b, value_t p);
value_t submod (value_t a, value_t b, value_t p);
value_t mulmod (greater_value_t a, value_t b, value_t c, value_t p);
value_t mulbeta (greater_value_t a, value_t c, value_t p);
value_t mulbeta_sk (value_t a);

template<size_t h>
std::array<std::array<value_t, h>, n>
rns_pol (NTL::ZZ x,
	 const std::array<value_t, h> &ps)
{
  std::array<std::array<value_t, h>, n> x_red;
  if (x > P/2) x -= P;

  for (size_t i = 0; i < n-1; i++) {
    NTL::ZZ xi = x % B1;
    if (xi > B1/2) xi -= B1;
    x = (x - xi) / B1;
    for (size_t j = 0; j < h; j++) {
      NTL::ZZ x_red_ij = xi % NTL::ZZ (ps[j]);
      x_red[i][j] = NTL::conv<value_t> (x_red_ij);
    }
  }

  //last digit may be larger than B1/2
  for (size_t j = 0; j < h; j++) {
    NTL::ZZ x_red_ij = x % NTL::ZZ (ps[j]);
    x_red[n-1][j] = NTL::conv<value_t> (x_red_ij);
  }

  return x_red;
}

template<size_t h1, size_t h2>
NTL::ZZ
rns_pol_inv (const std::array<std::array<value_t, h1>, n> &x1,
	     const std::array<std::array<value_t, h2>, n> &x2,
	     const std::array<value_t, h1> &ps1,
	     const std::array<value_t, h2> &ps2,
	     size_t k2)
{
  NTL::ZZ res (0);
  NTL::ZZ gamma_i (1);

  for (size_t i = 0; i < n; i++) {
    NTL::ZZ X;
    NTL::ZZ P;

    if (x1[i][0] > ps1[0]/2) {
      X = NTL::ZZ (x1[i][0]) - NTL::ZZ (ps1[0]);
    } else {
      X = NTL::ZZ (x1[i][0]);
    }

    P = ps1[0];

    for (size_t j = 1; j < h1; j++) {
      NTL::CRT (X, P, NTL::ZZ (x1[i][j]), NTL::ZZ (ps1[j]));
    }
    for (size_t j = 0; j < k2; j++) {
      NTL::CRT (X, P, NTL::ZZ (x2[i][j]), NTL::ZZ (ps2[j]));
    }

    res = res + X * gamma_i;
    gamma_i = gamma_i * B1;
  }

   return res;
}

HPR::HPR ()
{
  /* channel arithmetic */
  for (size_t i = 0; i < h1; i++) {
    c1[i] = (((greater_value_t)1)<<logw) - b1[i];
  }
  for (size_t i = 0; i < h2; i++) {
    c2[i] = (((greater_value_t)1)<<logw) - b2[i];
  }

  mask = (((greater_value_t)1)<<logw) - 1;
  mask_sk = (((greater_value_t)1)<<log_bsk) - 1;

  /* precomputations */
  NTL::ZZ psk (1);
  psk <<= log_bsk;

  NTL::ZZ B2 (1);
  for (size_t i = 0; i < h2; i++) {
    B2 *= NTL::ZZ (b2[i]);
  }
  NTL::ZZ B1 (1);
  for (size_t i = 0; i < h1; i++) {
    B1 *= NTL::ZZ (b1[i]);
  }
  
  for (size_t i = 0; i < h2; i++) {
    NTL::ZZ B2_over_b2i (B2 / NTL::ZZ (b2[i]));
    b2i_over_B2[i] = NTL::conv<value_t>
      (NTL::InvMod (B2_over_b2i % NTL::ZZ (b2[i]), NTL::ZZ (b2[i])));

    for (size_t j = 0; j < h1; j++) {
      B2_over_b2i_mod_b1i[j][i] = NTL::conv<value_t>
	(B2_over_b2i % NTL::ZZ (b1[j]));
    }

    B1_inv_mod_b2i[i] = NTL::conv<value_t>
      (NTL::InvMod (B1 % NTL::ZZ (b2[i]), NTL::ZZ (b2[i])));

    B1_mod_b2i[i] = NTL::conv<value_t>
      (B1 % NTL::ZZ (b2[i]));
  }
  
  for (size_t i = 0; i < h1; i++) {
    NTL::ZZ B1_over_b1i (B1 / NTL::ZZ (b1[i]));
    b1i_over_B1_mod_b1i[i] = NTL::conv<value_t>
      (NTL::InvMod (B1_over_b1i % NTL::ZZ (b1[i]), NTL::ZZ (b1[i])));
    
    for (size_t j = 0; j < h2; j++) {
      B1_over_b1j_mod_b2i[j][i] = NTL::conv<value_t>
	(B1_over_b1i % NTL::ZZ (b2[j]));
    }

    B1_over_b1i_sk[i] = NTL::conv<value_t>
      (B1_over_b1i % psk);

    B2_mod_b1i[i] = NTL::conv<value_t>
      (B2 % NTL::ZZ (b1[i]));
  }

  for (size_t i = 0; i < h2; i++) {
    b2i_inv_sk[i] = NTL::conv<value_t>
      (NTL::InvMod (NTL::ZZ (b2[i]) % psk, psk));
  }
  B1_inv_sk = NTL::conv<value_t>
    (NTL::InvMod (B1 % psk, psk));
  B2_inv_sk = NTL::conv<value_t>
    (NTL::InvMod (B2 % psk, psk));
}

value_t
reduce (greater_value_t x, value_t c, value_t p)
{
  const greater_value_t bound = ((greater_value_t)1)<<logw;
  
  while (x >= bound) {
    greater_value_t x0 = x & mask;
    greater_value_t x1 = x >> logw;
    x = x0 + c * x1;
  }

  if (x >= p) x -= p;

  return x;
}

value_t
addmod (greater_value_t a, value_t b, value_t p)
{
  greater_value_t c = a + b;
  if (c >= p) c -= p;
  return c;
}

value_t
submod (value_t a, value_t b, value_t p)
{
  return (b == 0 ? a : addmod (a, p-b, p));
}

value_t
mulmod (greater_value_t a, value_t b, value_t c, value_t p)
{
  greater_value_t d = a*b;
  return reduce (d, c, p);
}

template<typename T, size_t beta1, size_t bit, size_t i>
struct _mul_beta { };

template<typename T, size_t i>
struct _mul_beta<T, 0, 0, i>
{
  T apply (const T &x)
  {
    return 0;
  }
};

template<typename T, size_t beta1, size_t i>
struct _mul_beta<T, beta1, 0, i>
{
  T apply (const T &x)
  {
    return _mul_beta<T, (beta1>>1), (beta1&1), i+1>{}.apply (x);
  }
};

template<typename T, size_t beta1, size_t i>
struct _mul_beta<T, beta1, 1, i>
{
  T apply (const T &x)
  {
    return (x<<i) +
      _mul_beta<T, (beta1>>1), (beta1&1), i+1>{}.apply (x);
  }
};

value_t
mulbeta (greater_value_t a, value_t c, value_t p)
{
  greater_value_t d = _mul_beta<greater_value_t, (beta>>1), beta&1, 0>{}.apply (a);
  
  return reduce (d, c, p);
}

value_t
mulbeta_sk (value_t a)
{
  value_t d = _mul_beta<value_t, (beta>>1), beta&1, 0>{}.apply (a);
  
  return d & mask_sk;
}

HPR_ZZq
HPR::to_ZZq (NTL::ZZ x)
{
  auto a_1 = rns_pol (x, b1);
  auto a_2 = rns_pol (x, b2);

  std::array<value_t, n> a_sk;
  NTL::ZZ psk (1);
  psk <<= log_bsk;
  
  if (x > P/2) x -= P;
  for (size_t i = 0; i < n-1; i++) {
    NTL::ZZ xi = x % B1;
    if (xi > B1/2) xi -= B1;
    x = (x - xi) / B1;
    a_sk[i] = NTL::conv<value_t>
      (xi % psk);
  }
  a_sk[n-1] = NTL::conv<value_t> (x % psk);

  return HPR_ZZq
    {
      a_1, a_2, a_sk
    };
}

NTL::ZZ
HPR::to_ZZ (const HPR_ZZq &A)
{
  auto x = rns_pol_inv (A.a_1, A.a_2, b1, b2, 1);

  return x % P;
}

template<size_t h>
void
rns_pol_mul_beta (std::array<std::array<value_t, h>, n> &zs,
		  const std::array<std::array<value_t, h>, n> &xs,
		  const std::array<std::array<value_t, h>, n> &ys,
		  const std::array<value_t, h> &cs,
		  const std::array<value_t, h> &ps)
{
  for (size_t i = 0; i < h; i++) {
    for (size_t j = 0; j < n; j++) {
      zs[j][i] = 0;

      for (size_t k = j+1; k < n; k++) {
	value_t prod = mulmod (xs[k][i], ys[n+j-k][i], cs[i], ps[i]);
	zs[j][i] = addmod (zs[j][i], prod, ps[i]);
      }
      zs[j][i] = mulbeta (zs[j][i], cs[i], ps[i]);

      for (size_t k = 0; k <= j; k++) {
	value_t prod = mulmod (xs[k][i], ys[j-k][i], cs[i], ps[i]);
	zs[j][i] = addmod (zs[j][i], prod, ps[i]);
      }
    }
  }
}

void
rns_pol_mul_sk_beta (std::array<value_t, n> &zs,
		     const std::array<value_t, n> &xs,
		     const std::array<value_t, n> &ys)
{
  for (size_t j = 0; j < n; j++) {
    zs[j] = 0;

    for (size_t k = j+1; k < n; k++) {
      value_t prod = xs[k] * ys[n+j-k];
      zs[j] += prod;
    }
    zs[j] = mulbeta_sk (zs[j]);

    for (size_t k = 0; k <= j; k++) {
      value_t prod = xs[k] * ys[j-k];
      zs[j] += prod;
    }

    zs[j] &= mask_sk;
  }
}

template<size_t h1, size_t h2>
void
rns_basis_extension (std::array<value_t, h2> &zs,
		     const std::array<value_t, h1> &xs,
		     const std::array<std::array<value_t, h1>, h2> &bs,
		     const std::array<value_t, h2> cs2,
		     const std::array<value_t, h2> ps2,
		     const std::array<value_t, h1> ps1)
{
  for (size_t i = 0; i < h2; i++) {
    zs[i] = 0;
      
    for (size_t k = 0; k < h1; k++) {
      value_t pdiff = (ps1[k] > ps2[i] ? ps2[i]<<1 : ps2[i]) - ps1[k];
      value_t xskj = xs[k];
      bool negative = false;
      if (xskj > ps1[k]/2) {
	negative = true;
      }
      if (xskj > ps2[i]) {
	xskj -= ps2[i];
      }
      if (negative) {
	xskj = addmod (xskj, pdiff, ps2[i]);
      }
      value_t prod = mulmod (xskj, bs[i][k], cs2[i], ps2[i]);
      zs[i] = addmod (zs[i], prod, ps2[i]);
    }
  }
}
			 
template<size_t h1>
void
rns_basis_extension_sk (value_t &zs,
			const std::array<value_t, h1> &xs,
			const std::array<value_t, h1> &bs,
			const std::array<value_t, h1> &ps1)
{
  zs = 0;
      
  for (size_t k = 0; k < h1; k++) {
    value_t xskj = xs[k];
    if (xskj > ps1[k]/2) {
      xskj -= ps1[k];
    }
    value_t prod = xskj * bs[k];
    zs += prod;
  }

  zs &= mask_sk;
}

template<size_t h>
void
rns_scale (std::array<value_t, h> &zs,
	   const std::array<value_t, h> &xs,
	   const std::array<value_t, h> &ys,
	   const std::array<value_t, h> &cs,
	   const std::array<value_t, h> &ps)
{
  for (size_t i = 0; i < h; i++) {
    zs[i] = mulmod (xs[i], ys[i], cs[i], ps[i]);
  }
}

template<size_t h>
void
rns_scale_beta (std::array<value_t, h> &zs,
		const std::array<value_t, h> &xs,
		const std::array<value_t, h> &cs,
		const std::array<value_t, h> &ps)
{
  for (size_t i = 0; i < h; i++) {
    zs[i] = mulbeta (xs[i], cs[i], ps[i]);
  }
}

void
rns_scale_sk (value_t &zs,
	      const value_t &xs,
	      value_t y)
{
  zs = xs * y;
  zs &= mask_sk;
}

void
rns_scale_beta_sk (value_t &zs,
		   const value_t &xs)
{
  zs = mulbeta_sk (xs);
}

template<size_t h>
void
rns_add (std::array<value_t, h> &zs,
	 const std::array<value_t, h> &xs,
	 const std::array<value_t, h> &ys,
	 const std::array<value_t, h> &ps)
{
  for (size_t i = 0; i < h; i++) {
    zs[i] = addmod (xs[i], ys[i], ps[i]);
  }
}

template<size_t h>
void
rns_sub (std::array<value_t, h> &zs,
	 const std::array<value_t, h> &xs,
	 const std::array<value_t, h> &ys,
	 const std::array<value_t, h> &ps)
{
  for (size_t i = 0; i < h; i++) {
    zs[i] = submod (xs[i], ys[i], ps[i]);
  }
}

void
rns_add_sk (value_t &zs,
	    const value_t &xs,
	    const value_t &ys)
{
  zs = xs + ys;
  zs &= mask_sk;
}

void
rns_sub_sk (value_t &zs,
	    const value_t &xs,
	    const value_t &ys)
{
  zs = xs - ys;
  zs &= mask_sk;
}

void
rns_correct (std::array<value_t, h1> &a_1,
	     const value_t &alpha,
	     const std::array<value_t, h1> B2_mod_b1i)
{
  for (size_t i = 0; i < h1; i++) {
    value_t alpha_j = alpha;
    bool negative = false;
    if (alpha_j > (((value_t)1)<<(log_bsk-1))) {
      negative = true;
      alpha_j = -alpha_j;
      alpha_j &= mask_sk;
    }
    if (alpha_j > b1[i])
      alpha_j -= b1[i];
    value_t prod = mulmod (alpha_j, B2_mod_b1i[i], c1[i], b1[i]);
    if (negative) {
      prod = b1[i] - prod;
    }
    a_1[i] = submod (a_1[i], prod, b1[i]);
  }
}

template<size_t h>
void
rns_import_from_sk (std::array<value_t, h> &a_1,
		    value_t a_sk,
		    const std::array<value_t, h> &ps)
{
  bool negative = false;
  if (a_sk > (((value_t)1)<<(log_bsk-1))) {
    negative = true;
    a_sk = -a_sk;
    a_sk &= mask_sk;
  }

  for (size_t i = 0; i < h; i++) {
    value_t a_sk_i = a_sk;
    if (a_sk_i > ps[i])
      a_sk_i -= ps[i];
    if (negative) {
      a_1[i] = ps[i] - a_sk_i;
    } else {
      a_1[i] = a_sk_i;
    }
  }
}

void
HPR::mul (HPR_ZZq &c, const HPR_ZZq &a, const HPR_ZZq &b)
{
  rns_pol_mul_beta (c.a_1, a.a_1, b.a_1, c1, b1);
  rns_pol_mul_beta (c.a_2, a.a_2, b.a_2, c2, b2);
  rns_pol_mul_sk_beta (c.a_sk, a.a_sk, b.a_sk);

  for (size_t i = 0; i < n; i++) {
    std::array<value_t, h1> xi_1;
    rns_scale (xi_1, c.a_1[i], b1i_over_B1_mod_b1i, c1, b1);
    std::array<value_t, h2> q_2;
    rns_basis_extension (q_2, xi_1, B1_over_b1j_mod_b2i, c2, b2, b1);
    std::array<value_t, h2> carry_2;
    rns_sub (carry_2, c.a_2[i], q_2, b2);
    rns_scale (carry_2, carry_2, B1_inv_mod_b2i, c2, b2);
    std::array<value_t, h2> xi_2;
    rns_scale (xi_2, carry_2, b2i_over_B2, c2, b2);
    value_t q_sk;
    rns_basis_extension_sk (q_sk, xi_1, B1_over_b1i_sk, b1);
    value_t carry_sk;
    rns_sub_sk (carry_sk, c.a_sk[i], q_sk);
    rns_scale_sk (carry_sk, carry_sk, B1_inv_sk);
    value_t alpha_sk;
    rns_basis_extension_sk (alpha_sk, xi_2, b2i_inv_sk, b2);
    value_t tmp;
    rns_scale_sk (tmp, carry_sk, B2_inv_sk);
    rns_sub_sk (alpha_sk, alpha_sk, tmp);
    std::array<value_t, h1> carry_1;
    rns_basis_extension (carry_1, xi_2, B2_over_b2i_mod_b1i, c1, b1, b2);
    rns_correct (carry_1, alpha_sk, B2_mod_b1i);

    if (i < n-1) {
      rns_add (c.a_1[i+1], c.a_1[i+1], carry_1, b1);
      rns_add (c.a_2[i+1], c.a_2[i+1], carry_2, b2);
      rns_add_sk (c.a_sk[i+1], c.a_sk[i+1], carry_sk);
    } else {
        rns_scale_beta (carry_1, carry_1, c1, b1);
	rns_scale_beta (carry_2, carry_2, c2, b2);
	rns_scale_beta_sk (carry_sk, carry_sk);

	rns_add (c.a_1[0], c.a_1[0], carry_1, b1);
	rns_add (c.a_2[0], c.a_2[0], carry_2, b2);
	rns_add_sk (c.a_sk[0], c.a_sk[0], carry_sk);
    }

    c.a_2[i] = q_2;
    c.a_sk[i] = q_sk;
  }

  for (size_t i = 0; i < n-2; i++) {
    std::array<value_t, h1> xi_1;
    rns_scale (xi_1, c.a_1[i], b1i_over_B1_mod_b1i, c1, b1);
    value_t q_sk;
    rns_basis_extension_sk (q_sk, xi_1, B1_over_b1i_sk, b1);
    value_t carry_sk;
    rns_sub_sk (carry_sk, c.a_sk[i], q_sk);
    rns_scale_sk (carry_sk, carry_sk, B1_inv_sk);
    c.a_sk[i] = q_sk;
    std::array<value_t, h2> carry_2;
    rns_import_from_sk (carry_2, carry_sk, b2);
    std::array<value_t, h2> tmp;
    rns_scale (tmp, carry_2, B1_mod_b2i, c2, b2);
    rns_sub (c.a_2[i], c.a_2[i], tmp, b2);
    std::array<value_t, h1> carry_1;
    rns_import_from_sk (carry_1, carry_sk, b1);
    
    rns_add (c.a_1[i+1], c.a_1[i+1], carry_1, b1);
    rns_add (c.a_2[i+1], c.a_2[i+1], carry_2, b2);
    rns_add_sk (c.a_sk[i+1], c.a_sk[i+1], carry_sk);
  }
}

HPR_ZZq
HPR::mul (const HPR_ZZq &a, const HPR_ZZq &b)
{
  HPR_ZZq c;
  mul (c, a, b);
  return c;
}
