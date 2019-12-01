#include "Params.h"

struct HPR_ZZq {
  std::array<std::array<value_t, h1>, n> a_1;
  std::array<std::array<value_t, h2>, n> a_2;
  std::array<value_t, n> a_sk;
};

struct HPR {
  HPR ();
  HPR_ZZq to_ZZq (NTL::ZZ a);
  NTL::ZZ to_ZZ (const HPR_ZZq &A);
  void mul (HPR_ZZq &c, const HPR_ZZq &a, const HPR_ZZq &b);
  HPR_ZZq mul (const HPR_ZZq &a, const HPR_ZZq &b);
};
