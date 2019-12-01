#include "HPR.h"
#include <vector>
#include <chrono>
#include <assert.h>

int
main (int argc, char *argv[])
{
  size_t numtests = 10000;
  NTL::SetSeed (NTL::ZZ (0));
  HPR hpr;

  std::vector<HPR_ZZq> xs(numtests);
  std::vector<HPR_ZZq> ys(numtests);
  std::vector<HPR_ZZq> zs(numtests);
  std::vector<NTL::ZZ> exp_zs(numtests);

  for (size_t i = 0; i < numtests; i++) {
    NTL::ZZ x (RandomBnd (P));
    auto xi = hpr.to_ZZq (x);
    assert (hpr.to_ZZ (xi) == x);    
  }

  for (size_t i = 0; i < numtests; i++) {
    NTL::ZZ x (RandomBnd (P));
    NTL::ZZ y (RandomBnd (P));
    NTL::ZZ exp_z = (x*y) % P;

    xs[i] = hpr.to_ZZq (x);
    ys[i] = hpr.to_ZZq (y);
    exp_zs[i] = exp_z;
  }

  //warmup
  for (size_t i = 0; i < numtests; i++) {
    hpr.mul (zs[i], xs[i], ys[i]);
  }

  auto start = std::chrono::high_resolution_clock::now();
  for (size_t i = 0; i < numtests; i++) {
    hpr.mul (zs[i], xs[i], ys[i]);
  }
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;

  for (size_t i = 0; i < numtests; i++) {
    assert (hpr.to_ZZ (zs[i]) == exp_zs[i]);
  }

  std::cout << "average execution time = " << elapsed_seconds.count () / numtests << "\n";
}
