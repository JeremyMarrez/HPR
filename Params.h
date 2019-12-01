#include <NTL/ZZ.h>
#include <stddef.h>
#include <stdint.h>
#include <array>

#define P384

#ifdef P512
typedef uint32_t value_t;
typedef uint64_t greater_value_t;
  
constexpr size_t n = 4;
constexpr size_t h1 = 4;
constexpr size_t h2 = 4;
constexpr size_t k = 10;
constexpr size_t beta = 2;
constexpr size_t logw = 32;
  
const NTL::ZZ P (NTL::INIT_VAL, "13407807330566887969602960805627432681041964773857428423494554478465229832259908929071397198063729901922606423888174985480069630617194830561296715100750623");
const NTL::ZZ B1 (NTL::INIT_VAL, "340282363117986678532205813233593895785");

const std::array<value_t, h1> b1 = {4294967287, 4294967285, 4294967283, 4294967281};
const std::array<value_t, h2> b2 = {4294967293, 4294967291, 4294967279, 4294967273};
constexpr size_t log_bsk = 32;

/*
  num_ops = 383
*/
#endif

#ifdef P448
typedef uint32_t value_t;
typedef uint64_t greater_value_t;
  
constexpr size_t n = 2;
constexpr size_t h1 = 7;
constexpr size_t h2 = 7;
constexpr size_t k = 10;
constexpr size_t beta = 46;
constexpr size_t logw = 32;
  
const NTL::ZZ P (NTL::INIT_VAL, "726838691464923927763584728560397469892191696384669397203351854461062993756785094899260663546322497585003463546081470061602503689214083");
const NTL::ZZ B1 (NTL::INIT_VAL, "26959946058271777034864697610288475726886074446769625910911259768377");

const std::array<value_t, h1> b1 = {4294967293, 4294967291, 4294967287, 4294967281, 4294967279, 4294967273, 4294967271};
const std::array<value_t, h2> b2 = {4294967275, 4294967269, 4294967267, 4294967263, 4294967251, 4294967249, 4294967239};
constexpr size_t log_bsk = 32;

/*
  num_ops = 338
*/
#endif

#ifdef P384
typedef uint32_t value_t;
typedef uint64_t greater_value_t;
  
constexpr size_t n = 2;
constexpr size_t h1 = 6;
constexpr size_t h2 = 6;
constexpr size_t k = 10;
constexpr size_t beta = 28;
constexpr size_t logw = 32;
  
const NTL::ZZ P (NTL::INIT_VAL, "39402004288203672800299513266296106078442387709487952440882863429633482730313686776512926555284932454571101674530597");
const NTL::ZZ B1 (NTL::INIT_VAL, "6277101583390511993976998989794444305009475777960744148575");

const std::array<value_t, h1> b1 = {4294967287, 4294967281, 4294967279, 4294967277, 4294967275, 4294967273};
const std::array<value_t, h2> b2 = {4294967293, 4294967291, 4294967269, 4294967267, 4294967263, 4294967257};
constexpr size_t log_bsk = 32;

/*
  num_ops = 267
*/
#endif
