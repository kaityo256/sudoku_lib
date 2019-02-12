#pragma once
#include <iostream>
#include <string>
#ifdef _WIN32
#include <intrin.h>
#else
#include <x86intrin.h>
#endif

typedef unsigned __int128 mbit;

const static mbit mask81 = (mbit(1) << 81) - 1;
inline std::ostream &operator<<(std::ostream &dest, mbit v) {
  char buffer[81];
  for (int i = 0; i < 81; i++) {
    buffer[i] = (v & 1) + '0';
    v = v >> 1;
  }
  dest.rdbuf()->sputn(buffer, 81);
  return dest;
}

inline mbit str2mbit(const std::string &str) {
  mbit t = 0;
  for (unsigned int i = 0; i < 81 && i < str.size(); i++) {
    t ^= mbit((str[i] != '0')) << i;
  }
  return t;
}

inline std::string mbit2str(mbit m, const std::string &str) {
  assert(str.length() == 81);
  std::string str2 = str;
  for (int i = 0; i < 81; i++) {
    mbit v = 1;
    v = v << i;
    if ((m & v) == 0) str2[i] = '0';
  }
  return str2;
}

inline int popcnt_u128(const mbit &s) {
  const uint64_t *m = (uint64_t *)&s;
  return _mm_popcnt_u64(m[0]) + _mm_popcnt_u64(m[1]);
}

inline int bitpos(mbit s) {
  return popcnt_u128(s - 1);
}
