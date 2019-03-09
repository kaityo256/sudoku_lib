#pragma once
#include "grid.hpp"
#include <algorithm>
#include <random>
class AnsMaker {
private:
  bool search_answer(Grid &g2, const int index, std::string &ans) {
    if (index == 81) {
      ans = g2.str();
      return true;
    }
    int a[] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    std::shuffle(a, a + 9, mt);
    for (unsigned int i = 0; i < 9; i++) {
      Grid g(g2);
      if (!g.can_put(index, a[i]))
        continue;
      g.put(index, a[i]);
      if (search_answer(g, index + 1, ans))
        return true;
    }
    return false;
  }
  std::mt19937 mt;

public:
  AnsMaker() : mt(1) {}
  AnsMaker(const int seed) : mt(seed) {}
  std::string make(void) {
    Grid g;
    std::string ans;
    search_answer(g, 0, ans);
    return ans;
  }
};
