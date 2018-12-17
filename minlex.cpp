#include "minlex.hpp"
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <functional>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

void test(void) {
  std::string str =
      "207005000000340000150000009005000001040000320000016500000002084700000010010580000";
  std::string ans =
      "000000012000034005006007300001300007053080000080000100010005090200100000700400030";
  MinlexSearcher s;
  std::string min = s.search(str);
  std::cout << min << std::endl;
  if (min == ans) {
    std::cout << "OK" << std::endl;
  } else {
    std::cout << "NG" << std::endl;
  }
}

void input(const char *filename) {
  std::string line;
  MinlexSearcher s;
  std::ifstream ifs(filename);
  while (getline(ifs, line)) {
    std::cout << s.search(line.c_str()) << std::endl;
  }
}

int main(int argc, char **argv) {
  if (argc > 1) {
    input(argv[1]);
  } else {
    test();
  }
}
