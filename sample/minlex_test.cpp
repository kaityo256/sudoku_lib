#include "../minlex.hpp"
#include <fstream>
#include <iostream>

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
    std::cout << "Usage: ./minlex_test filename" << std::endl;
  }
}
