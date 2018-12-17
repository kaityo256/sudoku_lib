#include "../ansmaker.hpp"
#include <iostream>

int main() {
  // Making Answer
  std::cout << "Making Answer" << std::endl;
  AnsMaker am;
  for (int i = 0; i < 10; i++) {
    std::cout << am.make() << std::endl;
  }
}
