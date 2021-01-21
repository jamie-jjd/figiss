#include <fstream>
#include <iostream>
#include <stdexcept>

#include "benchmark.hpp"

int main (int argc, char **argv)
{
  if (argc != 2)
  {
    throw std::runtime_error(std::string("usage: ") + argv[0] + " [input path]");
  }
  gci::test_count(argv[1]);
  return 0;
}
