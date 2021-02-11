#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>

#include "benchmark.hpp"

namespace fs = std::filesystem;

int main (int argc, char **argv)
{
  if (argc != 2)
  {
    throw std::runtime_error(std::string("usage: ") + argv[0] + " [text path]");
  }
  fs::path text_path {argv[1]};
  // gci::generate_canonical_text(argv[1]);
  // gci::print_text_attributes(argv[1]);
  // gci::test_count(argv[1]);
  // gci::benchmark_count(argv[1], 100, 50000);
  // gci::benchmark_count(argv[1], 100, 100000);
  // gci::benchmark_count(argv[1], 100, 200000);
  return 0;
}
