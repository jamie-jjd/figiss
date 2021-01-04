#include <fstream>
#include <iostream>
#include <stdexcept>

#include <sdsl/suffix_trees.hpp>

#include "benchmark.hpp"
#include "grammar_compressed_index.hpp"
#include "utility.hpp"

int main (int argc, char **argv)
{
  if (argc != 2)
  {
    throw std::runtime_error
    (
      std::string("usage: ")
      + argv[0]
      + " [text filename]"
    );
  }
  gci::gc_index index;
  std::string index_path {"../input/index/gc_index_" + gci::util::basename(argv[1])};
  gci::construct(index, argv[1]);

  gci::load(index, index_path);
  sdsl::csa_wt<> fm_index;
  std::ifstream fm_index_input {"../input/index/fm_index_" + gci::util::basename(argv[1])};
  sdsl::load(fm_index, fm_index_input);

  gci::benchmark_count(index, fm_index, argv[1], 1000, 20000);

  return 0;
}
