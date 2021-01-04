#include <stdexcept>

#include "grammar_compressed_index.hpp"

int main (int argc, char **argv)
{
  if (argc != 2)
  {
    throw std::runtime_error(std::string("usage: ") + argv[0] + " [input path]");
  }

<<<<<<< HEAD
  sdsl::csa_wt<> fm_index;
  std::ifstream fm_index_input {"../input/index/fm_index_" + gci::util::basename(argv[1])};
  sdsl::load(fm_index, fm_index_input);
=======
  sdsl::int_vector<8> text;
  sdsl::load_vector_from_file(text, argv[1], 1);
  sdsl::append_zero_symbol(text);
>>>>>>> parent of 10fb762... driver for benchmarking

  gc_index gci(text);

  return 0;
}
