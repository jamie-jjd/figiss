#include <stdexcept>
#include <string>
#include <sdsl/suffix_trees.hpp>
#include "grammar_compressed_index.hpp"

int main (int argc, char **argv)
{
  if (argc != 2)
  {
    throw std::runtime_error(std::string("usage: ") + argv[0] + " [input path]");
  }

  Grammar_Compressed_Index gci {argv[1]};

  return 0;
}
