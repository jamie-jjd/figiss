#include <stdexcept>

#include "grammar_compressed_index.hpp"

int main (int argc, char **argv)
{
  if (argc != 2)
  {
    throw std::runtime_error(std::string("usage: ") + argv[0] + " [input path]");
  }

  sdsl::int_vector<8> text;
  sdsl::load_vector_from_file(text, argv[1], 1);
  sdsl::append_zero_symbol(text);

  gc_index gci(text);

  return 0;
}
