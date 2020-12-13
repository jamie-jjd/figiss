#include <algorithm>
#include <fstream>
#include <iostream>
#include <numeric>

#include <sdsl/construct.hpp>
// #include "parse_text.hpp"
// #include "grammar_compressed_index.hpp"
#include "sl_typed_string.hpp"

template <typename T>
void f (T)
{
  std::cout << __PRETTY_FUNCTION__ << '\n';
}

int main (int argc, char **argv)
{
  if (argc > 0 && argv[0][0])
  {
    sdsl::int_vector<8> text;
    sdsl::load_vector_from_file(text, argv[1], 1);
    sdsl::append_zero_symbol(text);
    SL_Typed_String SL_typed_string {text};
    std::cout << SL_typed_string << '\n';
  }
  return 0;
}
