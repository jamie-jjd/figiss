#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>

#include "index.h"

void Usage ()
{
  std::cout << "usage:\n"
  << "./figiss cs [text path] [index path]\n"
  << "\tconstruct index of (ASCII-encoded) text at [text path] and serialize it to [index path]\n"
  << "./figiss lc [index path] [pattern path]\n"
  << "\tload index from [index path] and report number of occurences of (ACSII-encoded) pattern at [pattern path]\n";
  return;
}

int main (int argc, char **argv)
{
  if (argc != 4)
  {
    Usage();
  }
  else
  {
    auto option {std::string(argv[1])};
    if (option == "cs")
    {
      auto byte_text_path {std::filesystem::path(argv[2])};
      auto index_path {std::filesystem::path(argv[3])};
      figiss::Index<> index {byte_text_path};
      index.Serialize(index_path);
    }
    else if (option == "lc")
    {
      auto index_path {std::filesystem::path(argv[2])};
      auto pattern_path {std::filesystem::path(argv[3])};
      sdsl::int_vector<8> pattern;
      sdsl::load_vector_from_file(pattern, pattern_path);
      figiss::Index<> index;
      index.Load(index_path);
      std::cout << index.Count(std::begin(pattern), std::end(pattern)) << "\n";
    }
    else
    {
      Usage();
    }
  }
  return 0;
}
