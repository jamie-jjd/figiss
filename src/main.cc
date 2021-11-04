#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>

#include "index.h"

void PrintUsage ()
{
  std::cout << "usage:\n"
  << "./figiss cs \u03bb \"byte-text path\" \"index path\"\n"
  << "\tconstruct index of (ASCII-encoded) text at \"byte-text path\"\n"
  << "\twith max factor size being \u03bb, which can only be an integer in [1..8]\n"
  << "\tand serialize it to \"index path\"\n"
  << "./figiss lc \"index path\" \"byte-pattern path\"\n"
  << "\tload index from \"index path\" and report number of occurences of (ACSII-encoded) pattern at \"byte-pattern path\"\n";
  return;
}

int main (int argc, char** argv)
{
  if (argc < 4 || argc > 5)
  {
    PrintUsage();
  }
  else
  {
    auto option {std::string(argv[1])};
    if (option == "cs")
    {
      if (argc != 5)
      {
        PrintUsage();
      }
      else
      {
        auto max_factor_size {std::stoull(argv[2])};
        auto byte_text_path {std::filesystem::path(argv[3])};
        auto index_path {std::filesystem::path(argv[4])};
        figiss::Index index {byte_text_path, static_cast<uint8_t>(max_factor_size)};
        index.Serialize(index_path);
      }
    }
    else if (option == "lc")
    {
      if (argc != 4)
      {
        PrintUsage();
      }
      else
      {
        auto index_path {std::filesystem::path(argv[2])};
        auto pattern_path {std::filesystem::path(argv[3])};
        sdsl::int_vector<8> pattern;
        sdsl::load_vector_from_file(pattern, pattern_path);
        figiss::Index index;
        index.Load(index_path);
        std::cout << index.Count(std::begin(pattern), std::end(pattern)) << "\n";
      }
    }
    else
    {
      PrintUsage();
    }
  }
  return 0;
}
