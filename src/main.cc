#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>

#include "benchmark.h"
#include "grammar_compressed_index.h"

int main (int argc, char **argv)
{
  if (argc != 2)
  {
    throw std::runtime_error(std::string("usage: ") + argv[0] + " [text path]");
  }
  std::filesystem::path text_path {argv[1]};
  {
    project::Index index;
    project::PrintSpace(index, text_path);
  }
  {
    project::Index index;
    project::TestCount(index, text_path);
  }
  return 0;
}
