#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>

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
    project::ConstructIndex(index, text_path);
    project::SerializeIndex(index, "_.sdsl");
  }
  {
    project::Index index;
    project::LoadIndex(index, "_.sdsl");
    sdsl::int_vector<8> text;
    sdsl::load_vector_from_file(text, text_path);
    std::cout << project::Count(index, std::begin(text), std::next(std::begin(text), 10)) << "\n";
  }
  return 0;
}
