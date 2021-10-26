#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>

#include "index.h"
#include "test.h"
#include "utility.h"

int main ()
{
  auto compressed_corpus_path {std::filesystem::path{"../data/compressed_corpus"}};
  auto corpus_path {std::filesystem::path{"../data/corpus"}};
  figiss::DecompressCompressedCorpus(compressed_corpus_path, corpus_path);
  std::cout << "test counting of index built from text in corpus (" << std::filesystem::canonical(corpus_path).string() << ")\n";
  for (auto const& entry : std::filesystem::directory_iterator(corpus_path))
  {
    if (entry.is_regular_file())
    {
      auto byte_text_path {entry.path()};
      auto index_path {std::filesystem::path{"../data/index/ours/" + entry.path().filename().string() + ".i"}};
      auto pattern_parent_path {std::filesystem::path{"../data/pattern"}};
      figiss::Index<> index {byte_text_path};
      index.Serialize(index_path);
      figiss::TestCounting(byte_text_path, index_path, pattern_parent_path);
    }
  }
  return 0;
}
