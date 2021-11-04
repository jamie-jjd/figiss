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
      auto pattern_parent_path {std::filesystem::path{"../data/pattern"}};
      for (uint8_t max_factor_size {1}; max_factor_size <= 8; ++max_factor_size)
      {
        auto index_path
        {
          std::filesystem::path
          {
            "../data/index/figiss/" +
            std::to_string(max_factor_size) + "/" +
            entry.path().filename().string()
            + ".figiss"
          }
        };
        figiss::Index index {byte_text_path, max_factor_size};
        index.Serialize(index_path);
        figiss::TestCounting(byte_text_path, index_path, pattern_parent_path, 256, 1, 256);
      }
    }
  }
  return 0;
}
