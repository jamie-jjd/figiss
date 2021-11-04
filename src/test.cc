#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>

#include "index.h"
#include "test.h"
#include "utility.h"

void PrintUsage ()
{
  std::cout
  << "./test \u03BB_l \u03BB_r b p_l p_r\n"
  << "\ttest counting of figiss,\n"
  << "\twhere \u03BB is from \u03BB_l to \u03BB_r and 1 <= \u03BB_l <= \u03BB_r <= 8\n"
  << "\tfor each p from p_l to p_r,\n"
  << "\tpattern collection of b random patterns of length 2^p is generated\n";
}

int main (int argc, char** argv)
{
  if (argc != 6)
  {
    PrintUsage();
  }
  else
  {
    auto lower_max_factor_size {static_cast<uint8_t>(std::stoull(argv[1]))};
    auto upper_max_factor_size {static_cast<uint8_t>(std::stoull(argv[2]))};
    auto amount {std::stoull(argv[3])};
    auto lower_power {std::stoull(argv[4])};
    auto upper_power {std::stoull(argv[5])};
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
        for (uint8_t max_factor_size {lower_max_factor_size}; max_factor_size <= upper_max_factor_size; ++max_factor_size)
        {
          auto index_path {std::filesystem::path {"../data/index/figiss/" + std::to_string(max_factor_size) + "/index"}};
          figiss::Index index {byte_text_path, max_factor_size};
          index.Serialize(index_path);
          figiss::TestCounting(byte_text_path, index_path, pattern_parent_path, amount, (1ULL << lower_power), (1ULL << upper_power));
        }
      }
    }
  }
  return 0;
}
