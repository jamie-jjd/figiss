#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>

#include "experiment.h"

void PrintUsage ()
{
  std::cout
  << "./experiment \u03BB_l \u03BB_r b p_l p_r\n"
  << "\trun experiment on figiss\n"
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
    auto compressed_corpus_path {std::filesystem::path{"../data/compressed_corpus"}};
    auto corpus_path {std::filesystem::path{"../data/corpus"}};
    figiss::DecompressCompressedCorpus(compressed_corpus_path, corpus_path);
    auto lower_max_factor_size {static_cast<uint8_t>(std::stoull(argv[1]))};
    auto upper_max_factor_size {static_cast<uint8_t>(std::stoull(argv[2]))};
    auto amount {std::stoull(argv[3])};
    auto lower_power {std::stoull(argv[4])};
    auto upper_power {std::stoull(argv[5])};
    for (auto const& entry : std::filesystem::directory_iterator(corpus_path))
    {
      if (entry.is_regular_file())
      {
        auto byte_text_path {entry.path()};
        std::cout << "\n\trun experiment w.r.t. " << std::filesystem::canonical(byte_text_path).string() << "\n\n";
        {
          std::cout << "--- calculate text parameters ---\n";
          figiss::PrintTextParameters(byte_text_path);
        }
        {
          for (uint8_t max_factor_size {lower_max_factor_size}; max_factor_size <= upper_max_factor_size; ++max_factor_size)
          {
            figiss::ConstructAndSerialize<figiss::Index>(byte_text_path, max_factor_size, "figiss");
          }
          figiss::ConstructAndSerialize<sdsl::csa_wt<wt_fbb<>, 0xFFFF'FFFF, 0xFFFF'FFFF>>(byte_text_path, "fm_fbb");
          figiss::ConstructAndSerialize<sdsl::csa_wt<sdsl::wt_rlmn<>, 0xFFFF'FFFF, 0xFFFF'FFFF>>(byte_text_path, "rlfm");
          figiss::ConstructAndSerialize<CSA::RLCSA>(byte_text_path, "rlcsa");
          figiss::ConstructAndSerialize<sdsl::csa_wt<sdsl::wt_blcd<>, 0xFFFF'FFFF, 0xFFFF'FFFF>>(byte_text_path, "fm");
          figiss::ConstructAndSerialize<lz77index::static_selfindex_lzend>(byte_text_path, "lz_end");
        }
        for (auto length {1ULL << lower_power}; length <= (1ULL << upper_power); length <<= 1)
        {
          auto pattern_path
          {
            std::filesystem::path
            {
              "../data/pattern/" +
              byte_text_path.filename().string() + "/" +
              std::to_string(amount) + "/" +
              std::to_string(length) + "/" +
              "pattern"
            }
          };
          figiss::PatternCollection pattern_collection;
          if (!std::filesystem::exists(pattern_path))
          {
            std::cout << "--- generate & serialize pattern collection ---\n";
            pattern_collection = figiss::PatternCollection{byte_text_path, "../data/pattern", amount, length};
            pattern_collection.Serialize(pattern_path);
          }
          else
          {
            std::cout << "--- load pattern collection ---\n";
            pattern_collection.Load(pattern_path);
          }
          {
            std::string const& index_name {"figiss"};
            for (uint8_t max_factor_size {lower_max_factor_size}; max_factor_size <= upper_max_factor_size; ++max_factor_size)
            {
              std::cout << "--- measure " << index_name << " counting time (\u03BB = " << std::to_string(max_factor_size) << ") ---\n";
              auto middle_path {byte_text_path.filename().string() + "/" + index_name + "/" + std::to_string(max_factor_size)};
              auto middle_pattern_path {middle_path / std::filesystem::path{std::to_string(amount) + "/" + std::to_string(length)}};
              auto index_path {std::filesystem::path{"../data/index"} / middle_path / std::filesystem::path{"index"}};
              figiss::Index index;
              index.Load(index_path);
              auto counting_time_path {std::filesystem::path{"../data/counting_time"} / middle_pattern_path / std::filesystem::path{"counting_time"}};
              figiss::MeasureCountingTime(pattern_path, index, counting_time_path);
            }
          }
          {
            std::string const& index_name {"fm_fbb"};
            auto middle_path {byte_text_path.filename().string() + "/" + index_name};
            auto middle_pattern_path {middle_path / std::filesystem::path{std::to_string(amount) + "/" + std::to_string(length)}};
            std::cout << "--- measure " << index_name << " counting time ---\n";
            auto index_path {std::filesystem::path{"../data/index"} / middle_path / std::filesystem::path{"index"}};
            std::ifstream in {index_path};
            std::cout << "load " << index_name << " from " << std::filesystem::canonical(index_path).string() << "\n";
            sdsl::csa_wt<wt_fbb<>, 0xFFFF'FFFF, 0xFFFF'FFFF> index;
            index.load(in);
            {
              auto counting_time_path {std::filesystem::path{"../data/counting_time"} / middle_pattern_path / std::filesystem::path{"counting_time"}};
              figiss::MeasureCountingTime(pattern_path, index, counting_time_path);
            }
          }
          {
            std::string const& index_name {"rlfm"};
            auto middle_path {byte_text_path.filename().string() + "/" + index_name};
            auto middle_pattern_path {middle_path / std::filesystem::path{std::to_string(amount) + "/" + std::to_string(length)}};
            std::cout << "--- measure " << index_name << " counting time ---\n";
            auto index_path {std::filesystem::path{"../data/index"} / middle_path / std::filesystem::path{"index"}};
            std::ifstream in {index_path};
            std::cout << "load " << index_name << " from " << std::filesystem::canonical(index_path).string() << "\n";
            sdsl::csa_wt<sdsl::wt_rlmn<>, 0xFFFF'FFFF, 0xFFFF'FFFF> index;
            index.load(in);
            {
              auto counting_time_path {std::filesystem::path{"../data/counting_time"} / middle_pattern_path / std::filesystem::path{"counting_time"}};
              figiss::MeasureCountingTime(pattern_path, index, counting_time_path);
            }
          }
          {
            std::string const& index_name {"rlcsa"};
            auto middle_path {byte_text_path.filename().string() + "/" + index_name};
            auto middle_pattern_path {middle_path / std::filesystem::path{std::to_string(amount) + "/" + std::to_string(length)}};
            std::cout << "--- measure " << index_name << " counting time ---\n";
            auto index_subpath {std::filesystem::path{"../data/index"} / middle_path / byte_text_path.filename()};
            std::cout << "load " << index_name << " from " << std::filesystem::canonical(std::filesystem::path{index_subpath.string() + ".rlcsa.array"}).string() << "\n";
            CSA::RLCSA index(index_subpath.string());
            auto counting_time_path {std::filesystem::path{"../data/counting_time"} / middle_pattern_path / std::filesystem::path{"counting_time"}};
            figiss::MeasureCountingTime(pattern_path, index, counting_time_path);
          }
          {
            std::string const& index_name {"fm"};
            auto middle_path {byte_text_path.filename().string() + "/" + index_name};
            auto middle_pattern_path {middle_path / std::filesystem::path{std::to_string(amount) + "/" + std::to_string(length)}};
            std::cout << "--- measure " << index_name << " counting time ---\n";
            auto index_path {std::filesystem::path{"../data/index"} / middle_path / std::filesystem::path{"index"}};
            std::ifstream in {index_path};
            std::cout << "load " << index_name << " from " << std::filesystem::canonical(index_path).string() << "\n";
            sdsl::csa_wt<sdsl::wt_blcd<>, 0xFFFF'FFFF, 0xFFFF'FFFF> index;
            index.load(in);
            {
              auto counting_time_path {std::filesystem::path{"../data/counting_time"} / middle_pattern_path / std::filesystem::path{"counting_time"}};
              figiss::MeasureCountingTime(pattern_path, index, counting_time_path);
            }
          }
          {
            std::string const& index_name {"lz_end"};
            auto middle_path {byte_text_path.filename().string() + "/" + index_name};
            auto middle_pattern_path {middle_path / std::filesystem::path{std::to_string(amount) + "/" + std::to_string(length)}};
            std::cout << "--- measure " << index_name << " counting time ---\n";
            auto index_path {std::filesystem::path{"../data/index"} / middle_path / std::filesystem::path{"index"}};
            std::cout << "load " << index_name << " from " << std::filesystem::canonical(index_path).string() << "\n";
            auto* index {lz77index::static_selfindex::load(index_path.string().c_str())};
            {
              auto counting_time_path {std::filesystem::path{"../data/counting_time"} / middle_pattern_path / std::filesystem::path{"counting_time"}};
              figiss::MeasureCountingTime(pattern_path, index, counting_time_path);
            }
          }
        }
      }
    }
  }
  return 0;
}
