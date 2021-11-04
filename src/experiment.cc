#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>

#include "experiment.h"

int main ()
{
  auto compressed_corpus_path {std::filesystem::path{"../data/compressed_corpus"}};
  auto corpus_path {std::filesystem::path{"../data/corpus"}};
  figiss::DecompressCompressedCorpus(compressed_corpus_path, corpus_path);
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
        for (uint8_t max_factor_size {1}; max_factor_size <= 8; ++max_factor_size)
        {
          figiss::ConstructAndSerialize<figiss::Index>(byte_text_path, max_factor_size, "figiss");
        }
        figiss::ConstructAndSerialize<sdsl::csa_wt<wt_fbb<>, 0xFFFF'FFFF, 0xFFFF'FFFF>>(byte_text_path, "faster-minuter");
        figiss::ConstructAndSerialize<sdsl::csa_wt<sdsl::wt_rlmn<>, 0xFFFF'FFFF, 0xFFFF'FFFF>>(byte_text_path, "rlfm");
        figiss::ConstructAndSerialize<CSA::RLCSA>(byte_text_path, "rlcsa");
      }
      uint64_t amount {1ULL << 12};
      for (auto length {1ULL << 0}; length <= (1ULL << 15); length <<= 1)
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
          for (uint8_t max_factor_size {1}; max_factor_size <= 8; ++max_factor_size)
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
          std::string const& index_name {"faster-minuter"};
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
      }
    }
  }
  return 0;
}
