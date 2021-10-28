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
      auto index_parent_path {std::filesystem::path{"../data/index"}};
      if (!std::filesystem::exists(index_parent_path))
      {
        std::filesystem::create_directories(index_parent_path);
      }
      {
        auto construction_parent_path {std::filesystem::path{"../data/construction"}};
        if (!std::filesystem::exists(construction_parent_path))
        {
          std::filesystem::create_directories(construction_parent_path);
        }
        figiss::ConstructAndSerialize<sdsl::csa_wt<sdsl::wt_rlmn<>, 0xFFFF'FFFF, 0xFFFF'FFFF>>(byte_text_path, construction_parent_path, index_parent_path);
        figiss::ConstructAndSerialize<sdsl::csa_wt<wt_fbb<>, 0xFFFF'FFFF, 0xFFFF'FFFF>>(byte_text_path, construction_parent_path, index_parent_path);
        figiss::ConstructAndSerialize<CSA::RLCSA>(byte_text_path, construction_parent_path, index_parent_path);
        figiss::ConstructAndSerialize<figiss::Index<>>(byte_text_path, construction_parent_path, index_parent_path);
        auto index {figiss::Index<>{byte_text_path, true}};
      }
      auto pattern_parent_path {std::filesystem::path{"../data/pattern"}};
      if (!std::filesystem::exists(pattern_parent_path))
      {
        std::filesystem::create_directories(pattern_parent_path);
      }
      auto time_parent_path {std::filesystem::path{"../data/counting"}};
      if (!std::filesystem::exists(time_parent_path))
      {
        std::filesystem::create_directories(time_parent_path);
      }
      uint64_t amount {1ULL << 12};
      for (auto length {1ULL << 0}; length <= (1ULL << 15); length <<= 1)
      {
        auto pattern_path
        {
          pattern_parent_path /
          byte_text_path.filename() /
          std::filesystem::path
          {
            std::to_string(amount) +
            "_" + std::to_string(length) +
            ".pattern"
          }
        };
        figiss::PatternCollection pattern_collection;
        if (!std::filesystem::exists(pattern_path))
        {
          std::cout << "--- generate & serialize pattern collection ---\n";
          pattern_collection = figiss::PatternCollection{byte_text_path, pattern_parent_path, amount, length};
          pattern_collection.Serialize(pattern_path);
        }
        else
        {
          std::cout << "--- load pattern collection ---\n";
          pattern_collection.Load(pattern_path);
        }
        {
          std::cout << "--- measure rlfm counting time ---\n";
          auto index_path {index_parent_path / std::filesystem::path{byte_text_path.filename().string() + ".rlfm"}};
          std::ifstream in {index_path};
          std::cout << "load rlfm from " << std::filesystem::canonical(index_path).string() << "\n";
          sdsl::csa_wt<sdsl::wt_rlmn<>, 0xFFFF'FFFF, 0xFFFF'FFFF> index;
          index.load(in);
          {
            auto time_path
            {
              time_parent_path /
              byte_text_path.filename() /
              std::filesystem::path
              {
                std::to_string(amount) +
                "_" + std::to_string(length) +
                ".rlfm.counting.time"
              }
            };
            figiss::MeasureCountingTime(pattern_path, index, time_path);
          }
        }
        {
          std::cout << "--- measure faster-minuter counting time ---\n";
          auto index_path {index_parent_path / std::filesystem::path{byte_text_path.filename().string() + ".faster-minuter"}};
          std::ifstream in {index_path};
          std::cout << "load faster-minuter from " << std::filesystem::canonical(index_path).string() << "\n";
          sdsl::csa_wt<wt_fbb<>, 0xFFFF'FFFF, 0xFFFF'FFFF> index;
          index.load(in);
          auto time_path
          {
            time_parent_path /
            byte_text_path.filename() /
            std::filesystem::path
            {
              std::to_string(amount) +
              "_" + std::to_string(length) +
              ".faster-minuter.counting.time"
            }
          };
          figiss::MeasureCountingTime(pattern_path, index, time_path);
        }
        {
          std::cout << "--- measure rlcsa counting time ---\n";
          auto index_subpath {index_parent_path / std::filesystem::path{byte_text_path.filename().string()}};
          std::cout << "load rlcsa from " << std::filesystem::canonical(std::filesystem::path{index_subpath.string() + ".rlcsa.array"}).string() << "\n";
          CSA::RLCSA index(index_subpath.string());
          auto time_path
          {
            time_parent_path /
            byte_text_path.filename() /
            std::filesystem::path
            {
              std::to_string(amount) +
              "_" + std::to_string(length) +
              ".rlcsa.counting.time"
            }
          };
          figiss::MeasureCountingTime(pattern_path, index, time_path);
        }
        {
          std::cout << "--- measure figiss counting time ---\n";
          auto index_path {index_parent_path / std::filesystem::path{byte_text_path.filename().string() + ".figiss"}};
          figiss::Index<> index;
          index.Load(index_path);
          auto time_path
          {
            time_parent_path /
            byte_text_path.filename() /
            std::filesystem::path
            {
              std::to_string(amount) +
              "_" + std::to_string(length) +
              ".figiss.counting.time"
            }
          };
          figiss::MeasureCountingTime(pattern_path, index, time_path);
        }
      }
    }
  }
  return 0;
}
