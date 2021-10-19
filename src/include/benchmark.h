#pragma once

#include <chrono>

#include <sdsl/csa_wt.hpp>
#include <sdsl/suffix_array_algorithm.hpp>

#include "index.h"
#include "pattern_collection.h"
#include "utility.h"

namespace figiss
{
void TestCounting
(
  std::filesystem::path const& byte_text_path,
  std::filesystem::path const& index_path,
  std::filesystem::path& pattern_parent_path,
  uint64_t const amount = 4096,
  uint64_t const min_length = 1,
  uint64_t const max_length = 32768
)
{
  sdsl::csa_wt<sdsl::wt_rlmn<>, 0xFFFF'FFFF, 0xFFFF'FFFF> rlfm;
  {
    sdsl::int_vector<8> byte_text;
    sdsl::load_vector_from_file(byte_text, byte_text_path);
    sdsl::construct_im(rlfm, byte_text);
  }
  Index index;
  {
    index.Load(index_path);
  }
  {
    for (uint64_t length {min_length}; length <= max_length; length <<= 1)
    {
      PatternCollection patterns;
      for (uint8_t is_mutated {}; is_mutated != 2; ++is_mutated)
      {
        auto metadata_path
        {
          pattern_parent_path
          / std::filesystem::path{byte_text_path.filename().string()}
          / std::filesystem::path{std::to_string(length) + ".meta"}
        };
        if (is_mutated)
        {
          metadata_path.replace_extension(".mutated.meta");
        }
        patterns = PatternCollection(byte_text_path, metadata_path, amount, length, is_mutated);
        if (patterns)
        {
          // auto pattern_path
          // {
          //   pattern_parent_path
          //   / std::filesystem::path{byte_text_path.filename().string()}
          //   / std::filesystem::path{std::to_string(length) + ".pattern"}
          // };
          if (is_mutated)
          {
            // metadata_path.replace_extension(".mutated.pattern");
            std::cout << "test " << amount << " random \"mutated\" existed patterns of length " << length << "\n";
          }
          else
          {
            std::cout << "test " << amount << " random existed patterns of length " << length << "\n";
          }
          // patterns.Serialize(pattern_path);
          auto begin {std::begin(patterns)};
          auto end {begin};
          for (uint64_t i {}; i != patterns.GetAmount(); ++i)
          {
            begin = end;
            end = std::next(begin, patterns.GetLength());
            auto rlfm_count {sdsl::count(rlfm, begin, end)};
            auto index_count {index.Count(begin, end)};
            if (rlfm_count != index_count)
            {
              std::string pattern {begin, end};
              throw std::runtime_error
              (
                "\033[31mfailed at pattern:" + pattern
                + ",genuine count:" + std::to_string(rlfm_count)
                + ",count:" + std::to_string(index_count)
                + "\033[0m");
            }
          }
          std::cout << "\033[32msuccess\033[0m\n";
        }
        else { return; }
      }
    }
  }
  return;
}

void PrintCountingTime
(
  std::filesystem::path const& output_path,
  std::vector<std::pair<uint64_t, double>> time,
  bool const is_proper_representation = false
)
{
  if (!std::filesystem::exists(output_path.parent_path()))
  {
    std::filesystem::create_directories(output_path.parent_path());
  }
  std::fstream fout {output_path, std::ios_base::out | std::ios_base::trunc};
  std::cout << "write time information to " << std::filesystem::canonical(output_path) << "\n";
  for (auto const& pair : time)
  {
    if (is_proper_representation)
    {
      fout << std::fixed << std::setprecision(2);
      fout << std::get<0>(pair) << "," << ProperTimeRepresentation(std::get<1>(pair)) << "\n";
    }
    else
    {
      fout << std::get<0>(pair) << "," << std::get<1>(pair) << "\n";
    }
  }
  return;
}

// template <typename Index>
// void MeasureIndexCountingTime
// (
//   std::filesystem::path const& byte_text_path,
//   Index& index,
//   bool const is_proper_representation = false
// )
// {
//   auto index_path
//   {
//     std::filesystem::path
//     (
//       "../data/index/"
//       + std::to_string(Index::kMaxFactorSize) + "/"
//       + byte_text_path.filename().string()
//       + ".gci"
//     )
//   };
//   if (!std::filesystem::exists(index_path))
//   {
//     index = Index{byte_text_path};
//     index.Serialize(index_path);
//   }
//   std::vector<std::pair<uint64_t, double>> time;
//   uint64_t amount {1ULL << 12};
//   for (uint64_t length {1ULL << 8}; length <= (1ULL << 15); length <<= 1)
//   {
//     std::cout << "pattern size: " << length << "\n";
//     auto patterns {PatternCollection(byte_text_path, amount, length)};
//     index.Load(index_path);
//     {
//       auto begin {std::begin(patterns)};
//       while (begin != std::end(patterns))
//       {
//         auto end {std::next(begin, patterns.GetLength())};
//         index.Count(begin, end);
//         begin = end;
//       }
//     }
//     {
//       auto begin {std::begin(patterns)};
//       auto begin_time {std::chrono::steady_clock::now()};
//       while (begin != std::end(patterns))
//       {
//         auto end {std::next(begin, patterns.GetLength())};
//         index.Count(begin, end);
//         begin = end;
//       }
//       auto duration {static_cast<double>((std::chrono::steady_clock::now() - begin_time).count())};
//       time.emplace_back(patterns.GetLength(), duration / patterns.GetAmount());
//     }
//   }
//   PrintCountingTime
//   (
//     std::filesystem::path
//     {
//       std::string("../data/time/")
//       + std::to_string(Index::kMaxFactorSize) + "/"
//       + byte_text_path.filename().string()
//     },
//     time,
//     is_proper_representation
//   );
//   return;
// }
//
// void MeasureSdslIndexCountingTime
// (
//   std::filesystem::path const& byte_text_path,
//   bool const is_proper_representation = false
// )
// {
//   auto index_path
//   {
//     std::filesystem::path
//     (
//       "../data/index/rlfm/"
//       + byte_text_path.filename().string()
//       + ".rlfm"
//     )
//   };
//   sdsl::csa_wt<sdsl::wt_rlmn<>, 0xFFFF'FFFF, 0xFFFF'FFFF> rlfm;
//   if (!std::filesystem::exists(index_path))
//   {
//     std::filesystem::create_directories(index_path.parent_path());
//     sdsl::int_vector<8> byte_text;
//     sdsl::load_vector_from_file(byte_text, byte_text_path);
//     std::cout << "construct rlfm of " << std::filesystem::canonical(byte_text_path) << "\n";
//     sdsl::construct_im(rlfm, byte_text);
//     std::ofstream fout {index_path, std::ios_base::out | std::ios_base::trunc};
//     sdsl::serialize(rlfm, fout);
//     std::cout << "serialize rlfm to " << std::filesystem::canonical(index_path) << "\n";
//   }
//   std::vector<std::pair<uint64_t, double>> time;
//   uint64_t amount {1ULL << 12};
//   for (uint64_t length {1ULL << 8}; length <= (1ULL << 15); length <<= 1)
//   {
//     std::cout << "pattern size: " << length << "\n";
//     auto patterns {PatternCollection(byte_text_path, amount, length)};
//     {
//       std::ifstream in {index_path};
//       rlfm.load(in);
//       std::cout << "load rlfm from " << std::filesystem::canonical(index_path) << "\n";
//     }
//     {
//       auto begin {std::begin(patterns)};
//       while (begin != std::end(patterns))
//       {
//         auto end {std::next(begin, patterns.GetLength())};
//         sdsl::count(rlfm, begin, end);
//         begin = end;
//       }
//     }
//     {
//       auto begin {std::begin(patterns)};
//       auto begin_time {std::chrono::steady_clock::now()};
//       while (begin != std::end(patterns))
//       {
//         auto end {std::next(begin, patterns.GetLength())};
//         sdsl::count(rlfm, begin, end);
//         begin = end;
//       }
//       auto duration {static_cast<double>((std::chrono::steady_clock::now() - begin_time).count())};
//       time.emplace_back(patterns.GetLength(), duration / patterns.GetAmount());
//     }
//   }
//   PrintCountingTime
//   (
//     std::filesystem::path
//     {
//       std::string("../data/time/")
//       + "rlfm/"
//       + byte_text_path.filename().string()
//     },
//     time,
//     is_proper_representation
//   );
//   return;
// }
}
