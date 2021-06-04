#pragma once

#include <sdsl/csa_wt.hpp>
#include <sdsl/suffix_array_algorithm.hpp>

#include "utility.h"

namespace project
{
// template
// <
//   typename PatternAmount,
//   typename PatternSize,
//   typename Patterns
// >
// void BenchmarkIndexCount
// (
//   Index const &index,
//   PatternAmount const pattern_amount,
//   PatternSize const unit_size,
//   Patterns &patterns
// )
// {
//   auto begin {std::begin(patterns)};
//   auto end {begin};
//   for (uint64_t i {0}; i != pattern_amount; ++i)
//   {
//     begin = end;
//     end = std::next(begin, unit_size);
//     Count(index, begin, end);
//   }
//   return;
// }
//
// template
// <
//   typename Fmindex,
//   typename PatternAmount,
//   typename PatternSize,
//   typename Patterns
// >
// void BenchmarkFmindexCount
// (
//   Fmindex const &rlfm,
//   PatternAmount const pattern_amount,
//   PatternSize const unit_size,
//   Patterns &patterns
// )
// {
//   auto begin {std::begin(patterns)};
//   auto end {begin};
//   for (uint64_t i {0}; i != pattern_amount; ++i)
//   {
//     begin = end;
//     end = std::next(begin, unit_size);
//     sdsl::count(rlfm, begin, end);
//   }
//   return;
// }
//
// void BenchmarkCount
// (
//   std::filesystem::path const &text_path,
//   uint64_t const pattern_amount,
//   uint64_t const unit_size
// )
// {
//   if (unit_size < min_unit_size)
//   {
//     throw std::runtime_error
//     (
//       std::string{"pattern size should be at least "}
//       + std::to_string(min_unit_size)
//       + " (characters) for this text"
//     );
//   }
//   Index index;
//   sdsl::csa_wt<> rlfm;
//   LoadIndexAndFmindex(index, rlfm, text_path);
//   sdsl::int_vector<8> patterns;
//   GeneratePatterns(text_path, pattern_amount, unit_size, patterns);
//   BenchmarkIndexCount(index, pattern_amount, unit_size, patterns);
//   BenchmarkFmindexCount(rlfm, pattern_amount, unit_size, patterns);
//   return;
// }

template <typename Index>
void TestCount (Index &index, std::filesystem::path const &text_path)
{
  auto parent_index_path {CreateParentDirectoryByCategory("index", text_path)};
  {
    auto index_path {CreatePath(parent_index_path, text_path.filename().string(), ".index")};
    // if (!std::filesystem::exists(index_path))
    {
      ConstructIndex(index, text_path);
      SerializeIndex(index, index_path);
    }
    // else
    // {
    //   LoadIndex(index, index_path);
    // }
  }
  sdsl::csa_wt
  <
    sdsl::wt_rlmn<>,
    std::numeric_limits<uint32_t>::max(),
    std::numeric_limits<uint32_t>::max()
  >
  rlfm;
  {
    auto rlfm_path {CreatePath(parent_index_path, text_path.filename().string(), ".rlfm")};
    // if (!std::filesystem::exists(rlfm_path))
    {
      sdsl::int_vector<8> text;
      sdsl::load_vector_from_file(text, text_path);
      std::cout << "construct rlfm of " << std::filesystem::canonical(text_path) << "\n";
      sdsl::construct_im(rlfm, text);
      std::fstream rlfm_file(rlfm_path, std::ios_base::out | std::ios_base::trunc);
      std::cout << "serialize rlfm to " << std::filesystem::canonical(rlfm_path) << "\n";
      sdsl::serialize(rlfm, rlfm_file);
    }
    // else
    // {
    //   std::ifstream rlfm_file {rlfm_path};
    //   std::cout << "load rlfm from " << std::filesystem::canonical(rlfm_path) << "\n";
    //   sdsl::load(rlfm, rlfm_file);
    // }
  }
  {
    sdsl::int_vector<8> text;
    sdsl::load_vector_from_file(text, text_path);
    uint64_t max_unit_size {std::min(std::size(text), (1UL << 20))};
    std::cout << "test pattern counting\namount\tsize\n";
    for (uint64_t unit_size {1}; unit_size <= max_unit_size; unit_size *= 2)
    {
      auto patterns {ConcatenatedPatterns{std::min(max_unit_size / unit_size, (1UL << 13)), unit_size}};
      GenerateConcatnatedPatterns(text_path, patterns);
      auto begin {std::begin(patterns.labels)};
      auto end {begin};
      for (uint64_t i {}; i != patterns.amount; ++i)
      {
        begin = end;
        end = std::next(begin, patterns.unit_size);
        auto rlfm_count {sdsl::count(rlfm, begin, end)};
        auto index_count {Count(index, begin, end)};
        if (rlfm_count != index_count)
        {
          std::string pattern(begin, end);
          throw std::runtime_error
          (
            "\033[31mfailed at pattern: " + pattern
            + ", real count: " + std::to_string(rlfm_count)
            + ", count: " + std::to_string(index_count)
            + "\033[0m");
        }
      }
      std::cout << std::to_string(patterns.amount) << "\t" << std::to_string(unit_size) << "\t\033[32msucceed\033[0m\n";
    }
  }
  return;
}
}
