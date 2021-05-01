#pragma once

#include "grammar_compressed_index.h"
#include "utility.h"

namespace project
{
template <typename Fmindex>
void LoadIndexAndFmindex
(
  Index &index,
  Fmindex &fm_index,
  std::filesystem::path const &text_path
)
{
  auto parent_index_path {CreateParentDirectoryByCategory("index", text_path)};
  auto index_path {CreatePath(parent_index_path, text_path.filename().string(), ".index")};
  // if (!std::filesystem::exists(index_path))
  {
    std::cout << "construct & serialize index ...\n";
    ConstructIndex(index, text_path);
    SerializeIndex(index, index_path);
  }
  // else
  // {
  //   std::cout << "load index ...\n";
  //   LoadIndex(index, index_path);
  // }
  auto fm_index_path {CreatePath(parent_index_path, text_path.filename().string(), ".fm_index")};
  // if (!std::filesystem::exists(fm_index_path))
  {
    sdsl::int_vector<8> text;
    sdsl::load_vector_from_file(text, text_path);
    std::cout << "construct & serialize FM-index ...\n";
    sdsl::construct_im(fm_index, text);
    std::fstream fm_index_file(fm_index_path, std::ios_base::out | std::ios_base::trunc);
    sdsl::serialize(fm_index, fm_index_file);
  }
  // else
  // {
  //   std::ifstream fm_index_file {fm_index_path};
  //   std::cout << "load FM-index ...\n";
  //   sdsl::load(fm_index, fm_index_file);
  // }
  return;
}

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
//   Fmindex const &fm_index,
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
//     sdsl::count(fm_index, begin, end);
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
//   sdsl::csa_wt<> fm_index;
//   LoadIndexAndFmindex(index, fm_index, text_path);
//   sdsl::int_vector<8> patterns;
//   GeneratePatterns(text_path, pattern_amount, unit_size, patterns);
//   BenchmarkIndexCount(index, pattern_amount, unit_size, patterns);
//   BenchmarkFmindexCount(fm_index, pattern_amount, unit_size, patterns);
//   return;
// }

void TestCount (std::filesystem::path const &text_path)
{
  Index index;
  sdsl::csa_wt<> fm_index;
  LoadIndexAndFmindex(index, fm_index, text_path);
  sdsl::int_vector<8> text;
  sdsl::load_vector_from_file(text, text_path);
  for (uint64_t unit_size {1}; (unit_size <= 100) && (unit_size < std::size(text) + 1); ++unit_size)
  {
    auto patterns {ConcatenatedPatterns{100, unit_size}};
    GenerateConcatnatedPatterns(text_path, patterns);
    auto begin {std::begin(patterns.labels)};
    auto end {begin};
    for (uint64_t i {}; i != patterns.amount; ++i)
    {
      begin = end;
      end = std::next(begin, patterns.unit_size);
      auto fm_index_count {sdsl::count(fm_index, begin, end)};
      auto index_count {Count(index, begin, end)};
      if (fm_index_count != index_count)
      {
        throw std::runtime_error
        (
          "failed at size: " +
          std::to_string(unit_size)
        );
      }
    }
    std::cout << "succeed at size: " << std::to_string(unit_size) << "\n";
  }
  for (uint64_t unit_size {std::size(text)}; unit_size > 100; unit_size /= 2)
  {
    auto patterns {ConcatenatedPatterns{1, unit_size}};
    GenerateConcatnatedPatterns(text_path, patterns);
    auto begin {std::begin(patterns.labels)};
    auto end {std::end(patterns.labels)};
    auto fm_index_count {sdsl::count(fm_index, begin, end)};
    auto index_count {Count(index, begin, end)};
    if (fm_index_count != index_count)
    {
      throw std::runtime_error
      (
        "failed at size: " +
        std::to_string(unit_size)
      );
    }
    std::cout << "succeed at size: " + std::to_string(unit_size) << "\n";
  }
  return;
}
}
