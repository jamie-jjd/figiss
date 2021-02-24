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
  auto parent_index_path
  {
    std::filesystem::path{"../data/index"}
    / text_path.parent_path().filename()
  };
  if (!std::filesystem::exists(parent_index_path))
  {
    std::filesystem::create_directories(parent_index_path);
  }
  auto index_path {parent_index_path / (text_path.filename().string() + ".index")};
  if (!std::filesystem::exists(index_path))
  {
    Construct(index, text_path);
    Serialize(index, index_path);
  }
  else
  {
    Load(index, index_path);
  }
  auto fm_index_path {parent_index_path / (text_path.filename().string() + ".fmindex")};
  if (!std::filesystem::exists(fm_index_path))
  {
    sdsl::int_vector<8> text;
    sdsl::load_vector_from_file(text, text_path);
    sdsl::construct_im(fm_index, text);
    std::ofstream fm_index_file {fm_index_path};
    sdsl::serialize(fm_index, fm_index_file);
  }
  else
  {
    std::ifstream fm_index_file {fm_index_path};
    sdsl::load(fm_index, fm_index_file);
  }
  return;
}

template
<
  typename PatternAmount,
  typename PatternSize,
  typename Patterns
>
void BenchmarkIndexCount
(
  Index const &index,
  PatternAmount const pattern_amount,
  PatternSize const pattern_size,
  Patterns &patterns
)
{
  auto pattern_begin {std::begin(patterns)};
  auto pattern_end {pattern_begin};
  for (uint64_t i {0}; i != pattern_amount; ++i)
  {
    pattern_begin = pattern_end;
    pattern_end = std::next(pattern_begin, pattern_size);
    Count(index, pattern_begin, pattern_end);
  }
  return;
}

template
<
  typename Fmindex,
  typename PatternAmount,
  typename PatternSize,
  typename Patterns
>
void BenchmarkFmindexCount
(
  Fmindex const &fm_index,
  PatternAmount const pattern_amount,
  PatternSize const pattern_size,
  Patterns &patterns
)
{
  auto pattern_begin {std::begin(patterns)};
  auto pattern_end {pattern_begin};
  for (uint64_t i {0}; i != pattern_amount; ++i)
  {
    pattern_begin = pattern_end;
    pattern_end = std::next(pattern_begin, pattern_size);
    sdsl::count(fm_index, pattern_begin, pattern_end);
  }
  return;
}

void BenchmarkCount
(
  std::filesystem::path const &text_path,
  uint64_t const pattern_amount,
  uint64_t const pattern_size
)
{
  auto min_pattern_size {CalculateMaxSlFactorSize(text_path)};
  if (pattern_size < min_pattern_size)
  {
    throw std::runtime_error
    (
      std::string{"pattern size should be at least "}
      + std::to_string(min_pattern_size)
      + " (characters) for this text"
    );
  }
  Index index;
  sdsl::csa_wt<> fm_index;
  LoadIndexAndFmindex(index, fm_index, text_path);
  sdsl::int_vector<8> patterns;
  GeneratePatterns(text_path, pattern_amount, pattern_size, patterns);
  BenchmarkIndexCount(index, pattern_amount, pattern_size, patterns);
  BenchmarkFmindexCount(fm_index, pattern_amount, pattern_size, patterns);
  return;
}

void TestCount (std::filesystem::path const &text_path)
{
  Index index;
  sdsl::csa_wt<> fm_index;
  LoadIndexAndFmindex(index, fm_index, text_path);
  sdsl::int_vector<8> text;
  sdsl::load_vector_from_file(text, text_path);
  auto max_sl_factor_size {CalculateMaxSlFactorSize(text)};
  for (uint64_t multiple {2}; multiple != 11; ++multiple)
  {
    uint64_t pattern_amount {1000 / multiple};
    auto pattern_size {static_cast<uint64_t>(max_sl_factor_size * multiple * 0.5)};
    sdsl::int_vector<8> patterns;
    GeneratePatterns(text_path, pattern_amount, pattern_size, patterns);
    auto pattern_begin {std::begin(patterns)};
    auto pattern_end {pattern_begin};
    for (uint64_t i {0}; i != pattern_amount; ++i)
    {
      pattern_begin = pattern_end;
      pattern_end = std::next(pattern_begin, pattern_size);
      auto fm_index_count {sdsl::count(fm_index, pattern_begin, pattern_end)};
      auto index_count {Count(index, pattern_begin, pattern_end)};
      if (fm_index_count != index_count)
      {
        throw std::runtime_error
        (
          "failed at pattern size: " +
          std::to_string(pattern_size)
        );
      }
    }
    std::cout << "succeed at pattern size: " << std::to_string(pattern_size) << "\n";
  }
  for (uint64_t divisor {5}; divisor != 1; --divisor)
  {
    uint64_t pattern_size {std::size(text) / divisor};
    if (pattern_size >= max_sl_factor_size)
    {
      uint64_t pattern_amount {1};
      sdsl::int_vector<8> patterns;
      GeneratePatterns(text_path, pattern_amount, pattern_size, patterns);
      auto pattern_begin {std::begin(patterns)};
      auto pattern_end {pattern_begin};
      for (uint64_t i {0}; i != pattern_amount; ++i)
      {
        pattern_begin = pattern_end;
        pattern_end = std::next(pattern_begin, pattern_size);
        auto fm_index_count {sdsl::count(fm_index, pattern_begin, pattern_end)};
        auto index_count {Count(index, pattern_begin, pattern_end)};
        if (fm_index_count != index_count)
        {
          throw std::runtime_error
          (
            "failed at pattern size: " +
            std::to_string(pattern_size)
          );
        }
      }
      std::cout << "succeed at pattern size: " + std::to_string(pattern_size) << "\n";
    }
  }
  auto fm_index_count {sdsl::count(fm_index, std::begin(text), std::end(text))};
  auto index_count {Count(index, std::begin(text), std::end(text))};
  if (fm_index_count != index_count)
  {
    throw std::runtime_error("failed at text size");
  }
  std::cout << "succeed at text size\n";
  return;
}
}
