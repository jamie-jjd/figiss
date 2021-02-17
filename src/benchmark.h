#pragma once

#include "grammar_compressed_index.h"
#include "utility.h"

namespace project
{
template
<
  typename Index,
  typename Fmindex
>
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
  Load(index, index_path);
  auto fm_index_path {parent_index_path / (text_path.filename().string() + ".fmindex")};
  if (!std::filesystem::exists(fm_index_path))
  {
    sdsl::int_vector<> text;
    std::ifstream text_file {text_path};
    text.load(text_file);
    sdsl::construct_im(fm_index, text);
    std::ofstream fm_index_file {fm_index_path};
    sdsl::serialize(fm_index, fm_index_file);
  }
  {
    std::ifstream fm_index_file {fm_index_path};
    sdsl::load(fm_index, fm_index_file);
  }
  return;
}

void BenchmarkIndexCount
(
  Index &index,
  std::filesystem::path const &pattern_path
)
{
  std::ifstream pattern_file {pattern_path};
  sdsl::int_vector<> pattern;
  uint64_t pattern_amount {0};
  pattern_file >> pattern_amount;
  for (uint64_t i {0}; i != pattern_amount; ++i)
  {
    pattern.load(pattern_file);
    Count(index, std::begin(pattern), std::end(pattern));
  }
  return;
}

template <typename Fmindex>
void BenchmarkFmindexCount
(
  Fmindex &fm_index,
  std::filesystem::path const &pattern_path
)
{
  std::ifstream pattern_file {pattern_path};
  sdsl::int_vector<> pattern;
  uint64_t pattern_amount {0};
  pattern_file >> pattern_amount;
  for (uint64_t i {0}; i != pattern_amount; ++i)
  {
    pattern.load(pattern_file);
    sdsl::count(fm_index, std::begin(pattern), std::end(pattern));
  }
  return;
}

template
<
  typename Index,
  typename Fmindex
>
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
  auto pattern_path {GeneratePattern(pattern_path, pattern_amount, pattern_size)};
  Index index;
  Fmindex fm_index;
  LoadIndexAndFmindex(index, fm_index, text_path);
  BenchmarkIndexCount(index, pattern_path);
  BenchmarkFmindexCount(fm_index, pattern_path);
  sdsl::remove(pattern_path);
  return;
}

template
<
  typename Index,
  typename Fmindex
>
void TestCount (std::filesystem::path const &text_path)
{
  Index index;
  Fmindex fm_index;
  LoadIndexAndFmindex(index, fm_index, text_path);
  sdsl::int_vector<> text;
  std::ifstream text_file {text_path};
  text.load(text_file);
  auto max_sl_factor_size {CalculateMaxSlFactorSize(text)};
  for (uint64_t multiple {2}; multiple != 11; ++multiple)
  {
    uint64_t pattern_amount {1000 / multiple};
    auto pattern_size {static_cast<uint64_t>(max_sl_factor_size * multiple * 0.5)};
    auto pattern_path {GeneratePattern(text_path, pattern_amount, pattern_size)};
    std::ifstream pattern_file {pattern_path};
    sdsl::int_vector<> pattern_amount_size;
    pattern_amount_size.load(pattern_file);
    pattern_amount = pattern_amount_size[0];
    pattern_size = pattern_amount_size[1];
    sdsl::int_vector<> patterns;
    patterns.load(pattern_file);
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
    std::filesystem::remove(pattern_path);
  }
  for (uint64_t divisor {5}; divisor != 1; --divisor)
  {
    uint64_t pattern_size {std::size(text) / divisor};
    if (pattern_size >= max_sl_factor_size)
    {
      uint64_t pattern_amount {1};
      auto pattern_path {GeneratePattern(text_path, pattern_amount, pattern_size)};
      std::ifstream pattern_file {pattern_path};
      sdsl::int_vector<> pattern_amount_size;
      pattern_amount_size.load(pattern_file);
      pattern_amount = pattern_amount_size[0];
      pattern_size = pattern_amount_size[1];
      sdsl::int_vector<> patterns;
      patterns.load(pattern_file);
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
      std::filesystem::remove(pattern_path);
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
