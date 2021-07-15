#pragma once

#include <chrono>

#include <sdsl/csa_wt.hpp>
#include <sdsl/suffix_array_algorithm.hpp>

#include "grammar_compressed_index.h"
#include "pattern.h"
#include "utility.h"

namespace gciis
{
template <typename Index>
void TestCount (std::filesystem::path const &text_path, Index &index)
{
  sdsl::csa_wt<sdsl::wt_rlmn<>, 0xFFFF'FFFF, 0xFFFF'FFFF> rlfm;
  {
    sdsl::int_vector<8> text;
    sdsl::load_vector_from_file(text, text_path);
    sdsl::construct_im(rlfm, text);
  }
  {
    index = Index {text_path};
  }
  {
    uint64_t amount {1ULL << 10};
    uint64_t max_unit_size {1ULL << 14};
    std::cout << "test pattern counting\namount\tsize\n";
    for (uint64_t unit_size {1}; unit_size <= max_unit_size; unit_size <<= 1)
    {
      Patterns patterns;
      for (uint8_t is_mutated {}; is_mutated != 2; ++is_mutated)
      {
        patterns = Patterns(text_path, amount, unit_size, is_mutated);
        auto begin {std::begin(patterns)};
        auto end {begin};
        for (uint64_t i {}; i != patterns.GetAmount(); ++i)
        {
          begin = end;
          end = std::next(begin, patterns.GetUnitSize());
          auto rlfm_count {sdsl::count(rlfm, begin, end)};
          auto index_count {index.Count(begin, end)};
          if (rlfm_count != index_count)
          {
            std::string pattern {begin, end};
            throw std::runtime_error
            (
              "\033[31mfailed at pattern:" + pattern
              + ", real count:" + std::to_string(rlfm_count)
              + ", count:" + std::to_string(index_count)
              + "\033[0m");
          }
        }
      }
      std::cout
      << std::to_string(patterns.GetAmount()) << "\t"
      << std::to_string(patterns.GetUnitSize()) << "\t\033[32msuccess\033[0m\n";
    }
  }
  return;
}

template <typename Index>
void MeasureIndexCountingTime
(
  std::filesystem::path const &text_path,
  std::vector<std::filesystem::path> const &pattern_paths,
  Index &index,
  std::vector<std::pair<uint64_t, double>> &time
)
{
  auto index_path {std::filesystem::path("_.index")};
  {
    index = Index{text_path};
    index.Serialize(index_path);
  }
  for (auto const &path : pattern_paths)
  {
    Patterns patterns;
    patterns.Load(path);
    index.Load(index_path);
    {
      auto begin {std::begin(patterns)};
      while (begin != std::end(patterns))
      {
        auto end {std::next(begin, patterns.GetUnitSize())};
        index.Count(begin, end);
        begin = end;
      }
    }
    {
      auto begin {std::begin(patterns)};
      auto begin_time {std::chrono::steady_clock::now()};
      while (begin != std::end(patterns))
      {
        auto end {std::next(begin, patterns.GetUnitSize())};
        index.Count(begin, end);
        begin = end;
      }
      time.emplace_back
      (
        patterns.GetUnitSize(),
        (
          static_cast<double>((std::chrono::steady_clock::now() - begin_time).count())
          / patterns.GetAmount()
          / patterns.GetUnitSize()
        )
      );
    }
  }
  std::filesystem::remove(index_path);
  return;
}

template <typename SdslIndex>
void MeasureSdslIndexCountingTime
(
  std::filesystem::path const &text_path,
  std::vector<std::filesystem::path> const &pattern_paths,
  SdslIndex &sdsl_index,
  std::vector<std::pair<uint64_t, double>> &time
)
{
  auto index_path {std::filesystem::path("_.sdsl")};
  {
    sdsl::int_vector<8> text;
    sdsl::load_vector_from_file(text, text_path);
    sdsl::construct_im(sdsl_index, text);
    std::ofstream fout {index_path, std::ios_base::out | std::ios_base::trunc};
    sdsl::serialize(sdsl_index, fout);
  }
  for (auto const &path : pattern_paths)
  {
    Patterns patterns;
    patterns.Load(path);
    {
      std::ifstream in {index_path};
      sdsl_index.load(in);
    }
    {
      auto begin {std::begin(patterns)};
      while (begin != std::end(patterns))
      {
        auto end {std::next(begin, patterns.GetUnitSize())};
        sdsl::count(sdsl_index, begin, end);
        begin = end;
      }
    }
    {
      auto begin {std::begin(patterns)};
      auto begin_time {std::chrono::steady_clock::now()};
      while (begin != std::end(patterns))
      {
        auto end {std::next(begin, patterns.GetUnitSize())};
        sdsl::count(sdsl_index, begin, end);
        begin = end;
      }
      time.emplace_back
      (
        patterns.GetUnitSize(),
        (
          static_cast<double>((std::chrono::steady_clock::now() - begin_time).count())
          / patterns.GetAmount()
          / patterns.GetUnitSize()
        )
      );
    }
  }
  std::filesystem::remove(index_path);
  return;
}

template <typename Index, typename SdslIndex>
void MeasureCountingTime
(
  std::filesystem::path const &text_path,
  Index &index,
  SdslIndex &sdsl_index
)
{
  std::vector<std::filesystem::path> pattern_paths;
  uint64_t amount {1ULL << 14};
  for (uint64_t unit_size {1}; unit_size <= (1ULL << 14); unit_size <<= 1)
  {
    auto patterns {Patterns(text_path, amount, unit_size)};
    pattern_paths.emplace_back
    (
      std::filesystem::path
      (
        text_path.filename().string() + "_" +
        std::to_string(patterns.GetAmount()) + "_" +
        std::to_string(patterns.GetUnitSize())
      )
    );
    patterns.Serialize(pattern_paths.back());
  }
  {
    std::vector<std::pair<uint64_t, double>> time;
    MeasureIndexCountingTime(text_path, pattern_paths, index, time);
    PrintCountingTime
    (
      std::filesystem::path
      {
        std::string("../data/time/")
        + std::to_string(Index::kMaxFactorSize) + "/"
        + text_path.filename().string()
      },
      time
    );
  }
  {
    std::vector<std::pair<uint64_t, double>> time;
    MeasureSdslIndexCountingTime(text_path, pattern_paths, sdsl_index, time);
    PrintCountingTime
    (
      std::filesystem::path
      {
        std::string("../data/time/")
        + "rlfm/"
        + text_path.filename().string()
      },
      time
    );
  }
  for (auto const path : pattern_paths)
  {
    std::filesystem::remove(path);
  }
  return;
}
}
