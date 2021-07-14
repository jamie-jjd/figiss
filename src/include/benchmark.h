#pragma once

#include <chrono>

#include <sdsl/csa_wt.hpp>
#include <sdsl/suffix_array_algorithm.hpp>

#include "grammar_compressed_index.h"
#include "pattern.h"
#include "utility.h"

namespace project
{
template <typename SdslIndex>
void SetUpSdslIndex
(
  std::filesystem::path const &text_path,
  SdslIndex &sdsl_index
)
{
  auto parent_path {CreateParentDirectoryByCategory("sdsl_index", text_path)};
  auto path {CreatePath(parent_path, text_path.filename().string(), ".sdsl")};
  // if (!std::filesystem::exists(path))
  {
    sdsl::int_vector<8> text;
    sdsl::load_vector_from_file(text, text_path);
    std::cout << "construct sdsl index of " << std::filesystem::canonical(text_path) << "\n";
    sdsl::construct_im(sdsl_index, text);
    std::fstream index_file(path, std::ios_base::out | std::ios_base::trunc);
    std::cout << "serialize sdsl index to " << std::filesystem::canonical(path) << "\n";
    sdsl::serialize(sdsl_index, index_file);
  }
  // else
  // {
  //   std::ifstream index_file {path};
  //   std::cout << "load sdsl index from " << std::filesystem::canonical(path) << "\n";
  //   sdsl::load(sdsl_index, index_file);
  // }
  return;
}

template <typename Index>
void SetUpIndex (std::filesystem::path const &text_path, Index &index)
{
  auto parent_path {CreateParentDirectoryByCategory("index", text_path)};
  auto path {CreatePath(parent_path, text_path.filename().string(), ".index")};
  // if (!std::filesystem::exists(path))
  {
    index = Index(text_path);
    index.Serialize(path);
  }
  // else
  // {
  //   index.Load(path);
  // }
  return;
}

template <typename Index>
void TestCount (std::filesystem::path const &text_path, Index &index)
{
  sdsl::csa_wt
  <
    sdsl::wt_rlmn<>,
    std::numeric_limits<uint32_t>::max(),
    std::numeric_limits<uint32_t>::max()
  >
  rlfm;
  SetUpSdslIndex(text_path, rlfm);
  SetUpIndex(text_path, index);
  {
    uint64_t amount {1UL << 12};
    uint64_t max_unit_size {1UL << 10};
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
      << std::to_string(patterns.GetUnitSize()) << "\t\033[32msucceed\033[0m\n";
    }
  }
  return;
}

template <typename SdslIndex>
void MeasureSdslIndexCountingTime
(
  std::vector<std::filesystem::path> const &pattern_paths,
  SdslIndex const &sdsl_index,
  std::vector<std::pair<uint64_t, double>> &time
)
{
  for (auto const &path : pattern_paths)
  {
    Patterns patterns;
    patterns.Load(path);
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
  return;
}

template <typename Index>
void MeasureIndexCountingTime
(
  std::vector<std::filesystem::path> const &pattern_paths,
  Index const &index,
  std::vector<std::pair<uint64_t, double>> &time
)
{
  for (auto const &path : pattern_paths)
  {
    Patterns patterns;
    patterns.Load(path);
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
  for (uint64_t unit_size {1}; unit_size <= (1ULL << 16); unit_size <<= 1)
  {
    auto patterns {Patterns(text_path, 1ULL << 14, unit_size)};
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
    SetUpIndex(text_path, index);
    MeasureIndexCountingTime(pattern_paths, index);
    PrintTime("time", text_path, time);
  }
  {
    std::vector<std::pair<uint64_t, double>> time;
    SetUpSdslIndex(text_path, sdsl_index);
    MeasureSdslIndexCountingTime(pattern_paths, sdsl_index);
    PrintTime("sdsl_time", text_path, time);
  }
  return;
}
}
