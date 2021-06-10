#pragma once

#include <chrono>

#include <sdsl/csa_wt.hpp>
#include <sdsl/suffix_array_algorithm.hpp>

#include "grammar_compressed_index.h"
#include "utility.h"

namespace project
{
template <typename SdslIndex>
void SetUpSdslIndex
(
  std::filesystem::path const &text_path,
  std::string const &index_name,
  SdslIndex &sdsl_index
)
{
  auto parent_path {CreateParentDirectoryByCategory("sdsl_index", text_path)};
  auto path {CreatePath(parent_path, text_path.filename().string(), "." + index_name)};
  if (!std::filesystem::exists(path))
  {
    sdsl::int_vector<8> text;
    sdsl::load_vector_from_file(text, text_path);
    std::cout << "construct " << index_name << " of " << std::filesystem::canonical(text_path) << "\n";
    sdsl::construct_im(sdsl_index, text);
    std::fstream index_file(path, std::ios_base::out | std::ios_base::trunc);
    std::cout << "serialize " << index_name << " to " << std::filesystem::canonical(path) << "\n";
    sdsl::serialize(sdsl_index, index_file);
  }
  else
  {
    std::ifstream index_file {path};
    std::cout << "load " << index_name << " from " << std::filesystem::canonical(path) << "\n";
    sdsl::load(sdsl_index, index_file);
  }
  return;
}

template <typename Index>
void SetUpIndex (std::filesystem::path const &text_path, Index &index)
{
  auto parent_path {CreateParentDirectoryByCategory("index", text_path)};
  {
    auto path {CreatePath(parent_path, text_path.filename().string(), ".index")};
    if (!std::filesystem::exists(path))
    {
      ConstructIndex(index, text_path);
      SerializeIndex(index, path);
    }
    else
    {
      LoadIndex(index, path);
    }
  }
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
  SetUpSdslIndex(text_path, "rlfm", rlfm);
  SetUpIndex(text_path, index);
  {
    sdsl::int_vector<8> text;
    sdsl::load_vector_from_file(text, text_path);
    uint64_t max_unit_size {std::min(std::size(text), (1UL << 20))};
    std::cout << "test pattern counting\namount\tsize\n";
    for (uint64_t unit_size {1}; unit_size <= max_unit_size; unit_size *= 2)
    {
      auto patterns {ConcatenatedPatterns{std::min(max_unit_size / unit_size, (1UL << 13)), unit_size}};
      GenerateConcatenatedPatterns(text_path, patterns);
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
      std::cout
      << std::to_string(patterns.amount) << "\t"
      << std::to_string(patterns.unit_size) << "\t\033[32msucceed\033[0m\n";
    }
  }
  return;
}

template <typename SdslIndex>
void MeasureSdslIndexCountingTime
(
  std::filesystem::path const &text_path,
  std::vector<std::filesystem::path> const &pattern_paths,
  SdslIndex &sdsl_index,
  std::string const &index_name
)
{
  std::vector<std::pair<uint64_t, double>> times;
  for (auto const &path : pattern_paths)
  {
    ConcatenatedPatterns patterns;
    LoadConcatenatedPatterns(path, patterns);
    SetUpSdslIndex(text_path, index_name, sdsl_index);
    auto begin {std::begin(patterns.labels)};
    auto begin_time {std::chrono::steady_clock::now()};
    while (begin != std::end(patterns.labels))
    {
      auto end {std::next(begin, patterns.unit_size)};
      sdsl::count(sdsl_index, begin, end);
      begin = end;
    }
    times.emplace_back
    (
      patterns.unit_size,
      (static_cast<double>((std::chrono::steady_clock::now() - begin_time).count())) / patterns.amount / patterns.unit_size
    );
  }
  PrintTime("sdsl_time", text_path, times);
  return;
}

template <typename Index>
void MeasureIndexCountingTime
(
  std::vector<std::filesystem::path> const &pattern_paths,
  std::filesystem::path const &text_path,
  Index &index
)
{
  std::vector<std::pair<uint64_t, double>> times;
  for (auto const &path : pattern_paths)
  {
    ConcatenatedPatterns patterns;
    LoadConcatenatedPatterns(path, patterns);
    SetUpIndex(text_path, index);
    auto begin {std::begin(patterns.labels)};
    auto begin_time {std::chrono::steady_clock::now()};
    while (begin != std::end(patterns.labels))
    {
      auto end {std::next(begin, patterns.unit_size)};
      Count(index, begin, end);
      begin = end;
    }
    times.emplace_back
    (
      patterns.unit_size,
      (static_cast<double>((std::chrono::steady_clock::now() - begin_time).count())) / patterns.amount / patterns.unit_size
    );
  }
  PrintTime("time", text_path, times);
  return;
}

void MeasureCoutingTime (std::filesystem::path const &text_path)
{
  std::vector<std::filesystem::path> pattern_paths;
  for (uint64_t pattern_size {1}; pattern_size <= (1ULL << 16); pattern_size <<= 1)
  {
    auto patterns {ConcatenatedPatterns{1ULL << 14, pattern_size}};
    pattern_paths.push_back(GenerateAndSerializeConcatenatedPatterns(text_path, patterns));
  }
  {
    Index index;
    MeasureIndexCountingTime(pattern_paths, text_path, index);
  }
  {
    sdsl::csa_wt
    <
      sdsl::wt_rlmn<>,
      std::numeric_limits<int32_t>::max(),
      std::numeric_limits<int32_t>::max()
    >
    rlfm;
    MeasureSdslIndexCountingTime(text_path, pattern_paths, rlfm, "rlfm");
  }
  return;
}
}
