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
    uint64_t amount {1ULL << 12};
    uint64_t max_unit_size {1ULL << 16};
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

void PrintCountingTime
(
  std::filesystem::path const &output_path,
  std::vector<std::pair<uint64_t, double>> time,
  bool const is_proper = true
)
{
  if (!std::filesystem::exists(output_path.parent_path()))
  {
    std::filesystem::create_directories(output_path.parent_path());
  }
  std::fstream fout {output_path, std::ios_base::out | std::ios_base::trunc};
  std::cout << "write time information to " << std::filesystem::canonical(output_path) << "\n";
  for (auto const &pair : time)
  {
    if (is_proper)
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

template <typename Index>
void MeasureIndexCountingTime
(
  std::filesystem::path const &text_path,
  Index &index
)
{
  auto index_path {std::filesystem::path("_.index")};
  {
    index = Index{text_path};
    index.Serialize(index_path);
  }
  std::vector<std::pair<uint64_t, double>> time;
  uint64_t amount {1ULL << 12};
  for (uint64_t unit_size {1ULL << 8}; unit_size <= (1ULL << 16); unit_size <<= 1)
  {
    std::cout << "pattern size: " << unit_size << "\n";
    auto patterns {Patterns(text_path, amount, unit_size)};
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
      auto duration {static_cast<double>((std::chrono::steady_clock::now() - begin_time).count())};
      time.emplace_back(patterns.GetUnitSize(), duration / patterns.GetAmount());
    }
  }
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
  std::filesystem::remove(index_path);
  return;
}

void MeasureRlfmCountingTime (std::filesystem::path const &text_path)
{
  auto index_path {std::filesystem::path("_.rlfm")};
  sdsl::csa_wt<sdsl::wt_rlmn<>, 0xFFFF'FFFF, 0xFFFF'FFFF> rlfm;
  {
    sdsl::int_vector<8> text;
    sdsl::load_vector_from_file(text, text_path);
    std::cout << "construct rlfm of " << std::filesystem::canonical(text_path) << "\n";
    sdsl::construct_im(rlfm, text);
    std::ofstream fout {index_path, std::ios_base::out | std::ios_base::trunc};
    sdsl::serialize(rlfm, fout);
    std::cout << "serialize rlfm to " << std::filesystem::canonical(index_path) << "\n";
  }
  std::vector<std::pair<uint64_t, double>> time;
  uint64_t amount {1ULL << 12};
  for (uint64_t unit_size {1ULL << 8}; unit_size <= (1ULL << 16); unit_size <<= 1)
  {
    std::cout << "pattern size: " << unit_size << "\n";
    auto patterns {Patterns(text_path, amount, unit_size)};
    {
      std::ifstream in {index_path};
      rlfm.load(in);
      std::cout << "load rlfm from " << std::filesystem::canonical(index_path) << "\n";
    }
    {
      auto begin {std::begin(patterns)};
      while (begin != std::end(patterns))
      {
        auto end {std::next(begin, patterns.GetUnitSize())};
        sdsl::count(rlfm, begin, end);
        begin = end;
      }
    }
    {
      auto begin {std::begin(patterns)};
      auto begin_time {std::chrono::steady_clock::now()};
      while (begin != std::end(patterns))
      {
        auto end {std::next(begin, patterns.GetUnitSize())};
        sdsl::count(rlfm, begin, end);
        begin = end;
      }
      auto duration {static_cast<double>((std::chrono::steady_clock::now() - begin_time).count())};
      time.emplace_back(patterns.GetUnitSize(), duration / patterns.GetAmount());
    }
  }
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
  std::filesystem::remove(index_path);
  return;
}
}
