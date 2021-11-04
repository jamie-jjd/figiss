#pragma once

#include <faster-minuter/wt_fbb.hpp>
#include <sdsl/csa_wt.hpp>
#include <sdsl/suffix_array_algorithm.hpp>

#include "index.h"
#include "pattern_collection.h"

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
  sdsl::csa_wt<wt_fbb<>, 0xFFFF'FFFF, 0xFFFF'FFFF> faster_minuter;
  {
    sdsl::int_vector<8> byte_text;
    sdsl::load_vector_from_file(byte_text, byte_text_path);
    sdsl::construct_im(faster_minuter, byte_text);
  }
  Index index;
  {
    index.Load(index_path);
  }
  {
    std::cout << "test counting on " << amount
    << " random (mutated) sample patterns of length between "
    << min_length << " and " << max_length << " from text\n";
    for (uint64_t length {min_length}; length <= max_length; length <<= 1)
    {
      for (uint8_t is_mutated {}; is_mutated != 2; ++is_mutated)
      {
        PatternCollection patterns;
        auto pattern_path
        {
          pattern_parent_path /
          byte_text_path.filename() /
          std::filesystem::path
          {
            std::to_string(amount) + "/" +
            std::to_string(length) + "/" +
            ((is_mutated) ? "mutated/" : "") +
            "pattern"
          }
        };
        if (!std::filesystem::exists(pattern_path))
        {
          std::filesystem::create_directories(pattern_path.parent_path());
          patterns = PatternCollection(byte_text_path, pattern_parent_path, amount, length, is_mutated);
          patterns.Serialize(pattern_path);
        }
        else
        {
          patterns.Load(pattern_path);
        }
        if (patterns)
        {
          {
            std::cout << "test " << amount << " random "
            << ((is_mutated) ? "mutated" : "")
            << "existed patterns of length " << length << "\n";
          }
          auto begin {std::begin(patterns)};
          auto end {begin};
          for (uint64_t i {}; i != patterns.GetAmount(); ++i)
          {
            begin = end;
            end = std::next(begin, patterns.GetLength());
            auto count {sdsl::count(faster_minuter, begin, end)};
            auto index_count {index.Count(begin, end)};
            if (count != index_count)
            {
              std::string pattern {begin, end};
              throw std::runtime_error
              (
                "\033[31mfailed at pattern:" + pattern
                + ",genuine count:" + std::to_string(count)
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
}
