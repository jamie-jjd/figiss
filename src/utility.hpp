#ifndef UTILITY_HPP_
#define UTILITY_HPP_

#include <chrono>
#include <fstream>
#include <functional>
#include <random>
#include <string>
#include <unordered_map>

namespace gci
{
namespace util
{
void generate_pattern
(
  std::string const text_path,
  std::string const pattern_path,
  uint64_t const pattern_number,
  uint64_t const pattern_size
)
{
  sdsl::int_vector<8> text;
  sdsl::load_vector_from_file(text, text_path);
  if (pattern_size < std::size(text))
  {
    std::ofstream pattern_output {pattern_path};
#ifndef LOG_VERBOSE_GCI_COUNT
    std::mt19937 engine {std::random_device{}()};
    std::uniform_int_distribution<uint64_t> distribution(0, std::size(text) - pattern_size);
    auto random_begin_dist {std::bind(distribution, engine)};
#endif
    sdsl::int_vector<8> pattern(pattern_size);
    pattern_output.write((char*)(&pattern_number), sizeof(pattern_number));
    for (uint64_t i {0}; i != pattern_number; ++i)
    {
#ifndef LOG_VERBOSE_GCI_COUNT
      auto text_it {std::next(std::begin(text), random_begin_dist())};
#else
      auto text_it {std::next(std::begin(text), (i * pattern_size) % (std::size(text) - pattern_size))};
#endif
      auto pattern_it {std::begin(pattern)};
      auto pattern_end {std::end(pattern)};
      while (pattern_it != pattern_end)
      {
        *pattern_it = *text_it;
        ++text_it;
        ++pattern_it;
      }
      pattern.serialize(pattern_output);
    }
  }
  else
  {
    std::cout << "pattern size >= text size\n";
    std::cout << "generated pattern is one copy of the text\n";
    std::ofstream pattern_output {pattern_path};
    uint64_t pattern_number_ {1};
    pattern_output.write((char*)(&pattern_number_), sizeof(pattern_number_));
    text.serialize(pattern_output);
  }
  return;
}

std::string basename (std::string const path)
{
  auto basename_begin_dist {path.find_last_of("/")};
  if (basename_begin_dist != std::string::npos)
  {
    ++basename_begin_dist;
  }
  else
  {
    basename_begin_dist = 0;
  }
  return path.substr(basename_begin_dist);
}

enum class category_id
{
#if defined(LOG_GCI_COUNT) || defined(LOG_VERBOSE_GCI_COUNT)
  GCI_COUNT,
#endif
#if defined(LOG_VERBOSE_GCI_COUNT)
  CALCULATE_SL_FACTOR,
  LOOKUP_GRAMMAR_RULE,
  WAVELET_TREE_OPERATIONS_L,
  WAVELET_TREE_OPERATIONS_S,
#endif
#ifdef LOG_FMI_COUNT
  FMI_COUNT,
#endif
  DUMMY
};

struct timer
{
  using inner_clock = std::chrono::high_resolution_clock;

  static std::unordered_map
  <
    category_id,
    std::tuple
    <
      std::string,
      std::chrono::nanoseconds,
      inner_clock::time_point
    >
  >
  records;

  static void register_category
  (
    category_id id,
    std::string const category
  )
  {
    records[id] =
    {
      category,
      std::chrono::nanoseconds::zero(),
      inner_clock::time_point::min()
    };
    return;
  }

  static void resume (category_id id)
  {
    std::get<2>(records[id]) = inner_clock::now();
    return;
  }

  static void pause (category_id id)
  {
    std::get<1>(records[id]) +=
    std::chrono::duration_cast<std::chrono::nanoseconds>
    (
      inner_clock::now()
      -
      std::get<2>(records[id])
    );
    return;
  }

  static void print
  (
    std::string const path,
    std::string const unit = "seconds"
  )
  {
    std::ofstream fout {path};
    fout << "category," << unit << "\n";
    for (auto const &record : records)
    {
      auto category {std::get<0>(std::get<1>(record))};
      auto duration {std::get<1>(std::get<1>(record))};
      fout << category << ",";
      if (unit == "minutes")
      {
        fout << std::chrono::duration_cast<std::chrono::minutes>(duration).count() << "\n";
      }
      else if (unit == "seconds")
      {
        fout << std::chrono::duration_cast<std::chrono::seconds>(duration).count() << "\n";
      }
      else if (unit == "milliseconds")
      {
        fout << std::chrono::duration_cast<std::chrono::milliseconds>(duration).count() << "\n";
      }
      else if (unit == "microseconds")
      {
        fout << std::chrono::duration_cast<std::chrono::microseconds>(duration).count() << "\n";
      }
      else
      {
        fout << duration.count() << "\n";
      }
    }
    records.clear();
    return;
  }
};

std::unordered_map
<
  category_id,
  std::tuple
  <
    std::string,
    std::chrono::nanoseconds,
    timer::inner_clock::time_point
  >
>
timer::records
{};

}
}

#endif
