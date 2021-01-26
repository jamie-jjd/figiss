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
    std::mt19937 engine {std::random_device{}()};
    std::uniform_int_distribution<uint64_t> distribution(0, std::size(text) - pattern_size);
    auto random_begin_dist {std::bind(distribution, engine)};
    sdsl::int_vector<8> pattern(pattern_size);
    pattern_output.write((char*)(&pattern_number), sizeof(pattern_number));
    for (uint64_t i {0}; i != pattern_number; ++i)
    {
      // auto text_it {std::next(std::begin(text), random_begin_dist())};
      auto text_it {std::next(std::begin(text), (i * pattern_size) % (std::size(text) - pattern_size))};
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

struct timer
{
  using inner_clock = std::chrono::high_resolution_clock;
  using inner_time_point = std::chrono::time_point<inner_clock>;

  std::unordered_map<std::string, std::pair<std::chrono::microseconds, inner_time_point>> time_records;

  void reset_timer (std::string const category)
  {
    time_records[category] = {std::chrono::microseconds{}, inner_time_point{}};
    return;
  }

  void resume_timer (std::string const category)
  {
    if (time_records.find(category) != time_records.end())
    {
      auto &begin_time {std::get<1>(time_records[category])};
      begin_time = inner_clock::now();
    }
    return;
  }

  void pause_timer (std::string const category)
  {
    if (time_records.find(category) != time_records.end())
    {
      auto &duration {std::get<0>(time_records[category])};
      auto &begin_time {std::get<1>(time_records[category])};
      duration = std::chrono::duration_cast<std::chrono::microseconds>(inner_clock::now() - begin_time);
    }
    return;
  }

  void print (std::string const time_path)
  {
    std::ofstream time_output {time_path};
    time_output << "category,milliseconds\n";
    for (auto const record : time_records)
    {
      auto &category {std::get<0>(record)};
      auto &microseconds {std::get<0>(std::get<1>(record))};
      time_output
      << category << ","
      << std::chrono::duration_cast<std::chrono::milliseconds>(microseconds).count() << "\n";
    }
    return;
  }
};
}
}

#endif
