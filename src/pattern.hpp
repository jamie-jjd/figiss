#ifndef GENERATE_PATTERN_HPP_
#define GENERATE_PATTERN_HPP_

#include <fstream>
#include <functional>
#include <random>

namespace gci
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
      auto text_it {std::next(std::begin(text), random_begin_dist())};
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
}

#endif
