#ifndef GENERATE_PATTERN_HPP_
#define GENERATE_PATTERN_HPP_

#include <fstream>
#include <functional>
#include <random>

void generate_patterns
(
  char *const input_text_file,
  char *const output_patterns_file,
  uint32_t const number_patterns,
  uint32_t const size_pattern
)
{
  sdsl::int_vector<8> text;
  load_vector_from_file(text, input_text_file);
  if (size_pattern <= std::size(text))
  {
    std::fstream fout (output_patterns_file, std::ios_base::out | std::ios_base::binary);
    std::mt19937 engine {std::random_device{}()};
    std::uniform_int_distribution<uint32_t> distribution(0, std::size(text) - size_pattern);
    auto random_begin_dist {std::bind(distribution, engine)};
    sdsl::int_vector<8> pattern(size_pattern);
    for (uint32_t i {0}; i != number_patterns; ++i)
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
      pattern.serialize(fout);
    }
  }
  return;
}

#endif
