#ifndef GENERATE_PATTERN_HPP_
#define GENERATE_PATTERN_HPP_

#include <fstream>
#include <functional>
#include <random>

namespace gci
{
void generate_patterns
(
  std::string const text_path,
  std::string const patterns_path,
  uint32_t const number_patterns,
  uint32_t const pattern_size
)
{
  sdsl::int_vector<8> text;
  sdsl::load_vector_from_file(text, text_path);
  if (pattern_size > std::size(text))
  {
    throw std::runtime_error("pattern size is larger than text size");
  }
  std::fstream fout (patterns_path, std::ios_base::out | std::ios_base::binary);
  std::mt19937 engine {std::random_device{}()};
  std::uniform_int_distribution<uint32_t> distribution(0, std::size(text) - pattern_size);
  auto random_begin_dist {std::bind(distribution, engine)};
  sdsl::int_vector<8> pattern(pattern_size);
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
  return;
}
}

#endif