#ifndef UTILITY_HPP_
#define UTILITY_HPP_

#include <algorithm>
#include <fstream>
#include <numeric>
#include <vector>

template <typename Character_Type>
int32_t get_effective_alphabet_size (std::vector<Character_Type> const &text)
{
  auto is_effective_characters {std::vector<int32_t>(static_cast<int32_t>(*std::max_element(text.begin(), text.end())) + 1, 0)};
  for (const auto &character : text) { is_effective_characters[character] = 1; }
  return std::accumulate(is_effective_characters.begin(), is_effective_characters.end(), static_cast<int32_t>(0));
}

template <typename Value_Type>
void output_comma_separated_vector
(
  std::ofstream &output_file_stream,
  std::vector<Value_Type> const &vector_
)
{
  if (!vector_.empty())
  {
    auto vector_iterator {vector_.cbegin()};
    output_file_stream << (*vector_iterator);

    while ((++vector_iterator) != vector_.cend())
    {
      output_file_stream << ',' << (*vector_iterator);
    }
    output_file_stream << '\n';
  }
}

#endif
