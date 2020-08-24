#ifndef UTILITY_HPP_
#define UTILITY_HPP_

#include <algorithm>
#include <fstream>
#include <numeric>
#include <vector>

template <typename Type>
class Type_Displayer;

template <typename Type>
void display_type (Type &parameter)
{
  Type_Displayer<Type> type_;
  Type_Displayer<decltype(parameter)> paramter_type_;
}

template <typename Character_Type>
int32_t get_effective_alphabet_size (std::vector<Character_Type> const &text)
{
  auto is_effective_characters {std::vector<int32_t>(static_cast<int32_t>(*std::max_element(text.begin(), text.end())) + 1, 0)};
  for (const auto &character : text) { is_effective_characters[character] = 1; }
  return std::accumulate(is_effective_characters.begin(), is_effective_characters.end(), static_cast<int32_t>(0));
}

template <typename Const_Iterator_Type>
void output_comma_separated_vector
(
  std::ofstream &output_file_stream,
  Const_Iterator_Type vector_cbegin,
  Const_Iterator_Type vector_cend
)
{
  if (vector_cbegin != vector_cend)
  {
    output_file_stream << (*vector_cbegin);
    for (auto it {++vector_cbegin}; it != vector_cend; ++it)
    {
      output_file_stream << ',' << (*it);
    }
    output_file_stream << '\n';
  }
}

#endif
