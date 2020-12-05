#ifndef UTILITY_HPP_
#define UTILITY_HPP_

#include <algorithm>
#include <fstream>
#include <numeric>
#include <vector>

template <typename Type>
class Type_Displayer;

template <typename Type>
void display_type (Type &&parameter)
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

template <typename Character_Type>
std::ofstream& operator<< (std::ofstream &output_file_stream, std::vector<Character_Type> const &vector_)
{
  if (!vector_.empty())
  {
    output_file_stream << *(vector_.cbegin());
    for (auto vector_iterator {++(vector_.cbegin())}; vector_iterator != vector_.cend(); ++vector_iterator)
    {
      output_file_stream << ',' << *vector_iterator;
    }
    output_file_stream << '\n';
  }
  return output_file_stream;
}

#endif
