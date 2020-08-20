#ifndef UTILITY_HPP_
#define UTILITY_HPP_

#include <fstream>
#include <vector>

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
