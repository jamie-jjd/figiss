#ifndef PARSE_TEXT_HPP_
#define PARSE_TEXT_HPP_

#include <vector>
#include <sdsl/int_vector.hpp>
#include "sl_type_vector.hpp"

void bucket_sort_rightmost_L_type_characters
(
  sdsl::int_vector<8> const &text,
  SL_Type_Vector const &SL_type_vector,
  sdsl::int_vector<32> &text_position_vector
)
{
  return;
}

void induce_sort_L_type_substrings
(
  sdsl::int_vector<8> const &text,
  SL_Type_Vector const &SL_type_vector,
  sdsl::int_vector<32> &text_position_vector
)
{

}

void induce_sort_S_type_substrings
(
  sdsl::int_vector<8> const &text,
  SL_Type_Vector const &SL_type_vector,
  sdsl::int_vector<32> &text_position_vector
)
{
  return;
}

void induce_sort_terminal_substrings
(
  sdsl::int_vector<8> const &text,
  SL_Type_Vector const &SL_type_vector,
  sdsl::int_vector<32> &text_position_vector
)
{
  bucket_sort_rightmost_L_type_characters(text, SL_type_vector, text_position_vector);
  induce_sort_L_type_substrings(text, SL_type_vector, text_position_vector);
  induce_sort_S_type_substrings(text, SL_type_vector, text_position_vector);
  return;
}

void parse_text
(
  sdsl::int_vector<8> const &text,
  std::vector<sdsl::int_vector<8>> &terminal_substring_vector,
  sdsl::int_vector<32> &non_terminal_text
)
{
  auto SL_type_vector {SL_Type_Vector<sdsl::bit_vector, sdsl::int_vector<8>>(text)};
  auto text_position_vector {sdsl::int_vector<32>(text.size(), text.size())};

  induce_sort_terminal_substrings(text, SL_type_vector, text_position_vector);

  // TODO: terminal_substring_vector and non_terminal_text
  return;
}

#endif
