#ifndef PARSE_TEXT_HPP_
#define PARSE_TEXT_HPP_

#include <algorithm>
#include <vector>
#include <sdsl/int_vector.hpp>
#include "character_bucket_vector.hpp"
#include "sl_type_vector.hpp"

template <typename T>
void f (T)
{
  std::cout << __PRETTY_FUNCTION__ << '\n';
  return;
}

template
<
  typename Text_Type,
  typename SL_Type_Vector,
  typename Text_Position_Vector_Type,
  typename Unique_Grammar_Rule_Vector_Type,
  typename Grammar_Compressed_Text_Type
>
void calculate_unique_grammar_rule_vector_and_grammar_compressed_text
(
  Text_Type const                 &text,
  SL_Type_Vector const            &SL_type_vector,
  Text_Position_Vector_Type const &text_position_vector,
  Unique_Grammar_Rule_Vector_Type &unique_grammar_rule_vector,
  Grammar_Compressed_Text_Type    &grammar_compressed_text
)
{
  f(text);
  f(SL_type_vector);
  f(text_position_vector);
  f(unique_grammar_rule_vector);
  f(grammar_compressed_text);
  return;
}

template
<
  typename Text_Type,
  typename Text_Position_Vector_Type
>
void collect_grammar_rules
(
  Text_Type const                &text,
  Text_Position_Vector_Type &text_position_vector
)
{
  text_position_vector.resize
  (
    std::distance
    (
      text_position_vector.begin(),
      std::stable_partition
      (
        text_position_vector.begin(),
        text_position_vector.end(),
        [&] (auto const &position)
        {
          return (position < text.size());
        }
      )
    )
  );
  return;
}

template
<
  typename Text_Type,
  typename SL_Type_Vector_Type,
  typename Character_Bucket_Vector_Type,
  typename Text_Position_Vector_Type
>
void induce_sort_S_type_substrings
(
  Text_Type const               &text,
  SL_Type_Vector_Type const     &SL_type_vector,
  Character_Bucket_Vector_Type  &character_bucket_vector,
  Text_Position_Vector_Type     &text_position_vector
)
{
  using Size_Type = typename Text_Position_Vector_Type::size_type;
  character_bucket_vector.calculate_cumulative_bucket_end(text);
  for (Size_Type rank {text_position_vector.size() - 1}; rank > 0; --rank)
  {
    auto position {text_position_vector[rank]};
    if
    (
      (position < text.size())
      &&
      (position > 0)
      &&
      SL_type_vector.is_S_type(position - 1)
    )
    {
      auto character {text[position - 1]};
      auto current_bucket_end {character_bucket_vector.get_current_bucket_end(character)};
      text_position_vector[current_bucket_end] = (position - 1);
      text_position_vector[rank] = text.size();
    }
  }
  return;
}

template
<
  typename Text_Type,
  typename SL_Type_Vector_Type,
  typename Character_Bucket_Vector_Type,
  typename Text_Position_Vector_Type
>
void induce_sort_L_type_substrings
(
  Text_Type const               &text,
  SL_Type_Vector_Type const     &SL_type_vector,
  Character_Bucket_Vector_Type  &character_bucket_vector,
  Text_Position_Vector_Type     &text_position_vector
)
{
  using Size_Type = typename Text_Position_Vector_Type::size_type;
  for (Size_Type rank {0}; rank < text_position_vector.size(); ++rank)
  {
    auto position {text_position_vector[rank]};
    if
    (
      (position < text.size())
      &&
      (position > 0)
      &&
      SL_type_vector.is_L_type(position - 1)
    )
    {
      auto character {text[position - 1]};
      auto current_bucket_begin {character_bucket_vector.get_current_bucket_begin(character)};
      text_position_vector[current_bucket_begin] = (position - 1);
      text_position_vector[rank] = text.size();
    }
  }
  return;
}

template
<
  typename Text_Type,
  typename SL_Type_Vector_Type,
  typename Character_Bucket_Vector_Type,
  typename Text_Position_Vector_Type
>
void bucket_sort_rightmost_L_type_characters
(
  Text_Type const               &text,
  SL_Type_Vector_Type const     &SL_type_vector,
  Character_Bucket_Vector_Type  &character_bucket_vector,
  Text_Position_Vector_Type     &text_position_vector
)
{
  using Size_Type = typename Text_Type::size_type;
  character_bucket_vector.calculate_cumulative_bucket_begin(text);
  for (Size_Type position {0}; position < (text.size() - 1); ++position)
  {
    if (SL_type_vector.is_rightmost_L_type(position))
    {
      auto character {text[position]};
      auto current_bucket_begin {character_bucket_vector.get_current_bucket_begin(character)};
      text_position_vector[current_bucket_begin] = position;
    }
  }
  return;
}

template
<
  typename Text_Type,
  typename Unique_Grammar_Rule_Vector_Type,
  typename Grammar_Compressed_Text_Type
>
void parse_text
(
  Text_Type const &text,
  Unique_Grammar_Rule_Vector_Type &unique_grammar_rule_vector,
  Grammar_Compressed_Text_Type &grammar_compressed_text
)
{
  SL_Type_Vector<> SL_type_vector {text};
  Character_Bucket_Vector<> character_bucket_vector {text};
  sdsl::int_vector<> text_position_vector(text.size(), text.size());

  bucket_sort_rightmost_L_type_characters
  (
    text,
    SL_type_vector,
    character_bucket_vector,
    text_position_vector
  );

  induce_sort_L_type_substrings
  (
    text,
    SL_type_vector,
    character_bucket_vector,
    text_position_vector
  );

  induce_sort_S_type_substrings
  (
    text,
    SL_type_vector,
    character_bucket_vector,
    text_position_vector
  );

  collect_grammar_rules
  (
    text,
    text_position_vector
  );

  calculate_unique_grammar_rule_vector_and_grammar_compressed_text
  (
    text,
    SL_type_vector,
    text_position_vector,
    unique_grammar_rule_vector,
    grammar_compressed_text
  );

  return;
}

#endif
