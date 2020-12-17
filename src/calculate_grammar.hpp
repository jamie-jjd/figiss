#ifndef CALCULATE_GRAMMAR_HPP_
#define CALCULATE_GRAMMAR_HPP_

#include <sdsl/int_vector.hpp>

#include "sl_type.hpp"

template
<
  typename text_type,
  typename character_bucket_range_vector_type
>
void calculate_cumulative_character_bucket_end
(
  text_type const &text,
  character_bucket_range_vector_type &character_bucket_range_vector
)
{
  auto begin {std::begin(character_bucket_range_vector)};
  auto end {std::end(character_bucket_range_vector)};
  std::fill(begin, end, 0);
  for (auto const &character : text)
  {
    ++character_bucket_range_vector[character];
  }
  std::partial_sum(begin, end, begin);
  return;
}

template
<
  typename text_type,
  typename character_bucket_range_vector_type
>
void calculate_cumulative_character_bucket_begin
(
  text_type const &text,
  character_bucket_range_vector_type &character_bucket_range_vector
)
{
  calculate_cumulative_character_bucket_end
  (
    text,
    character_bucket_range_vector
  );
  auto rfirst {std::prev(std::end(character_bucket_range_vector))};
  auto rlast {std::begin(character_bucket_range_vector)};
  for (auto rit {rfirst}; rit != rlast; --rit)
  {
    *rit = *std::prev(rit);
  }
  *rlast = 0;
  return;
}

template
<
  typename text_type,
  typename sl_type_vector_type,
  typename character_bucket_range_vector_type,
  typename text_distance_vector_type
>
void bucket_sort_rightmost_l_type_characters
(
  text_type const &text,
  sl_type_vector_type const &sl_type_vector,
  character_bucket_range_vector_type &character_bucket_range_vector,
  text_distance_vector_type &text_distance_vector
)
{
  for (auto text_it {std::begin(text)}; text_it != std::end(text); ++text_it)
  {
    auto text_distance {std::distance(std::begin(text), text_it)};
    if (is_rightmost_l_type(std::next(std::begin(sl_type_vector), text_distance)))
    {
      auto bucket_begin_distance {character_bucket_range_vector[*text_it]++};
      text_distance_vector[bucket_begin_distance] = text_distance;
    }
  }
  return;
}

template
<
  typename text_type,
  typename sl_type_vector_type,
  typename character_bucket_range_vector_type,
  typename text_distance_vector_type
>
void induce_sort_l_type_characters
(
  text_type const &text,
  sl_type_vector_type const &sl_type_vector,
  character_bucket_range_vector_type &character_bucket_range_vector,
  text_distance_vector_type &text_distance_vector
)
{
  auto text_distance_begin {std::begin(text_distance_vector)};
  auto text_distance_end {std::end(text_distance_vector)};
  auto invalid_text_distance {std::size(text)};
  for
  (
    auto text_distance_it {std::next(text_distance_begin)};
    text_distance_it != text_distance_end;
    ++text_distance_it
  )
  {
    auto text_distance {*text_distance_it};
    if
    (
      (text_distance != invalid_text_distance)
      &&
      (text_distance != 0)
      &&
      (sl_type_vector[text_distance - 1] == L_TYPE)
    )
    {
      auto bucket_begin_distance {character_bucket_range_vector[text[text_distance - 1]]++};
      text_distance_vector[bucket_begin_distance] = (text_distance - 1);
      *text_distance_it = invalid_text_distance;
    }
  }
  return;
}

template
<
  typename text_type,
  typename sl_type_vector_type,
  typename character_bucket_range_vector_type,
  typename text_distance_vector_type
>
void induce_sort_s_type_characters
(
  text_type const &text,
  sl_type_vector_type const &sl_type_vector,
  character_bucket_range_vector_type &character_bucket_range_vector,
  text_distance_vector_type &text_distance_vector
)
{
  auto text_distance_rfirst {std::prev(std::end(text_distance_vector))};
  auto text_distance_rlast {std::begin(text_distance_vector)};
  auto invalid_text_distance {std::size(text)};
  for (auto text_distance_rit {text_distance_rfirst}; text_distance_rit != text_distance_rlast; --text_distance_rit)
  {
    auto text_distance {*text_distance_rit};
    if
    (
      (text_distance != invalid_text_distance)
      &&
      (text_distance != 0)
      &&
      (sl_type_vector[text_distance - 1] == S_TYPE)
    )
    {
      auto bucket_end_distance {--character_bucket_range_vector[text[text_distance - 1]]};
      text_distance_vector[bucket_end_distance] = (text_distance - 1);
      *text_distance_rit = invalid_text_distance;
    }
  }
  return;
}

template
<
  typename text_type,
  typename sl_type_vector_type,
  typename text_distance_vector_type
>
void induce_sort_grammar_rule
(
  text_type const &text,
  sl_type_vector_type &sl_type_vector,
  text_distance_vector_type &text_distance_vector
)
{
  auto invalid_text_distance {text.size()};
  text_distance_vector.resize(text.size());
  std::fill
  (
    std::begin(text_distance_vector),
    std::end(text_distance_vector),
    invalid_text_distance
  );
  auto character_upper_bound {*std::max_element(std::begin(text), std::end(text)) + 1};
  sdsl::int_vector<> character_bucket_range_vector(character_upper_bound);
  calculate_cumulative_character_bucket_begin
  (
    text,
    character_bucket_range_vector
  );
  bucket_sort_rightmost_l_type_characters
  (
    text,
    sl_type_vector,
    character_bucket_range_vector,
    text_distance_vector
  );
  induce_sort_l_type_characters
  (
    text,
    sl_type_vector,
    character_bucket_range_vector,
    text_distance_vector
  );
  calculate_cumulative_character_bucket_end
  (
    text,
    character_bucket_range_vector
  );
  induce_sort_s_type_characters
  (
    text,
    sl_type_vector,
    character_bucket_range_vector,
    text_distance_vector
  );
  return;
}

template
<
  typename text_type,
  typename text_distance_vector_type
>
auto calculate_text_distance_vector_boundary
(
  text_type const &text,
  text_distance_vector_type &text_distance_vector
)
{
  auto invalid_text_distance {std::size(text)};
  return std::stable_partition
  (
    std::begin(text_distance_vector),
    std::end(text_distance_vector),
    [&] (auto const &text_distance)
    {
      return (text_distance != invalid_text_distance);
    }
  );
}

template
<
  typename text_type,
  typename sl_type_vector_type,
  typename grammar_rule_begin_distance_vector_iterator_type,
  typename temporary_grammar_compressed_text_iterator_type
>
void filter_grammar_rule_begin_distance_vector_and_calculate_temporary_grammar_compressed_text
(
  text_type const &text,
  sl_type_vector_type const &sl_type_vector,
  grammar_rule_begin_distance_vector_iterator_type begin_distance_begin,
  grammar_rule_begin_distance_vector_iterator_type begin_distance_end,
  temporary_grammar_compressed_text_iterator_type temporary_grammar_compressed_text_begin
)
{
  using begin_distance_type = typename std::iterator_traits<grammar_rule_begin_distance_vector_iterator_type>::value_type;
  using non_terminal_type = typename std::iterator_traits<temporary_grammar_compressed_text_iterator_type>::value_type;
  auto invalid_text_distance {std::size(text)};
  non_terminal_type non_terminal {0};
  auto prev_grammar_rule_it {std::prev(std::end(text))};
  auto prev_sl_type_vector_it {std::prev(std::end(sl_type_vector))};
  for
  (
    auto begin_distance_it {begin_distance_begin};
    begin_distance_it != begin_distance_end;
    ++begin_distance_it
  )
  {
    begin_distance_type begin_distance {*begin_distance_it};
    auto grammar_rule_it {std::next(std::begin(text), begin_distance)};
    auto sl_type_vector_it {std::next(std::begin(sl_type_vector), begin_distance)};
    if
    (
      (*prev_grammar_rule_it == *grammar_rule_it)
      &&
      (*prev_sl_type_vector_it == *sl_type_vector_it)
    )
    {
      do
      {
        ++prev_grammar_rule_it;
        ++grammar_rule_it;
        ++prev_sl_type_vector_it;
        ++sl_type_vector_it;
      }
      while
      (
        !is_leftmost_s_type(prev_sl_type_vector_it)
        &&
        !is_leftmost_s_type(sl_type_vector_it)
        &&
        (*prev_grammar_rule_it == *grammar_rule_it)
        &&
        (*prev_sl_type_vector_it == *sl_type_vector_it)
      );
      if
      (
        is_leftmost_s_type(prev_sl_type_vector_it)
        &&
        is_leftmost_s_type(sl_type_vector_it)
      )
      {
        --non_terminal;
        *begin_distance_it = invalid_text_distance;
      }
    }
    *std::next(temporary_grammar_compressed_text_begin, (begin_distance + 1) / 2) = ++non_terminal;
    prev_grammar_rule_it = std::next(std::begin(text), begin_distance);
    prev_sl_type_vector_it = std::next(std::begin(sl_type_vector), begin_distance);
  }
  return;
}

template
<
  typename text_type,
  typename vector_iterator_type
>
auto collect_valid_entries
(
  text_type const &text,
  vector_iterator_type begin,
  vector_iterator_type end
)
{
  auto invalid_value {std::size(text)};
  return std::stable_partition
  (
    begin,
    end,
    [&] (auto const &value)
    {
      return (value != invalid_value);
    }
  );
}

template
<
  typename text_type,
  typename sl_type_vector_type,
  typename text_distance_vector_type
>
auto calculate_grammar_rule_begin_distance_vector_and_temporary_grammar_compressed_text
(
  text_type const &text,
  sl_type_vector_type const &sl_type_vector,
  text_distance_vector_type &text_distance_vector
)
{
  auto text_distance_vector_boundary {calculate_text_distance_vector_boundary(text, text_distance_vector)};
  auto grammar_rule_begin_distance_vector_begin {std::begin(text_distance_vector)};
  auto grammar_rule_begin_distance_vector_end {text_distance_vector_boundary};
  auto temporary_grammar_compressed_text_begin {text_distance_vector_boundary};
  auto temporary_grammar_compressed_text_end {std::end(text_distance_vector)};
  filter_grammar_rule_begin_distance_vector_and_calculate_temporary_grammar_compressed_text
  (
    text,
    sl_type_vector,
    grammar_rule_begin_distance_vector_begin,
    grammar_rule_begin_distance_vector_end,
    temporary_grammar_compressed_text_begin
  );
  grammar_rule_begin_distance_vector_end =
  collect_valid_entries
  (
    text,
    grammar_rule_begin_distance_vector_begin,
    grammar_rule_begin_distance_vector_end
  );
  temporary_grammar_compressed_text_end =
  collect_valid_entries
  (
    text,
    temporary_grammar_compressed_text_begin,
    temporary_grammar_compressed_text_end
  );
  return std::make_tuple
  (
    grammar_rule_begin_distance_vector_begin,
    grammar_rule_begin_distance_vector_end,
    temporary_grammar_compressed_text_begin,
    temporary_grammar_compressed_text_end
  );
}

template
<
  typename text_type,
  typename sl_type_vector_type,
  typename text_distance_vector_type
>
auto calculate_grammar
(
  text_type const &text,
  sl_type_vector_type &sl_type_vector,
  text_distance_vector_type &text_distance_vector
)
{
  calculate_sl_type_vector
  (
    text,
    sl_type_vector
  );
  induce_sort_grammar_rule
  (
    text,
    sl_type_vector,
    text_distance_vector
  );
  return calculate_grammar_rule_begin_distance_vector_and_temporary_grammar_compressed_text
  (
    text,
    sl_type_vector,
    text_distance_vector
  );
}

#endif
