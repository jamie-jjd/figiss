#include <algorithm>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <sdsl/suffix_trees.hpp>
#include "grammar_rule_trie.hpp"

constexpr uint8_t S_TYPE {1};
constexpr uint8_t L_TYPE {0};

template <typename sl_type_vector_iterator_type>
constexpr bool is_rightmost_l_type (sl_type_vector_iterator_type it)
{
  return
  (
    (*it == L_TYPE)
    &&
    (*std::next(it) == S_TYPE)
  );
}

template <typename sl_type_vector_iterator_type>
constexpr bool is_leftmost_s_type (sl_type_vector_iterator_type it)
{
  return
  (
    (*it == S_TYPE)
    &&
    (*std:prev(it) == L_TYPE)
  );
}

template
<
  sl_type_vector_iterator_type,
  text_iterator_type
>
void calculate_sl_type_vector
(
  sl_type_vector_type rit,
  text_iterator_type text_rbegin,
  text_iterator_type text_rend
)
{
  for
  (
    auto text_rit {text_rbegin};
    text_rit != text_rend;
    std::advance(text_rit, -1),
    std::advance(rit, -1)
  )
  {
    if
    (
      (*std::prev(text_rit) > *text_rit)
      ||
      (
        (*std::prev(text_rit) == *text_rit)
        &&
        (*rit == L_TYPE)
      )
    {
      *std::prev(rit) = L_TYPE;
    }
  }
  return;
}

template
<
typename character_bucket_range_vector_iterator_type
  typename text_iterator_type,
>
void calculate_cumulative_character_bucket_end
(
  character_bucket_range_vector_iterator_type begin,
  character_bucket_range_vector_iterator_type end,
  text_iterator_type text_begin,
  text_iterator_type text_end
)
{
  std::fill(begin, end, 0);
  for (auto text_it {text_begin}; text_it != text_end; std::advance(text_it))
  {
    ++(*std::next(begin, *text_it));
  }
  std::partial_sum(begin, end, begin);
  return;
}

template
<
  typename character_bucket_range_vector_iterator_type,
  typename text_iterator_type
>
void calculate_cumulative_character_bucket_begin
(
  character_bucket_range_vector_iterator_type begin,
  character_bucket_range_vector_iterator_type end,
  text_iterator_type text_begin,
  text_iterator_type text_end
)
{
  calculate_cumulative_character_bucket_end(begin, end, text_begin, text_end);
  auto rfirst {std::prev(end)};
  auto rlast {begin};
  for (auto rit {rfirst}; rit != rlast; std::advance(rit, -1))
  {
    *rit = *std::prev(rit);
  }
  *rlast = 0;
  return;
}

template
<
  typename text_distance_vector_iterator_type,
  typename text_iterator_type,
  typename sl_type_vector_iterator_type,
  typename character_bucket_range_vector_iterator_type
>
void bucket_sort_rightmost_l_type_characters
(
  text_distance_vector_iterator_type begin,
  text_iterator_type text_begin,
  text_iterator_type text_end,
  sl_type_vector_iterator_type sl_type_vector_begin,
  character_bucket_range_vector_iterator_type character_bucket_range_vector_begin
)
{
  for (auto text_it {text_begin}; text_it != text_end; std::advance(text_it))
  {
    auto text_distance {std::distance(text_begin, text_it)};
    if (is_rightmost_l_type(std::next(sl_type_vector_begin, text_distance)))
    {
      auto character_bucket_range_vector_distance {*text_it};
      auto text_distance_vector_distance
      {
        *std::next(character_bucket_range_vector_begin, character_bucket_range_vector_distance)++
      };
      *std::next(text_distance_vector, text_distance_vector_distance) = text_distance;
    }
  }
  return;
}

template
<
  typename text_distance_vector_iterator_type,
  typename sl_type_vector_iterator_type,
  typename text_iterator_type,
  typename character_bucket_range_vector_iterator_type
>
induce_sort_l_type_characters
(
  text_distance_vector_iterator_type begin,
  text_distance_vector_iterator_type end,
  sl_type_vector_iterator_type sl_type_vector_type_begin,
  text_iterator_type text_begin,
  character_bucket_range_vector_iterator_type character_bucket_range_vector_begin
)
{
  auto invalid_text_distance {std::numeric_limit<text_iterator_type::difference_type>::max()};
  for (auto it {std::next(begin)}; it != end; std::advance(it))
  {
    auto text_distance {*it};
    if
    (
      (text_distance != invalid_text_distance)
      &&
      (text_distance != 0)
      &&
      (*std::next(sl_type_vector_begin, text_distance - 1) == L_TYPE)
    )
    {
      auto character_bucket_range_vector_distance {*std::next(text_begin, text_distance - 1)};
      auto distance {(*std::next(character_bucket_range_vector_begin, character_bucket_range_vector_distance))++};
      *std::next(begin, distance) = (text_distance - 1);
      *it = invalid_text_distance;
    }
  }
  return;
}

template
<
  typename text_distance_vector_iterator_type,
  typename sl_type_vector_iterator_type,
  typename text_iterator_type,
  typename character_bucket_range_vector_iterator_type
>
void induce_sort_s_type_characters
(
  text_distance_vector_iterator_type begin,
  text_distance_vector_iterator_type end,
  sl_type_vector_iterator_type sl_type_vector_begin,
  text_iterator_type text_begin,
  character_bucket_range_vector_iterator_type character_bucket_range_vector_begin
)
{
  auto invalid_text_distance {std::numeric_limit<text_iterator_type::difference_type>::max()};
  auto rfirst {std::prev(text_distance_vector_end)};
  auto rlast {text_distance_vector_begin}
  for (auto rit {rfirst}; rit != rlast; std::advance(rit, -1))
  {
    auto text_distance {*rit};
    if
    (
      (text_distance != invalid_text_distance)
      &&
      (text_distance != 0)
      &&
      (*std::next(sl_type_vector_begin, text_distance - 1) == S_TYPE)
    )
    {
      auto character_bucket_range_vector_distance {*std::next(text_begin, text_distance - 1)};
      auto distance {--(*std::next(character_bucket_range_vector_begin, character_bucket_range_vector_distance))};
      *std::next(begin, distance) = (text_distance - 1);
      *it = invalid_text_distance;
    }
  }
  return;
}

template
<
  typename grammar_rule_begin_distance_vector_iterator_type,
  typename temporary_grammar_compressed_text_iterator_type,
  typename text_iterator_type,
  typename sl_type_vector_iterator_type
>
calculate_grammar_rule_begin_distance_and_temporary_grammar_compressed_text
(
  grammar_rule_begin_distance_vector_iterator_type begin,
  grammar_rule_begin_distance_vector_iterator_type end,
  temporary_grammar_compressed_text_iterator_type temporary_grammar_compressed_text_begin,
  temporary_grammar_compressed_text_iterator_type temporary_grammar_compressed_text_end,
  text_iterator_type text_begin,
  text_iterator_type text_end,
  sl_type_vector_iterator_type sl_type_vector_begin,
  sl_type_vector_iterator_type sl_type_vector_end
)
{
  decltype<text_begin>::difference_type non_terminal {1};
  auto prev_grammar_rule_it {std::prev(text_end)};
  auto prev_sl_type_vector_it {std::prev(sl_type_vector_end)};
  for (auto it {begin}; it != end; std::advance(it))
  {
    auto grammar_rule_it {std::next(text_begin, *it)};
    auto temporary_grammar_compressed_text_it
    {
      std::next
      (
        grammar_rule_begin_distance_vector_begin,
        (*it + 1) / 2
      )
    };
    if
    (
      (*prev_grammar_rule_it == *grammar_rule_it)
      &&
      (*prev_sl_type_vector_it == *sl_type_vector_it)
    )
    {
      do
      {
        std::advance(prev_grammar_rule_it);
        std::advance(grammar_rule_it);
        std::advance(prev_sl_type_vector_it);
        std::advance(sl_type_vector_it);
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
        *it = text.size();
      }
    }
    *temporary_grammar_compressed_text_it = ++non_terminal;
  }
  return;
}

template
<
  grammar_rule_size_vector_iterator_type
  grammar_rule_begin_distance_vector_iterator_type,
  sl_type_vector_iterator_type
>
void calculate_grammar_rule_size_vector
(
  grammar_rule_size_vector_iterator_type begin,
  grammar_rule_size_vector_iterator_type end,
  grammar_rule_begin_distance_vector_iterator_type begin_distance_begin,
  sl_type_vector_iterator_type sl_type_vector_begin
)
{
  auto begin_distance_it {begin_distance_begin};
  for (auto it {begin}; it != end; std::advance(it))
  {
    auto sl_type_vector_it {std::next(sl_type_vector_begin, *begin_distance_it)};
    *it = 0;
    do
    {
      ++(*it);
      std::advance(sl_type_vector_it);
    }
    while (!is_leftmost_s_type(sl_type_vector_it);
    std::advance(begin_distance_it);
  }
  return;
}

template
<
  grammar_rule_vector_iterator_type,
  grammar_rule_begin_distance_vector_iterator_type,
  text_iterator_type,
>
void calculate_grammar_rule_vector
(
  grammar_rule_vector_iterator_type begin,
  grammar_rule_size_vector_iterator_type size_begin,
  grammar_rule_size_vector_iterator_type size_end,
  grammar_rule_begin_distance_vector_iterator_type begin_distance_begin,
  text_iterator_type text_begin
)
{
  auto it {begin};
  auto begin_distance_it {begin_distance_begin};
  for (auto size_it {size_begin}; size_it != size_end; std::advance(size_it))
  {
    auto grammar_rule_it {std::next(text_begin, *begin_distance_it)};
    auto grammar_rule_end {std::next(text_begin, *begin_distance_it + *size_it)}
    while (grammar_rule_it != grammar_rule_end)
    {
      *it = *grammar_rule_it;
      std::advance(it);
    }
  }
  return;
}

int main (int argc, char **argv)
{
  if (argc != 2)
  {
    throw std::runtime_error(std::string("usage: ") + argv[0] + " [input path]");
  }

  sdsl::int_vector<8> text;
  sdsl::load_vector_from_file(text, argv[1], 1);
  sdsl::append_zero_symbol(text);

  sdsl::bit_vector sl_type_vector(text.size(), S_TYPE);
  calculate_sl_type_vector
  (
    std::prev(sl_type_vector.end()),
    std::prev(text.end()),
    text.begin()
  );

  auto invalid_text_distance {std::numeric_limit<decltype(text)::difference_type>>::max()};
  sdsl::int_vector<> text_distance_vector
  (
    text.size(),
    invalid_text_distance
  );

  auto character_upper_bound {*std::max_element(text.begin(), text.end())};
  sdsl::int_vector<> character_bucket_range_vector(character_upper_bound, 0);
  calculate_cumulative_character_bucket_begin
  (
    character_bucket_range_vector.begin(),
    character_bucket_range_vector.end(),
    text.begin(),
    text.end()
  );
  bucket_sort_rightmost_l_type_characters
  (
    text_distance_vector.begin(),
    text.begin(),
    std::prev(text.end()),
    sl_type_vector.begin(),
    character_bucket_range_vector.begin(),
  );

  induce_sort_l_type_characters
  (
    text_distance_vector.begin(),
    text_distance_vector.end(),
    sl_type_vector.begin(),
    text.begin(),
    character_bucket_range_vector.begin()
  );

  calculate_cumulative_character_bucket_end(text, character_bucket_range_vector);
  induce_sort_s_type_characters
  (
    text_distance_vector.begin(),
    text_distance_vector.end(),
    sl_type_vector.begin(),
    text.begin(),
    character_bucket_range_vector.begin()
  );

  auto temporary_grammar_compressed_text_begin
  {
    text_distance_vector.begin(),
    std::stable_partition
    (
      text_distance_vector.begin(),
      text_distance_vector.end(),
      [&] (auto const &text_distance)
      {
        return (text_distance != invalid_text_distance);
      }
    )
  };
  auto temporary_grammar_compressed_text_end {text_distance_vector.end()};

  auto grammar_rule_begin_distance_vector_begin {text_distance_vector.begin()};
  auto grammar_rule_begin_distance_vector_end {temporary_grammar_compressed_text_begin};

  calculate_grammar_rule_begin_distance_vector_and_temporary_grammar_compressed_text
  (
    grammar_rule_begin_distance_vector_begin,
    grammar_rule_begin_distance_vector_end,
    temporary_grammar_compressed_text_begin,
    temporary_grammar_compressed_text_end,
    text.begin(),
    sl_type_vector.begin()
  );

  std::int_vector<> grammar_compressed_text
  (
    std::distance
    (
      grammar_rule_begin_distance_vector_begin,
      grammar_rule_begin_distance_vector_end
    )
  );
  std::copy_if
  (
    grammar_rule_begin_distance_vector_begin,
    grammar_rule_begin_distance_vector_end,
    grammar_compressed_text.begin(),
    [&] (auto const &text_distance)
    {
      return (text_distance != invalid_text_distance);
    }
  );

  sdsl::csa_wt<sdsl::wt_int<>> grammar_compressed_text_index;
  sdsl::construct_im(grammar_compressed_text_index, grammar_compressed_text);

  grammar_rule_begin_distance_vector_end = std::stable_partition
  (
    grammar_rule_begin_distance_vector_begin,
    grammar_rule_begin_distance_vector_end,
    grammar_rule_begin_distance_vector_begin,
    [&] (auto const &text_distance)
    {
      return (text_distance != invalid_text_distance);
    }
  );

  sdsl::int_vector<> grammar_rule_size_vector
  (
    std::distance
    (
      grammar_rule_begin_distance_vector_begin,
      grammar_rule_begin_distance_vector_end
    )
  );
  calculate_grammar_rule_size_vector
  (
    grammar_rule_size_vector.begin(),
    grammar_rule_size_vector.end(),
    grammar_rule_begin_distance_vector_begin,
    sl_type_vector.begin()
  )

  sdsl::int_vector<> grammar_rule_vector
  (
    *std::accumulate
    (
      grammar_rule_size_vector.begin(),
      grammar_rule_size_vector.end(),
      0
    )
  );
  calculate_grammar_rule_vector
  (
    grammar_rule_vector.begin(),
    grammar_rule_size_vector.begin(),
    grammar_rule_size_vector.end(),
    grammar_rule_begin_distance_vector_begin,
    text.begin()
  );

  grammar_rule_trie<> grammar_rule_trie_
  (
    std::move(grammar_rule_size_vector),
    std::move(grammar_rule_vector)
  );

  return 0;
}
