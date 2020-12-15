#include <algorithm>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <sdsl/suffix_trees.hpp>

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
  typename text_distance_for_grammar_rule_begin_vector_iterator_type,
  typename temporary_grammar_compressed_text_iterator_type,
  typename text_iterator_type,
  typename sl_type_vector_iterator_type
>
calculate_text_distance_for_grammar_rule_begin_and_temporary_grammar_compressed_text
(
  text_distance_for_grammar_rule_begin_vector_iterator_type begin,
  text_distance_for_grammar_rule_begin_vector_iterator_type end,
  temporary_grammar_compressed_text_iterator_type temporary_grammar_rule_vector_begin,
  temporary_grammar_compressed_text_iterator_type temporary_grammar_rule_vector_end,
  text_iterator_type text_begin,
  text_iterator_type text_end,
  sl_type_vector_iterator_type sl_type_vector_begin,
  sl_type_vector_iterator_type sl_type_vector_end
)
{
  decltype<text_begin>::difference_type non_terminal {1};
  auto prev_grammar_rule_it {std::prev(text_end)};
  auto prev_sl_type_it {std::prev(sl_type_vector_end)};
  for (auto it {begin}; it != end; std::advance(it))
  {
    auto grammar_rule_it {std::next(text_begin, *it)};
    if
    (
      (*prev_grammar_rule_it == *grammar_rule_it)
      &&
      (*prev_sl_type_it == *sl_type_it)
    )
    {
      std::advance
      while
      (
        !is_leftmost_s_type(sl_type_vector, distance_prev_grammar_rule)
        &&
        !is_leftmost_s_type(sl_type_vector, distance_grammar_rule)
        &&
        (text[distance_prev_grammar_rule] == text[distance_grammar_rule])
        &&
        (sl_type_vector[distance_prev_grammar_rule] == sl_type_vector[distance_grammar_rule])
      )
      {
        ++distance_prev_grammar_rule;
        ++distance_grammar_rule;
      }
      if // same grammar rule
      (
        is_leftmost_s_type(sl_type_vector, distance_prev_grammar_rule)
        &&
        is_leftmost_s_type(sl_type_vector, distance_grammar_rule)
      )
      {
        --non_terminal;
        text_distance_vector[rank] = text.size();
      }
    }
    auto temporary_grammar_rule_vector_distance { + ((begin_grammar_rule + 1) / 2)};
    *std::next(temporary_grammar_rule_vector_begin, temporary_grammar_rule_vector_distance) = ++non_terminal;

  }
  return;
}

template
<
  text_distance_for_grammar_rule_begin_vector_iterator_type,
  sl_type_vector_iterator_type
>
auto calculate_concatenated_grammar_rule_size
(
  text_distance_for_grammar_rule_begin_vector_iterator_type begin,
  text_distance_for_grammar_rule_begin_vector_iterator_type end,
  text_iterator_type text_begin,
  sl_type_vector_iterator_type sl_type_vector_begin
)
{
  decltype(text_begin)::difference_type size {0};
  for (auto it {begin}; it != end; std::advance(it))
  {
    auto sl_type_it {std::next(sl_type_vector_begin, *it + 1)};
    while (!is_leftmost_s_type(sl_type_it)
    {
      ++sl_type_it;
    }
    size += std::distance(std::next(sl_type_begin, *it), sl_type_it);
  }
  return size;
}

template
<
  character_type,
  size_type
>
struct grammar_rule_trie_node
{
  using range_type = std::pair<size_type, size_type>;

  range_type edge_range;
  std::map<character_type, std::shared_ptr<grammar_rule_trie_node>> branch;
  range_type non_terminal_range;

  grammar_rule_trie_node
  (
    range_type const &edge_range_
  ) : edge_range {edge_range_}
  {
  }

  grammar_rule_trie_node
  (
    range_type const &edge_range_
    size_type const &non_terminal
  ) : edge_range {edge_range_},
      non_terminal_range {{non_terminal, non_terminal}}
  {
  }
};

template
<
  typename grammar_rule_trie_node_type,
  typename grammar_rule_vector_type,
  typename cumulative_distance_type,
  typename non_terminal_type
>
void insert_grammar_rule
(
  grammar_rule_trie_node_type root,
  grammar_rule_vector_type const &grammar_rule_vector,
  cumulative_distance_type const begin,
  cumulative_distance_type const end,
  non_terminal_type const non_terminal
)
{
  auto node {root};
  auto distance {begin};
  while (true)
  {
    auto character {grammar_rule_vector[distance]};
    auto child {node->branch[character]};
    if (child == node->branch.end())
    {
      node->branch[character] = std::make_shared<grammar_rule_trie_node_type>({distance, end}, non_terminal);
      return;
    }
    auto distance_edge {std::get<0>(child->edge_range)};
    auto end_edge {std::get<1>(child->edge_range)};
    while
    (
      (distance_edge < end_edge)
      &&
      (distance < end)
      &&
      (grammar_rule_vector[distance_edge] == grammar_rule_vector[distance])
    )
    {
      ++distance_child;
      ++distance;
    }
    if (distance_edge < end_edge)
    {
      return;
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

  /* read text from file */
  sdsl::int_vector<8> text;
  sdsl::load_vector_from_file(text, argv[1], 1);
  sdsl::append_zero_symbol(text);

  /* initialize SL-type of characters of text */
  sdsl::bit_vector sl_type_vector(text.size(), S_TYPE);
  calculate_sl_type_vector
  (
    std::prev(sl_type_vector.end()),
    std::prev(text.end()),
    text.begin()
  )

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

  auto text_distance_for_grammar_rule_begin_vector_begin {text_distance_vector.begin()};
  auto text_distance_for_grammar_rule_begin_vector_end {temporary_grammar_compressed_text_begin};

  calculate_text_distance_for_grammar_rule_begin_and_temporary_grammar_compressed_text
  (
    text_distance_for_grammar_rule_begin_vector_begin,
    text_distance_for_grammar_rule_begin_vector_end,
    temporary_grammar_compressed_text_begin,
    temporary_grammar_compressed_text_end,
    text.begin(),
    sl_type_vector.begin()
  );

  std::int_vector<> grammar_compressed_text
                    (
                      std::distance
                      (
                        text_distance_for_grammar_rule_begin_vector_begin,
                        text_distance_for_grammar_rule_begin_vector_end
                      )
                    );
  std::copy_if
  (
    temporary_grammar_rule_vector_begin,
    temporary_grammar_rule_vector_end,
    grammar_compressed_text.begin(),
    [&] (auto const &text_distance)
    {
      return (text_distance != invalid_text_distance);
    }
  );

  sdsl::csa_wt<sdsl::wt_int<>> grammar_compressed_text_fm_index;
  sdsl::construct_im(grammar_compressed_text_fm_index, grammar_compressed_text);

  text_distance_for_grammar_rule_begin_vector_end = std::stable_partition
                                            (
                                              temporary_grammar_rule_vector_begin,
                                              temporary_grammar_rule_vector_end,
                                              temporary_grammar_rule_vector_begin,
                                              [&] (auto const &text_distance)
                                              {
                                                return (text_distance != invalid_text_distance);
                                              }
                                            );

  sdsl::int_vector<> contatenated_grammar_rule_vector
  (
    calculate_concatenated_grammar_rule_size
    (
      text_distance_for_grammar_rule_begin_vector_begin,
      text_distance_for_grammar_rule_begin_vector_end,
      text.begin(),
      sl_type_vector.begin()
    )
  );

  /* constrcut compact trie of grammar rule set */
  std::shared_ptr<grammar_rule_trie_node> grammar_rule_trie_root;
  std::shared_ptr<grammar_rule_trie_node> reverse_grammar_rule_trie_root;
  auto cumulative_begin {0ULL};
  for (auto rank {0ULL}; rank < end_grammar_rule_begin_vector; ++rank)
  {
    auto begin_grammar_rule {static_cast<uint64_t>(text_distance_vector[rank])};
    grammar_rule_vector[cumulative_begin] = text[begin_grammar_rule];
    auto offset {1ULL};
    while (!is_leftmost_s_type(sl_type_vector, begin_grammar_rule + offset))
    {
      grammar_rule_vector[cumulative_begin + offset] = text[begin_grammar_rule + offset];
      ++offset;
    }
    grammar_rule_trie_node = insert_grammar_rule
    (
      grammar_rule_trie_root,
      grammar_rule_vector,
      cumulative_begin,
      (cumulative_begin + offset),
      (rank + 1)
    );
    reverse_grammar_rule_trie_root = insert_reverse_grammar_rule
    (
      reverse_grammar_rule_trie_root,
      grammar_rule_vector,
      cumulative_begin,
      (cumulative_begin + offset),
      (rank + 1)
    )
    cumulative_begin += offset;
  }

  return 0;
}
