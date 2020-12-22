#ifndef GRAMMAR_COMPRESSED_INDEX_HPP_
#define GRAMMAR_COMPRESSED_INDEX_HPP_

#include <map>
#include <memory>

#include <sdsl/wavelet_trees.hpp>

constexpr uint8_t S {1};
constexpr uint8_t L {0};

template <typename sl_types_iterator_type>
constexpr bool is_rightmost_l (sl_types_iterator_type it)
{
  return
  (
    (*it == L)
    &&
    (*std::next(it) == S)
  );
}

template <typename sl_types_iterator_type>
constexpr bool is_leftmost_s (sl_types_iterator_type it)
{
  return
  (
    (*it == S)
    &&
    (*std::prev(it) == L)
  );
}

template
<
  typename string_type,
  typename sl_types_type
>
void calculate_sl_types
(
  string_type const &string,
  sl_types_type &sl_types
)
{
  sl_types.resize(string.size());
  std::fill
  (
    std::begin(sl_types),
    std::end(sl_types),
    S
  );
  auto string_rit {std::prev(std::end(string))};
  auto string_rend {std::begin(string)};
  auto rit {std::prev(std::end(sl_types))};
  while (string_rit != string_rend)
  {
    if
    (
      (*std::prev(string_rit) > *string_rit)
      ||
      (
        (*std::prev(string_rit) == *string_rit)
        &&
        (*rit == L)
      )
    )
    {
      *std::prev(rit) = L;
    }
    --rit;
    --string_rit;
  }
  return;
}

template
<
  typename text_type,
  typename character_buckets_type
>
void calculate_character_bucket_ends
(
  text_type const &text,
  character_buckets_type &character_buckets
)
{
  auto begin {std::begin(character_buckets)};
  auto end {std::end(character_buckets)};
  std::fill(begin, end, 0);
  for (auto const &character : text)
  {
    ++character_buckets[character];
  }
  std::partial_sum(begin, end, begin);
  return;
}

template
<
  typename text_type,
  typename character_buckets_type
>
void calculate_character_bucket_begins
(
  text_type const &text,
  character_buckets_type &character_buckets
)
{
  calculate_character_bucket_ends(text, character_buckets);
  auto rfirst {std::prev(std::end(character_buckets))};
  auto rlast {std::begin(character_buckets)};
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
  typename sl_types_type,
  typename character_buckets_type,
  typename text_dists_type
>
void bucket_sort_rightmost_l_type_characters
(
  text_type const &text,
  sl_types_type const &sl_types,
  character_buckets_type &character_buckets,
  text_dists_type &text_dists
)
{
  for (auto text_it {std::begin(text)}; text_it != std::end(text); ++text_it)
  {
    auto text_dist {std::distance(std::begin(text), text_it)};
    if (is_rightmost_l(std::next(std::begin(sl_types), text_dist)))
    {
      auto bucket_begin_dist {character_buckets[*text_it]++};
      text_dists[bucket_begin_dist] = text_dist;
    }
  }
  return;
}

template
<
  typename text_type,
  typename sl_types_type,
  typename character_buckets_type,
  typename text_dists_type
>
void induce_sort_l_type_characters
(
  text_type const &text,
  sl_types_type const &sl_types,
  character_buckets_type &character_buckets,
  text_dists_type &text_dists
)
{
  auto text_dist_begin {std::begin(text_dists)};
  auto text_dist_end {std::end(text_dists)};
  auto invalid_text_dist {std::size(text)};
  for
  (
    auto text_dist_it {std::next(text_dist_begin)};
    text_dist_it != text_dist_end;
    ++text_dist_it
  )
  {
    auto text_dist {*text_dist_it};
    if
    (
      (text_dist != invalid_text_dist)
      &&
      (text_dist != 0)
      &&
      (sl_types[text_dist - 1] == L)
    )
    {
      auto bucket_begin_dist {character_buckets[text[text_dist - 1]]++};
      text_dists[bucket_begin_dist] = (text_dist - 1);
      *text_dist_it = invalid_text_dist;
    }
  }
  return;
}

template
<
  typename text_type,
  typename sl_types_type,
  typename character_buckets_type,
  typename text_dists_type
>
void induce_sort_s_type_characters
(
  text_type const &text,
  sl_types_type const &sl_types,
  character_buckets_type &character_buckets,
  text_dists_type &text_dists
)
{
  auto text_dist_rfirst {std::prev(std::end(text_dists))};
  auto text_dist_rlast {std::begin(text_dists)};
  auto invalid_text_dist {std::size(text)};
  for (auto text_dist_rit {text_dist_rfirst}; text_dist_rit != text_dist_rlast; --text_dist_rit)
  {
    auto text_dist {*text_dist_rit};
    if
    (
      (text_dist != invalid_text_dist)
      &&
      (text_dist != 0)
      &&
      (sl_types[text_dist - 1] == S)
    )
    {
      auto bucket_end_dist {--character_buckets[text[text_dist - 1]]};
      text_dists[bucket_end_dist] = (text_dist - 1);
      *text_dist_rit = invalid_text_dist;
    }
  }
  return;
}

template
<
  typename text_type,
  typename text_dists_type
>
auto calculate_text_dists_boundary
(
  text_type const &text,
  text_dists_type &text_dists
)
{
  auto invalid_text_dist {std::size(text)};
  return std::stable_partition
  (
    std::begin(text_dists),
    std::end(text_dists),
    [&] (auto const &text_dist)
    {
      return (text_dist != invalid_text_dist);
    }
  );
}

template
<
  typename text_type,
  typename sl_types_type,
  typename grammar_rule_begin_dists_iterator_type,
  typename temporary_gc_text_iterator_type
>
void filter_grammar_rule_begin_dists_and_calculate_temporary_gc_text
(
  text_type const &text,
  sl_types_type const &sl_types,
  grammar_rule_begin_dists_iterator_type begin_dists_begin,
  grammar_rule_begin_dists_iterator_type begin_dists_end,
  temporary_gc_text_iterator_type temporary_gc_text_begin
)
{
  auto invalid_text_dist {std::size(text)};
  uint64_t lex_rank {0};
  auto prev_grammar_rule_it {std::prev(std::end(text))};
  auto prev_sl_types_it {std::prev(std::end(sl_types))};
  auto begin_dists_it {begin_dists_begin};
  while (begin_dists_it != begin_dists_end)
  {
    uint64_t begin_dist {*begin_dists_it};
    auto grammar_rule_it {std::next(std::begin(text), begin_dist)};
    auto sl_types_it {std::next(std::begin(sl_types), begin_dist)};
    if
    (
      (*prev_grammar_rule_it == *grammar_rule_it)
      &&
      (*prev_sl_types_it == *sl_types_it)
    )
    {
      do
      {
        ++prev_grammar_rule_it;
        ++grammar_rule_it;
        ++prev_sl_types_it;
        ++sl_types_it;
      }
      while
      (
        !is_leftmost_s(prev_sl_types_it)
        &&
        !is_leftmost_s(sl_types_it)
        &&
        (*prev_grammar_rule_it == *grammar_rule_it)
        &&
        (*prev_sl_types_it == *sl_types_it)
      );
      if
      (
        is_leftmost_s(prev_sl_types_it)
        &&
        is_leftmost_s(sl_types_it)
      )
      {
        --lex_rank;
        *begin_dists_it = invalid_text_dist;
      }
    }
    *std::next(temporary_gc_text_begin, (begin_dist + 1) / 2) = ++lex_rank;
    prev_grammar_rule_it = std::next(std::begin(text), begin_dist);
    prev_sl_types_it = std::next(std::begin(sl_types), begin_dist);
    ++begin_dists_it;
  }
  return;
}

template
<
  typename text_type,
  typename random_access_iterator_type
>
auto collect_valid_entries
(
  text_type const &text,
  random_access_iterator_type begin,
  random_access_iterator_type end
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
  typename grammar_rule_sizes_type,
  typename sl_types_type,
  typename grammar_rule_begin_dists_iterator_type
>
void calculate_grammar_rule_sizes
(
  grammar_rule_sizes_type &grammar_rule_sizes,
  sl_types_type const &sl_types,
  grammar_rule_begin_dists_iterator_type begin_dists_begin,
  grammar_rule_begin_dists_iterator_type begin_dists_end
)
{
  grammar_rule_sizes.resize(std::distance(begin_dists_begin, begin_dists_end));
  auto grammar_rule_sizes_it {std::begin(grammar_rule_sizes)};
  auto begin_dists_it {begin_dists_begin};
  while (begin_dists_it != begin_dists_end)
  {
    auto sl_types_first {std::next(std::begin(sl_types), *begin_dists_it)};
    auto sl_types_it {sl_types_first};
    do
    {
      ++sl_types_it;
    }
    while (!is_leftmost_s(sl_types_it));
    *grammar_rule_sizes_it = std::distance(sl_types_first, sl_types_it);
    ++grammar_rule_sizes_it;
    ++begin_dists_it;
  }
  return;
}

struct trie_node
{
  int32_t begin_dist;
  int32_t end_dist;
  std::map<uint8_t, std::shared_ptr<trie_node>> branch;
  uint32_t leftmost_rank;
  uint32_t rightmost_rank;

  trie_node () = default;

  trie_node
  (
    int32_t begin_dist_,
    int32_t end_dist_,
    uint32_t rank
  )
  : begin_dist {begin_dist_},
    end_dist {end_dist_},
    leftmost_rank {rank},
    rightmost_rank {rank}
  {
  }
};

template
<
  typename text_type,
  typename grammar_rule_sizes_type,
  typename grammar_rules_type,
  typename grammar_rule_begin_dists_iterator_type
>
void calculate_grammar_rules
(
  text_type const &text,
  grammar_rule_sizes_type const &sizes,
  grammar_rules_type &rules,
  grammar_rule_begin_dists_iterator_type begin_dists_begin
)
{
  rules.resize(std::accumulate(std::begin(sizes), std::end(sizes), 0));
  auto rules_it {std::begin(rules)};
  auto begin_dists_it {begin_dists_begin};
  auto sizes_it {std::begin(sizes)};
  auto sizes_end {std::end(sizes)};
  while (sizes_it != sizes_end)
  {
    auto rule_it {std::next(std::begin(text), *begin_dists_it)};
    auto rule_end {std::next(std::begin(text), *begin_dists_it + *sizes_it)};
    while (rule_it != rule_end)
    {
      *rules_it = *rule_it;
      ++rules_it;
      ++rule_it;
    }
    ++begin_dists_it;
    ++sizes_it;
  }
  return;
}

// template
// <
//   typename text_type,
//   typename grammar_rule_begin_dists_iterator_type,
//   typename grammar_rule_sizes_type,
//   typename grammar_rules_type,
//   typename trie_node_pointer_type
// >
// void insert_grammar_rules
// (
//   text_type const &text,
//   grammar_rule_begin_dists_iterator_type grammar_rule_begin_dists_begin,
//   grammar_rule_begin_dists_iterator_type grammar_rule_begin_dists_end,
//   grammar_rule_sizes_type const &grammar_rule_sizes,
//   grammar_rules_uint32_t &grammar_rules,
//   trie_node_pointer_type lex_trie_root,
//   trie_node_pointer_type colex_trie_root
// )
// {
//   grammar_rules.resize(*std::accumulate(std::begin(grammar_rule_sizes), std::end(grammar_rule_sizes), 0));
//
//   return;
// }

struct gc_index
{
  sdsl::int_vector<8> grammar_rules;
  std::shared_ptr<trie_node> lex_trie_root;
  std::shared_ptr<trie_node> colex_trie_root;
  sdsl::int_vector<> lex_gc_character_begins;
  sdsl::wt_int<> lex_gc_bwt;
  sdsl::wt_int<> colex_gc_bwt;

  template <typename text_type>
  gc_index (text_type const &text)
  {
    sdsl::bit_vector sl_types;
    calculate_sl_types(text, sl_types);

    sdsl::int_vector<> text_dists(text.size());
    auto invalid_text_dist {text.size()};
    std::fill(std::begin(text_dists), std::end(text_dists), invalid_text_dist);

    sdsl::int_vector<> character_buckets(256);
    calculate_character_bucket_begins(text, character_buckets);
    bucket_sort_rightmost_l_type_characters(text, sl_types, character_buckets, text_dists);
    induce_sort_l_type_characters(text, sl_types, character_buckets, text_dists);
    calculate_character_bucket_ends(text, character_buckets);
    induce_sort_s_type_characters(text, sl_types, character_buckets, text_dists);

    auto text_dists_boundary {calculate_text_dists_boundary(text, text_dists)};
    auto grammar_rule_begin_dists_begin {std::begin(text_dists)};
    auto grammar_rule_begin_dists_end {text_dists_boundary};
    auto temporary_gc_text_begin {text_dists_boundary};
    auto temporary_gc_text_end {std::end(text_dists)};

    filter_grammar_rule_begin_dists_and_calculate_temporary_gc_text
    (
      text,
      sl_types,
      grammar_rule_begin_dists_begin,
      grammar_rule_begin_dists_end,
      temporary_gc_text_begin
    );

    grammar_rule_begin_dists_end =
    collect_valid_entries
    (
      text,
      grammar_rule_begin_dists_begin,
      grammar_rule_begin_dists_end
    );

    sdsl::int_vector<> grammar_rule_sizes;
    calculate_grammar_rule_sizes
    (
      grammar_rule_sizes,
      sl_types,
      grammar_rule_begin_dists_begin,
      grammar_rule_begin_dists_end
    );

    calculate_grammar_rules
    (
      text,
      grammar_rule_sizes,
      grammar_rules,
      grammar_rule_begin_dists_begin
    );

    temporary_gc_text_end =
    collect_valid_entries
    (
      text,
      temporary_gc_text_begin,
      temporary_gc_text_end
    );

  }
};

#endif
