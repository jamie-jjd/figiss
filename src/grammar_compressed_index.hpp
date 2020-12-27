#ifndef GRAMMAR_COMPRESSED_INDEX_HPP_
#define GRAMMAR_COMPRESSED_INDEX_HPP_

#include <map>
#include <memory>

#include <sdsl/suffix_trees.hpp>

template <typename T>
void show_type (T)
{
  std::cout << __PRETTY_FUNCTION__ << '\n';
  return;
}

namespace gci
{

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
  typename text_type,
  typename sl_types_type
>
void calculate_sl_types
(
  text_type const &text,
  sl_types_type &sl_types
)
{
  sl_types.resize(std::size(text));
  std::fill(std::begin(sl_types), std::end(sl_types), S);
  auto text_rit {std::prev(std::end(text))};
  auto text_rlast {std::begin(text)};
  auto sl_types_rit {std::prev(std::end(sl_types))};
  while (text_rit != text_rlast)
  {
    if
    (
      (*std::prev(text_rit) > *text_rit)
      ||
      (
        (*std::prev(text_rit) == *text_rit)
        &&
        (*sl_types_rit == L)
      )
    )
    {
      *std::prev(sl_types_rit) = L;
    }
    --sl_types_rit;
    --text_rit;
  }
  return;
}

template
<
  typename text_type,
  typename character_bucket_dists_type
>
void calculate_character_bucket_end_dists
(
  text_type const &text,
  character_bucket_dists_type &character_bucket_dists
)
{
  auto begin {std::begin(character_bucket_dists)};
  auto end {std::end(character_bucket_dists)};
  std::fill(begin, end, 0);
  for (auto const &character : text)
  {
    ++character_bucket_dists[character];
  }
  std::partial_sum(begin, end, begin);
  return;
}

template
<
  typename text_type,
  typename character_bucket_dists_type
>
void calculate_character_bucket_begin_dists
(
  text_type const &text,
  character_bucket_dists_type &character_bucket_dists
)
{
  calculate_character_bucket_end_dists(text, character_bucket_dists);
  auto rit {std::prev(std::end(character_bucket_dists))};
  auto rlast {std::begin(character_bucket_dists)};
  while (rit != rlast)
  {
    *rit = *std::prev(rit);
    --rit;
  }
  *rlast = 0;
  return;
}

template
<
  typename text_type,
  typename sl_types_type,
  typename character_bucket_dists_type,
  typename text_dists_type
>
void bucket_sort_rightmost_l_type_characters
(
  text_type const &text,
  sl_types_type const &sl_types,
  character_bucket_dists_type &character_bucket_dists,
  text_dists_type &text_dists
)
{
  auto text_it {std::begin(text)};
  auto text_end {std::end(text)};
  while (text_it != text_end)
  {
    auto text_dist {std::distance(std::begin(text), text_it)};
    if (is_rightmost_l(std::next(std::begin(sl_types), text_dist)))
    {
      text_dists[character_bucket_dists[*text_it]++] = text_dist;
    }
    ++text_it;
  }
  return;
}

template
<
  typename text_type,
  typename sl_types_type,
  typename character_bucket_dists_type,
  typename text_dists_type
>
void induce_sort_l_type_characters
(
  text_type const &text,
  sl_types_type const &sl_types,
  character_bucket_dists_type &character_bucket_dists,
  text_dists_type &text_dists
)
{
  auto text_dists_it {std::next(std::begin(text_dists))};
  auto text_dists_last {std::end(text_dists)};
  auto invalid_text_dist {std::size(text)};
  while (text_dists_it != text_dists_last)
  {
    auto text_dist {*text_dists_it};
    if
    (
      (text_dist != invalid_text_dist)
      &&
      (text_dist != 0)
      &&
      (sl_types[text_dist - 1] == L)
    )
    {
      text_dists[character_bucket_dists[text[text_dist - 1]]++] = (text_dist - 1);
      *text_dists_it = invalid_text_dist;
    }
    ++text_dists_it;
  }
  return;
}

template
<
  typename text_type,
  typename sl_types_type,
  typename character_bucket_dists_type,
  typename text_dists_type
>
void induce_sort_s_type_characters
(
  text_type const &text,
  sl_types_type const &sl_types,
  character_bucket_dists_type &character_bucket_dists,
  text_dists_type &text_dists
)
{
  auto text_dists_rit {std::prev(std::end(text_dists))};
  auto text_dists_rlast {std::begin(text_dists)};
  auto invalid_text_dist {std::size(text)};
  while (text_dists_rit != text_dists_rlast)
  {
    auto text_dist {*text_dists_rit};
    if
    (
      (text_dist != invalid_text_dist)
      &&
      (text_dist != 0)
      &&
      (sl_types[text_dist - 1] == S)
    )
    {
      text_dists[--character_bucket_dists[text[text_dist - 1]]] = (text_dist - 1);
      *text_dists_rit = invalid_text_dist;
    }
    --text_dists_rit;
  }
  return;
}

template <typename random_access_iterator_type>
auto collect_valid_entries
(
  random_access_iterator_type first,
  random_access_iterator_type last,
  uint64_t const invalid_value
)
{
  return std::stable_partition
  (
    first,
    last,
    [&] (auto const &value)
    {
      return (value != invalid_value);
    }
  );
}

template
<
  typename text_type,
  typename sl_types_type,
  typename grammar_rule_begin_dists_iterator_type,
  typename temp_gc_text_iterator_type
>
void calculate_grammar_rule_begin_dists_and_temp_gc_text
(
  text_type const &text,
  sl_types_type const &sl_types,
  grammar_rule_begin_dists_iterator_type begin_dists_begin,
  grammar_rule_begin_dists_iterator_type &begin_dists_end,
  temp_gc_text_iterator_type temp_gc_text_begin,
  temp_gc_text_iterator_type &temp_gc_text_end
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
    *std::next(temp_gc_text_begin, (begin_dist + 1) / 2) = ++lex_rank;
    prev_grammar_rule_it = std::next(std::begin(text), begin_dist);
    prev_sl_types_it = std::next(std::begin(sl_types), begin_dist);
    ++begin_dists_it;
  }
  begin_dists_end = collect_valid_entries(begin_dists_begin, begin_dists_end, invalid_text_dist);
  temp_gc_text_end = collect_valid_entries(temp_gc_text_begin, temp_gc_text_end, invalid_text_dist);
  return;
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
  int32_t edge_begin_dist;
  int32_t edge_end_dist;
  std::map<uint8_t, std::shared_ptr<trie_node>> branches;
  uint32_t leftmost_rank;
  uint32_t rightmost_rank;

  trie_node () = default;

  trie_node
  (
    int32_t edge_begin_dist_,
    int32_t edge_end_dist_,
    uint32_t rank
  )
  : edge_begin_dist {edge_begin_dist_},
    edge_end_dist {edge_end_dist_},
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

template
<
  typename grammar_rules_iterator_type,
  typename trie_node_pointer_type
>
void print_trie
(
  grammar_rules_iterator_type const &rules_begin,
  trie_node_pointer_type node,
  int32_t dist
)
{
  static size_t depth {0};
  if (node != nullptr)
  {
    std::cout << depth << ':';
    if (node->edge_begin_dist != node->edge_end_dist)
    {
      auto it {std::next(rules_begin, node->edge_begin_dist)};
      auto end {std::next(rules_begin, node->edge_end_dist)};
      while (it != end)
      {
        if (*it != '\n')
        {
          std::cout << *it;
        }
        else
        {
          std::cout << "\\n";
        }
        it += dist;
      }
    }
    std::cout << "(" << node->leftmost_rank << ',' << node->rightmost_rank << ")\n";
    for (auto character_child_pair : node->branches)
    {
      ++depth;
      print_trie(rules_begin, std::get<1>(character_child_pair), dist);
      --depth;
    }
  }
  return;
}

template
<
  typename trie_node_pointer_type,
  typename grammar_rules_iterator_type
>
void insert_grammar_rule
(
  trie_node_pointer_type node,
  grammar_rules_iterator_type rules_begin,
  grammar_rules_iterator_type rule_it,
  grammar_rules_iterator_type rule_end,
  int32_t dist,
  uint32_t lex_rank
)
{
  while (rule_it != rule_end)
  {
    auto character {*rule_it};
    if (node->branches.find(character) == std::end(node->branches))
    {
      node->branches[character] = std::make_shared<trie_node>
      (
        std::distance(rules_begin, rule_it),
        std::distance(rules_begin, rule_end),
        lex_rank
      );
      return;
    }
    else
    {
      auto child {node->branches[character]};
      auto edge_it {std::next(rules_begin, child->edge_begin_dist)};
      auto edge_end {std::next(rules_begin, child->edge_end_dist)};
      while
      (
        (rule_it != rule_end)
        &&
        (edge_it != edge_end)
        &&
        (*rule_it == *edge_it)
      )
      {
        rule_it += dist;
        edge_it += dist;
      }
      if (edge_it == edge_end)
      {
        if (rule_it != rule_end)
        {
          node = child;
        }
        else
        {
          child->leftmost_rank = lex_rank;
          child->rightmost_rank = lex_rank;
        }
      }
      else
      {
        auto edge_it_dist {std::distance(rules_begin, edge_it)};
        auto internal_node {std::make_shared<trie_node>(child->edge_begin_dist, edge_it_dist, 0)};
        child->edge_begin_dist = edge_it_dist;
        internal_node->branches[*edge_it] = child;
        if (rule_it != rule_end)
        {
          internal_node->branches[*rule_it] = std::make_shared<trie_node>
          (
            std::distance(rules_begin, rule_it),
            std::distance(rules_begin, rule_end),
            lex_rank
          );
        }
        else
        {
          internal_node->leftmost_rank = lex_rank;
          internal_node->rightmost_rank = lex_rank;
        }
        node->branches[character] = internal_node;
        return;
      }
    }
  }
  return;
}

template
<
  typename grammar_rule_sizes_type,
  typename grammar_rules_type,
  typename trie_node_pointer_type
>
void insert_grammar_rules
(
  grammar_rule_sizes_type const &grammar_rule_sizes,
  grammar_rules_type const &grammar_rules,
  trie_node_pointer_type &lex_trie_root,
  trie_node_pointer_type &colex_trie_root
)
{
  lex_trie_root = std::make_shared<trie_node>();
  colex_trie_root = std::make_shared<trie_node>();
  uint32_t lex_rank {1};
  auto sizes_it {std::begin(grammar_rule_sizes)};
  auto rules_begin {std::begin(grammar_rules)};
  auto rules_it {rules_begin};
  auto rules_end {std::end(grammar_rules)};
  while (rules_it != rules_end)
  {
    auto rule_begin {rules_it};
    auto rule_end {std::next(rule_begin, *sizes_it)};
    insert_grammar_rule
    (
      lex_trie_root,
      rules_begin,
      rule_begin,
      rule_end,
      1,
      lex_rank
    );
    insert_grammar_rule
    (
      colex_trie_root,
      rules_begin,
      std::prev(rule_end),
      std::prev(rule_begin),
      -1,
      lex_rank
    );
    ++lex_rank;
    ++sizes_it;
    rules_it = rule_end;
  }
  return;
}

template <typename trie_node_pointer_type>
void calculate_lex_trie_rank_ranges (trie_node_pointer_type node)
{
  if (!node->branches.empty())
  {
    auto branches_begin {std::begin(node->branches)};
    auto branches_it {branches_begin};
    auto branches_end {std::end(node->branches)};
    while (branches_it != branches_end)
    {
      calculate_lex_trie_rank_ranges(std::get<1>(*branches_it));
      ++branches_it;
    }
    auto first_child {std::get<1>(*branches_begin)};
    auto last_child {std::get<1>(*std::prev(branches_end))};
    if (node->leftmost_rank == 0)
    {
      node->leftmost_rank = first_child->leftmost_rank;
    }
    node->rightmost_rank = last_child->rightmost_rank;
  }
  return;
}

template
<
  typename trie_node_pointer_type,
  typename lex_colex_permutation_type
>
void calculate_colex_trie_rank_ranges_and_lex_colex_permutation
(
  trie_node_pointer_type node,
  lex_colex_permutation_type &lex_colex_permutation
)
{
  static uint32_t colex_rank {1};
  if (node->leftmost_rank != 0)
  {
    lex_colex_permutation[node->leftmost_rank] = colex_rank;
    node->leftmost_rank = node->rightmost_rank = colex_rank;
    ++colex_rank;
  }
  if (!node->branches.empty())
  {
    auto branches_begin {std::begin(node->branches)};
    auto branches_it {branches_begin};
    auto branches_end {std::end(node->branches)};
    while (branches_it != branches_end)
    {
      calculate_colex_trie_rank_ranges_and_lex_colex_permutation
      (
        std::get<1>(*branches_it),
        lex_colex_permutation
      );
      ++branches_it;
    }
    auto first_child {std::get<1>(*branches_begin)};
    auto last_child {std::get<1>(*std::prev(branches_end))};
    if (node->leftmost_rank == 0)
    {
      node->leftmost_rank = first_child->leftmost_rank;
    }
    node->rightmost_rank = last_child->rightmost_rank;
  }
  return;
}

template
<
  typename gc_text_type,
  typename temp_gc_text_iterator_type
>
void caculate_gc_text
(
  gc_text_type &gc_text,
  temp_gc_text_iterator_type temp_gc_text_begin,
  temp_gc_text_iterator_type temp_gc_text_end
)
{
  gc_text.resize(std::distance(temp_gc_text_begin, temp_gc_text_end) + 1);
  gc_text[std::size(gc_text) - 1] = 0;
  std::copy(temp_gc_text_begin, temp_gc_text_end, std::begin(gc_text));
  return;
}

template
<
  typename gc_text_type,
  typename temp_sa_bwt_type,
  typename lex_gc_bwt_wt_type
>
void calculate_lex_gc_bwt_wt
(
  gc_text_type const &gc_text,
  temp_sa_bwt_type &temp_sa_bwt,
  lex_gc_bwt_wt_type &lex_gc_bwt_wt
)
{
  sdsl::qsufsort::construct_sa(temp_sa_bwt, gc_text);
  auto temp_sa_bwt_it {std::begin(temp_sa_bwt)};
  auto temp_sa_bwt_end {std::end(temp_sa_bwt)};
  while (temp_sa_bwt_it != temp_sa_bwt_end)
  {
    if (*temp_sa_bwt_it != 0)
    {
      *temp_sa_bwt_it = gc_text[(*temp_sa_bwt_it) - 1];
    }
    ++temp_sa_bwt_it;
  }
  construct_im(lex_gc_bwt_wt, temp_sa_bwt);
  return;
}

template
<
  typename gc_text_type,
  typename temp_sa_bwt_type,
  typename lex_colex_permutation_type,
  typename colex_gc_bwt_wt_type
>
void calculate_colex_gc_bwt_wt
(
  gc_text_type const &gc_text,
  temp_sa_bwt_type &temp_sa_bwt,
  lex_colex_permutation_type const &lex_colex_permutation,
  colex_gc_bwt_wt_type &colex_gc_bwt_wt
)
{
  sdsl::qsufsort::construct_sa(temp_sa_bwt, gc_text);
  auto temp_sa_bwt_it {std::begin(temp_sa_bwt)};
  auto temp_sa_bwt_end {std::end(temp_sa_bwt)};
  while (temp_sa_bwt_it != temp_sa_bwt_end)
  {
    if (*temp_sa_bwt_it != 0)
    {
      *temp_sa_bwt_it = lex_colex_permutation[gc_text[(*temp_sa_bwt_it) - 1]];
    }
    ++temp_sa_bwt_it;
  }
  construct_im(colex_gc_bwt_wt, temp_sa_bwt);
  return;
}

template <typename pattern_iterator_type>
bool exists_another_factorization
(
  pattern_iterator_type rit,
  pattern_iterator_type rend
)
{
  while
  (
    (std::next(rend) != rit)
    &&
    (*std::prev(rit) == *rit)
  )
  {
    --rit;
  }
  if (std::next(rend) != rit)
  {
    return (*std::prev(rit) > *rit);
  }
  return false;
}

struct gc_index
{
  using text_type = sdsl::int_vector<8>;

  sdsl::int_vector<> grammar_rule_sizes;
  sdsl::int_vector<8> grammar_rules;
  std::shared_ptr<trie_node> lex_trie_root;
  std::shared_ptr<trie_node> colex_trie_root;
  sdsl::int_vector<> lex_gc_character_bucket_end_dists;
  sdsl::wt_int<> lex_gc_bwt_wt;
  sdsl::wt_int<> colex_gc_bwt_wt;

  gc_index () = default;

  gc_index (text_type &text)
  {
    sdsl::bit_vector sl_types;
    calculate_sl_types(text, sl_types);

    sdsl::int_vector<> text_dists(std::size(text));
    auto invalid_text_dist {std::size(text)};
    std::fill(std::begin(text_dists), std::end(text_dists), invalid_text_dist);

    sdsl::int_vector<> character_bucket_dists(256);
    calculate_character_bucket_begin_dists(text, character_bucket_dists);
    bucket_sort_rightmost_l_type_characters(text, sl_types, character_bucket_dists, text_dists);
    induce_sort_l_type_characters(text, sl_types, character_bucket_dists, text_dists);
    calculate_character_bucket_end_dists(text, character_bucket_dists);
    induce_sort_s_type_characters(text, sl_types, character_bucket_dists, text_dists);

    auto text_dists_boundary {collect_valid_entries(std::begin(text_dists), std::end(text_dists), invalid_text_dist)};
    auto grammar_rule_begin_dists_begin {std::begin(text_dists)};
    auto grammar_rule_begin_dists_end {text_dists_boundary};
    auto temp_gc_text_begin {text_dists_boundary};
    auto temp_gc_text_end {std::end(text_dists)};

    calculate_grammar_rule_begin_dists_and_temp_gc_text
    (
      text,
      sl_types,
      grammar_rule_begin_dists_begin,
      grammar_rule_begin_dists_end,
      temp_gc_text_begin,
      temp_gc_text_end
    );

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
    insert_grammar_rules
    (
      grammar_rule_sizes,
      grammar_rules,
      lex_trie_root,
      colex_trie_root
    );
    calculate_lex_trie_rank_ranges(lex_trie_root);

    sdsl::int_vector<> lex_colex_permutation(std::size(grammar_rule_sizes) + 1);
    lex_colex_permutation[0] = 0;
    calculate_colex_trie_rank_ranges_and_lex_colex_permutation
    (
      colex_trie_root,
      lex_colex_permutation
    );

    sdsl::int_vector<> gc_text;
    caculate_gc_text(gc_text, temp_gc_text_begin, temp_gc_text_end);

    lex_gc_character_bucket_end_dists.resize(std::size(grammar_rule_sizes) + 1);
    calculate_character_bucket_end_dists(gc_text, lex_gc_character_bucket_end_dists);

    text_dists.resize(std::size(gc_text));
    auto &temp_sa_bwt {text_dists};
    calculate_lex_gc_bwt_wt(gc_text, temp_sa_bwt, lex_gc_bwt_wt);
    calculate_colex_gc_bwt_wt(gc_text, temp_sa_bwt, lex_colex_permutation, colex_gc_bwt_wt);
  }

  template <typename pattern_iterator_type>
  void calculate_rlast
  (
    pattern_iterator_type &rlast,
    pattern_iterator_type const &rend,
    uint8_t prev_sl_type
  )
  {
    while
    (
      (rlast != rend)
      &&
      ! (
          (prev_sl_type == S)
          &&
          (*rlast > *std::next(rlast))
        )
    )
    {
      if
      (
        (prev_sl_type == L)
        &&
        (*rlast < *std::next(rlast))
      )
      {
        prev_sl_type = S;
      }
      --rlast;
    }
    return;
  }

  template
  <
    typename trie_node_pointer_type,
    typename string_iterator_type
  >
  void search_grammar_rule
  (
    trie_node_pointer_type node,
    string_iterator_type it,
    string_iterator_type end,
    int32_t dist,
    uint32_t &leftmost_rank,
    uint32_t &rightmost_rank
  )
  {
    while (node->branches.find(*it) != std::end(node->branches))
    {
      node = node->branches[*it];
      auto edge_it {std::next(std::begin(grammar_rules), node->edge_begin_dist)};
      auto edge_end {std::next(std::begin(grammar_rules), node->edge_end_dist)};
      while
      (
        (it != end)
        &&
        (edge_it != edge_end)
        &&
        (*it == *edge_it)
      )
      {
        it += dist;
        edge_it += dist;
      }
      if (it != end)
      {
        if (edge_it != edge_end)
        {
          return;
        }
      }
      else
      {
        rightmost_rank = node->rightmost_rank;
        leftmost_rank = node->leftmost_rank;
        return;
      }
    }
    return;
  }

  template <typename pattern_iterator_type>
  void backward_search_suffix
  (
    uint32_t &begin_dist,
    uint32_t &end_dist,
    pattern_iterator_type begin,
    pattern_iterator_type end
  )
  {
    uint32_t leftmost_lex_rank {0};
    uint32_t rightmost_lex_rank {0};
    search_grammar_rule
    (
      lex_trie_root,
      begin,
      end,
      1,
      leftmost_lex_rank,
      rightmost_lex_rank
    );
    if (leftmost_lex_rank != 0)
    {
      begin_dist = lex_gc_character_bucket_end_dists[leftmost_lex_rank - 1];
      end_dist = lex_gc_character_bucket_end_dists[rightmost_lex_rank];
    }
    return;
  }

  template <typename pattern_iterator_type>
  void backward_search_pivot
  (
    uint32_t &begin_dist,
    uint32_t &end_dist,
    pattern_iterator_type begin,
    pattern_iterator_type end
  )
  {
    uint32_t lex_rank {0};
    search_grammar_rule
    (
      lex_trie_root,
      begin,
      end,
      1,
      lex_rank,
      lex_rank
    );
    if (lex_rank != 0)
    {
      auto character_begin_dist {lex_gc_character_bucket_end_dists[lex_rank - 1]};
      begin_dist = character_begin_dist + lex_gc_bwt_wt.rank(begin_dist, lex_rank);
      end_dist = character_begin_dist + lex_gc_bwt_wt.rank(end_dist, lex_rank);
    }
    else
    {
      begin_dist = end_dist;
    }
    return;
  }

  template <typename pattern_iterator_type>
  uint32_t backward_search_prefix
  (
    uint32_t begin_dist,
    uint32_t end_dist,
    pattern_iterator_type rbegin,
    pattern_iterator_type rend
  )
  {
    uint32_t leftmost_colex_rank {0};
    uint32_t rightmost_colex_rank {0};
    search_grammar_rule
    (
      colex_trie_root,
      rbegin,
      rend,
      -1,
      leftmost_colex_rank,
      rightmost_colex_rank
    );
    if (leftmost_colex_rank != 0)
    {
      return std::get<0>
      (
        colex_gc_bwt_wt.range_search_2d
        (
          begin_dist,
          end_dist - 1,
          leftmost_colex_rank,
          rightmost_colex_rank
        )
      );
    }
    return 0;
  }

  template <typename pattern_iterator_type>
  auto count
  (
    pattern_iterator_type rbegin,
    pattern_iterator_type rend,
    uint8_t prev_sl_type
  )
  {
    uint32_t begin_dist {0};
    uint32_t end_dist {0};
    auto rfirst {rbegin};
    auto rlast {std::prev(rfirst)};
    calculate_rlast(rlast, rend, prev_sl_type);
    backward_search_suffix(begin_dist, end_dist, std::next(rlast), std::next(rfirst));
    while (begin_dist != end_dist)
    {
      rfirst = rlast--;
      calculate_rlast(rlast, rend, L);
      if (rlast != rend)
      {
        backward_search_pivot(begin_dist, end_dist, std::next(rlast), std::next(rfirst));
      }
      else
      {
        return backward_search_prefix(begin_dist, end_dist, rfirst, rend);
      }
    }
    return (end_dist - begin_dist);
  }

  template <typename pattern_iterator_type>
  uint32_t count
  (
    pattern_iterator_type pattern_begin,
    pattern_iterator_type pattern_end
  )
  {
    uint32_t count_ {0};
    if (std::distance(pattern_begin, pattern_end) != 0)
    {
      count_ = count(std::prev(pattern_end), std::prev(pattern_begin), S);
      if (exists_another_factorization(std::prev(pattern_end), std::prev(pattern_begin)))
      {
        count_ += count(std::prev(pattern_end), std::prev(pattern_begin), L);
      }
    }
    return count_;
  }

};

void construct
(
  gc_index &index,
  char const *input_text_file
)
{
  sdsl::int_vector<8> text;
  sdsl::load_vector_from_file(text, input_text_file);
  sdsl::append_zero_symbol(text);

  sdsl::bit_vector sl_types;
  calculate_sl_types(text, sl_types);

  auto invalid_text_dist {std::size(text)};
  sdsl::int_vector<> text_dists(std::size(text), invalid_text_dist, sdsl::bits::hi(std::size(text)) + 1);

  sdsl::int_vector<> character_bucket_dists(256, 0, sdsl::bits::hi(std::size(text)) + 1);
  calculate_character_bucket_begin_dists(text, character_bucket_dists);
  bucket_sort_rightmost_l_type_characters(text, sl_types, character_bucket_dists, text_dists);
  induce_sort_l_type_characters(text, sl_types, character_bucket_dists, text_dists);
  calculate_character_bucket_end_dists(text, character_bucket_dists);
  induce_sort_s_type_characters(text, sl_types, character_bucket_dists, text_dists);

  auto text_dists_boundary {collect_valid_entries(std::begin(text_dists), std::end(text_dists), invalid_text_dist)};
  auto grammar_rule_begin_dists_begin {std::begin(text_dists)};
  auto grammar_rule_begin_dists_end {text_dists_boundary};
  auto temp_gc_text_begin {text_dists_boundary};
  auto temp_gc_text_end {std::end(text_dists)};

  calculate_grammar_rule_begin_dists_and_temp_gc_text
  (
    text,
    sl_types,
    grammar_rule_begin_dists_begin,
    grammar_rule_begin_dists_end,
    temp_gc_text_begin,
    temp_gc_text_end
  );

  // calculate_grammar_rule_sizes
  // (
  //   grammar_rule_sizes,
  //   sl_types,
  //   grammar_rule_begin_dists_begin,
  //   grammar_rule_begin_dists_end
  // );
  //
  // calculate_grammar_rules
  // (
  //   text,
  //   grammar_rule_sizes,
  //   grammar_rules,
  //   grammar_rule_begin_dists_begin
  // );
  // insert_grammar_rules
  // (
  //   grammar_rule_sizes,
  //   grammar_rules,
  //   lex_trie_root,
  //   colex_trie_root
  // );
  // calculate_lex_trie_rank_ranges(lex_trie_root);
  //
  // sdsl::int_vector<> lex_colex_permutation(std::size(grammar_rule_sizes) + 1);
  // lex_colex_permutation[0] = 0;
  // calculate_colex_trie_rank_ranges_and_lex_colex_permutation
  // (
  //   colex_trie_root,
  //   lex_colex_permutation
  // );
  //
  // sdsl::int_vector<> gc_text;
  // caculate_gc_text(gc_text, temp_gc_text_begin, temp_gc_text_end);
  //
  // lex_gc_character_bucket_end_dists.resize(std::size(grammar_rule_sizes) + 1);
  // calculate_character_bucket_end_dists(gc_text, lex_gc_character_bucket_end_dists);
  //
  // text_dists.resize(std::size(gc_text));
  // auto &temp_sa_bwt {text_dists};
  // calculate_lex_gc_bwt_wt(gc_text, temp_sa_bwt, lex_gc_bwt_wt);
  // calculate_colex_gc_bwt_wt(gc_text, temp_sa_bwt, lex_colex_permutation, colex_gc_bwt_wt);
  show_type(index);
  return;
}

}

#endif
