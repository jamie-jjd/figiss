#ifndef GRAMMAR_COMPRESSED_INDEX_HPP_
#define GRAMMAR_COMPRESSED_INDEX_HPP_

#include <map>

#include <sdsl/suffix_trees.hpp>
#include <tudocomp_stat/StatPhase.hpp>

#include "utility.hpp"

// #define IS_LOGGED

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
  sdsl::util::set_to_value(sl_types, S);
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
  sdsl::util::_set_zero_bits(character_bucket_dists);
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
  auto it {first};
  auto new_last {first};
  while (it != last)
  {
    if (*it != invalid_value)
    {
      *new_last = *it;
      if (new_last != it)
      {
        *it = invalid_value;
      }
      ++new_last;
    }
    ++it;
  }
  return new_last;
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
    auto sl_types_last {std::next(sl_types_first)};
    while (!is_leftmost_s(sl_types_last))
    {
      ++sl_types_last;
    }
    *grammar_rule_sizes_it = std::distance(sl_types_first, sl_types_last);
    ++grammar_rule_sizes_it;
    ++begin_dists_it;
  }
  return;
}

struct trie_node
{
  int64_t edge_begin_dist;
  int64_t edge_end_dist;
  std::map<uint8_t, trie_node*> branches;
  uint64_t leftmost_rank;
  uint64_t rightmost_rank;

  trie_node () = default;

  trie_node
  (
    int64_t edge_begin_dist_,
    int64_t edge_end_dist_,
    uint64_t rank
  )
  : edge_begin_dist {edge_begin_dist_},
    edge_end_dist {edge_end_dist_},
    leftmost_rank {rank},
    rightmost_rank {rank}
  {
  }

  ~trie_node ()
  {
    auto branches_it {std::begin(branches)};
    auto branches_end {std::end(branches)};
    while (branches_it != branches_end)
    {
      delete std::get<1>(*branches_it);
      ++branches_it;
    }
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

template <typename grammar_rules_iterator_type>
void print_trie
(
  grammar_rules_iterator_type rules_begin,
  trie_node *node,
  int64_t const dist,
  uint64_t depth = 0
)
{
  if (node != nullptr)
  {
    std::cout << depth << ':';
    if (node->edge_begin_dist != node->edge_end_dist)
    {
      auto it {std::next(rules_begin, node->edge_begin_dist)};
      auto end {std::next(rules_begin, node->edge_end_dist)};
      while (it != end)
      {
        std::cout << *it;
        it += dist;
      }
    }
    std::cout << "(" << node->leftmost_rank << ',' << node->rightmost_rank << ")\n";
    auto branches_it {std::begin(node->branches)};
    auto branches_end {std::end(node->branches)};
    while (branches_it != branches_end)
    {
      print_trie(rules_begin, std::get<1>(*branches_it), dist, depth + 1);
      ++branches_it;
    }
  }
  return;
}

template <typename grammar_rules_iterator_type>
void insert_grammar_rule
(
  trie_node *node,
  grammar_rules_iterator_type rules_begin,
  grammar_rules_iterator_type rule_it,
  grammar_rules_iterator_type rule_end,
  int64_t const dist,
  uint64_t const lex_rank
)
{
  while (rule_it != rule_end)
  {
    auto character {*rule_it};
    if (node->branches.find(character) == std::end(node->branches))
    {
      node->branches[character] = new trie_node
      {
        std::distance(rules_begin, rule_it),
        std::distance(rules_begin, rule_end),
        lex_rank
      };
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
        auto internal_node {new trie_node {child->edge_begin_dist, edge_it_dist, 0}};
        child->edge_begin_dist = edge_it_dist;
        internal_node->branches[*edge_it] = child;
        if (rule_it != rule_end)
        {
          internal_node->branches[*rule_it] = new trie_node
          {
            std::distance(rules_begin, rule_it),
            std::distance(rules_begin, rule_end),
            lex_rank
          };
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
  trie_node_pointer_type &colex_trie_root,
  uint64_t const rank_inc
)
{
  lex_trie_root = new trie_node {};
  colex_trie_root = new trie_node {};
  uint64_t rank {1};
  auto sizes_it {std::begin(grammar_rule_sizes)};
  auto rules_it {std::begin(grammar_rules)};
  while (rules_it != std::end(grammar_rules))
  {
    auto rule_begin {rules_it};
    auto rule_end {std::next(rule_begin, *sizes_it)};
    insert_grammar_rule
    (
      lex_trie_root,
      std::begin(grammar_rules),
      rule_begin,
      rule_end,
      1,
      rank
    );
    insert_grammar_rule
    (
      colex_trie_root,
      std::begin(grammar_rules),
      std::prev(rule_end),
      std::prev(rule_begin),
      -1,
      rank
    );
    rank += rank_inc;
    ++sizes_it;
    rules_it = rule_end;
  }
  return;
}

void calculate_lex_trie_rank_ranges (trie_node *node)
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

template <typename lex_colex_permutation_type>
void calculate_colex_trie_rank_ranges_and_lex_colex_permutation
(
  trie_node *node,
  lex_colex_permutation_type &lex_colex_permutation,
  uint64_t &colex_rank
)
{
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
        lex_colex_permutation,
        colex_rank
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
  typename lex_colex_permutation_type,
  typename colex_lex_permutation_type
>
void calculate_colex_lex_permutation
(
  lex_colex_permutation_type const &lex_colex_permutation,
  colex_lex_permutation_type &colex_lex_permutation
)
{
  for (uint64_t lex_rank {0}; lex_rank != std::size(lex_colex_permutation); ++lex_rank)
  {
    colex_lex_permutation[lex_colex_permutation[lex_rank]] = lex_rank;
  }
  return;
}

template
<
  typename gc_text_type,
  typename temp_sa_bwt_type,
  typename lex_colex_permutation_type,
  typename gc_bwt_type
>
void calculate_colex_gc_bwt
(
  gc_text_type const &gc_text,
  temp_sa_bwt_type &temp_sa_bwt,
  lex_colex_permutation_type const &lex_colex_permutation,
  gc_bwt_type &colex_gc_bwt
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
  sdsl::construct_im(colex_gc_bwt, temp_sa_bwt);
  return;
}

void calculate_trie_rank_ranges
(
  trie_node *node,
  uint64_t &rank
)
{
  if (node->leftmost_rank != 0)
  {
    node->leftmost_rank = node->rightmost_rank = rank++;
  }
  if (!node->branches.empty())
  {
    auto branches_begin {std::begin(node->branches)};
    auto branches_it {branches_begin};
    auto branches_end {std::end(node->branches)};
    while (branches_it != branches_end)
    {
      calculate_trie_rank_ranges(std::get<1>(*branches_it), rank);
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

struct gc_index
{
  sdsl::int_vector<> grammar_rule_sizes;
  sdsl::int_vector<8> grammar_rules;
  trie_node *lex_trie_root;
  trie_node *colex_trie_root;
  sdsl::int_vector<> colex_lex_permutation;
  sdsl::int_vector<> lex_gc_character_bucket_end_dists;
  sdsl::wt_int<> colex_gc_bwt;

  ~gc_index ()
  {
    delete lex_trie_root;
    delete colex_trie_root;
  }
};

void print_gc_index (gc_index &index)
{
  std::cout << index.grammar_rule_sizes << '\n';
  std::cout << index.grammar_rules << '\n';
  print_trie(std::begin(index.grammar_rules), index.lex_trie_root, 1);
  print_trie(std::begin(index.grammar_rules), index.colex_trie_root, -1);
  std::cout << index.colex_lex_permutation << '\n';
  std::cout << index.lex_gc_character_bucket_end_dists << '\n';
  std::cout << index.colex_gc_bwt << '\n';
  return;
}

void construct
(
  gc_index &index,
  std::string const text_path
)
{
  std::ofstream json_output {"../output/construct/" + util::basename(text_path) + ".gci.json"};
  tdc::StatPhase phases {"construct_gc_index"};
  {
    sdsl::int_vector<8> text;
    tdc::StatPhase::wrap
    (
      "load_text",
      [&] ()
      {
        sdsl::load_vector_from_file(text, text_path);
        sdsl::append_zero_symbol(text);
      }
    );
    sdsl::bit_vector sl_types;
    tdc::StatPhase::wrap
    (
      "calculate_sl_types",
      [&] ()
      {
        calculate_sl_types(text, sl_types);
      }
    );
    auto invalid_text_dist {std::size(text)};
    auto text_size_width {sdsl::bits::hi(std::size(text)) + 1};
    sdsl::int_vector<> text_dists;
    tdc::StatPhase::wrap
    (
      "init_text_dists",
      [&] ()
      {
        text_dists.width(text_size_width);
        text_dists.resize(std::size(text));
        sdsl::util::set_to_value(text_dists, invalid_text_dist);
      }
    );
    sdsl::int_vector<> character_bucket_dists;
    tdc::StatPhase::wrap
    (
      "init_character_bucket_dists",
      [&] ()
      {
        character_bucket_dists.width(text_size_width);
        character_bucket_dists.resize(256);
      }
    );
    tdc::StatPhase::wrap
    (
      "calculate_character_bucket_begin_dists",
      [&] ()
      {
        calculate_character_bucket_begin_dists(text, character_bucket_dists);
      }
    );
    tdc::StatPhase::wrap
    (
      "bucket_sort_rightmost_l_type_characters",
      [&] ()
      {
        bucket_sort_rightmost_l_type_characters(text, sl_types, character_bucket_dists, text_dists);
      }
    );
    tdc::StatPhase::wrap
    (
      "induce_sort_l_type_characters",
      [&] ()
      {
        induce_sort_l_type_characters(text, sl_types, character_bucket_dists, text_dists);
      }
    );
    tdc::StatPhase::wrap
    (
      "calculate_character_bucket_end_dists",
      [&] ()
      {
        calculate_character_bucket_end_dists(text, character_bucket_dists);
      }
    );
    tdc::StatPhase::wrap
    (
      "induce_sort_s_type_characters",
      [&] ()
      {
        induce_sort_s_type_characters(text, sl_types, character_bucket_dists, text_dists);
      }
    );
    auto text_dists_boundary {std::begin(text_dists)};
    tdc::StatPhase::wrap
    (
      "calculate_text_dists_boundary",
      [&] ()
      {
        text_dists_boundary = collect_valid_entries
        (
          std::begin(text_dists),
          std::end(text_dists),
          invalid_text_dist
        );
      }
    );
    auto grammar_rule_begin_dists_begin {std::begin(text_dists)};
    auto grammar_rule_begin_dists_end {text_dists_boundary};
    auto temp_gc_text_begin {text_dists_boundary};
    auto temp_gc_text_end {std::end(text_dists)};
    tdc::StatPhase::wrap
    (
      "calculate_grammar_rule_begin_dists_and_temp_gc_text",
      [&] ()
      {
        calculate_grammar_rule_begin_dists_and_temp_gc_text
        (
          text,
          sl_types,
          grammar_rule_begin_dists_begin,
          grammar_rule_begin_dists_end,
          temp_gc_text_begin,
          temp_gc_text_end
        );
      }
    );
    tdc::StatPhase::wrap
    (
      "calculate_grammar_rule_sizes",
      [&] ()
      {
        calculate_grammar_rule_sizes
        (
          index.grammar_rule_sizes,
          sl_types,
          grammar_rule_begin_dists_begin,
          grammar_rule_begin_dists_end
        );
      }
    );
    tdc::StatPhase::wrap
    (
      "calculate_grammar_rules",
      [&] ()
      {
        calculate_grammar_rules
        (
          text,
          index.grammar_rule_sizes,
          index.grammar_rules,
          grammar_rule_begin_dists_begin
        );
      }
    );
    tdc::StatPhase::wrap
    (
      "clear_sl_types_and_text",
      [&] ()
      {
        sdsl::util::clear(sl_types);
        sdsl::util::clear(text);
      }
    );
    tdc::StatPhase::wrap
    (
      "insert_grammar_rules",
      [&] ()
      {
        insert_grammar_rules
        (
          index.grammar_rule_sizes,
          index.grammar_rules,
          index.lex_trie_root,
          index.colex_trie_root,
          1
        );
      }
    );
    tdc::StatPhase::wrap
    (
      "calculate_lex_trie_rank_ranges",
      [&] ()
      {
        calculate_lex_trie_rank_ranges(index.lex_trie_root);
      }
    );
    auto gc_text_sigma {std::size(index.grammar_rule_sizes) + 1};
    auto gc_text_width {sdsl::bits::hi(gc_text_sigma) + 1};
    sdsl::int_vector<> lex_colex_permutation;
    tdc::StatPhase::wrap
    (
      "init_lex_colex_permutation",
      [&] ()
      {
        lex_colex_permutation.width(gc_text_width);
        lex_colex_permutation.resize(gc_text_sigma);
        lex_colex_permutation[0] = 0;
      }
    );
    tdc::StatPhase::wrap
    (
      "calculate_colex_trie_rank_ranges_and_lex_colex_permutation",
      [&] ()
      {
        uint64_t colex_rank {1};
        calculate_colex_trie_rank_ranges_and_lex_colex_permutation
        (
          index.colex_trie_root,
          lex_colex_permutation,
          colex_rank
        );
      }
    );
    tdc::StatPhase::wrap
    (
      "calculate_colex_lex_permutation",
      [&] ()
      {
        index.colex_lex_permutation.width(gc_text_width);
        index.colex_lex_permutation.resize(gc_text_sigma);
        calculate_colex_lex_permutation
        (
          lex_colex_permutation,
          index.colex_lex_permutation
        );
      }
    );
    auto gc_text_size {std::distance(temp_gc_text_begin, temp_gc_text_end) + 1};
    auto gc_text_size_width {sdsl::bits::hi(gc_text_size) + 1};
    sdsl::int_vector<> gc_text;
    tdc::StatPhase::wrap
    (
      "calculate_gc_text",
      [&] ()
      {
        gc_text.width(gc_text_width);
        gc_text.resize(gc_text_size);
        std::copy(temp_gc_text_begin, temp_gc_text_end, std::begin(gc_text));
        gc_text[std::size(gc_text) - 1] = 0;
      }
    );
    tdc::StatPhase::wrap
    (
      "calculate_lex_gc_character_bucket_end_dists",
      [&] ()
      {
        index.lex_gc_character_bucket_end_dists.width(gc_text_size_width);
        index.lex_gc_character_bucket_end_dists.resize(gc_text_sigma);
        calculate_character_bucket_end_dists(gc_text, index.lex_gc_character_bucket_end_dists);
      }
    );
    tdc::StatPhase::wrap
    (
      "calculate_colex_gc_bwt",
      [&] ()
      {
        text_dists.resize(std::size(gc_text));
        auto &temp_sa_bwt {text_dists};
        calculate_colex_gc_bwt
        (
          gc_text,
          temp_sa_bwt,
          lex_colex_permutation,
          index.colex_gc_bwt
        );
      }
    );
  }
  phases.to_json().str(json_output);
  return;
}

void serialize
(
  gc_index &index,
  std::string const index_path
)
{
  std::ofstream index_output {index_path};
  std::ofstream json_output {"../output/serialize/" + util::basename(index_path) + ".json"};
  tdc::StatPhase phases {"serialize_gc_index"};
  tdc::StatPhase::wrap
  (
    "serialize_grammar_rule_sizes",
    [&] ()
    {
      index.grammar_rule_sizes.serialize(index_output);
    }
  );
  tdc::StatPhase::wrap
  (
    "serialize_grammar_rules",
    [&] ()
    {
      index.grammar_rules.serialize(index_output);
    }
  );
  tdc::StatPhase::wrap
  (
    "serialize_colex_lex_permutation",
    [&] ()
    {
      index.colex_lex_permutation.serialize(index_output);
    }
  );
  tdc::StatPhase::wrap
  (
    "serialize_lex_gc_character_bucket_end_dists",
    [&] ()
    {
      index.lex_gc_character_bucket_end_dists.serialize(index_output);
    }
  );
  tdc::StatPhase::wrap
  (
    "serialize_colex_gc_bwt",
    [&] ()
    {
      index.colex_gc_bwt.serialize(index_output);
    }
  );
  phases.to_json().str(json_output);
  return;
}

void load
(
  gc_index &index,
  std::string const index_path
)
{
  std::ifstream index_input {index_path};
  std::ofstream json_output {"../output/load/" + util::basename(index_path) + ".json"};
  tdc::StatPhase phases {"load_gc_index"};
  {
    tdc::StatPhase::wrap
    (
      "load_grammar_rule_sizes",
      [&] ()
      {
        sdsl::util::clear(index.grammar_rule_sizes);
        index.grammar_rule_sizes.load(index_input);
      }
    );
    tdc::StatPhase::wrap
    (
      "load_grammar_rules",
      [&] ()
      {
        sdsl::util::clear(index.grammar_rules);
        index.grammar_rules.load(index_input);
      }
    );
    tdc::StatPhase::wrap
    (
      "load_colex_lex_permutation",
      [&] ()
      {
        sdsl::util::clear(index.colex_lex_permutation);
        index.colex_lex_permutation.load(index_input);
      }
    );
    tdc::StatPhase::wrap
    (
      "load_lex_gc_character_bucket_end_dists",
      [&] ()
      {
        sdsl::util::clear(index.lex_gc_character_bucket_end_dists);
        index.lex_gc_character_bucket_end_dists.load(index_input);
      }
    );
    tdc::StatPhase::wrap
    (
      "load_colex_gc_bwt",
      [&] ()
      {
        sdsl::util::clear(index.colex_gc_bwt);
        index.colex_gc_bwt.load(index_input);
      }
    );
    tdc::StatPhase::wrap
    (
      "insert_grammar_rules",
      [&] ()
      {
        insert_grammar_rules
        (
          index.grammar_rule_sizes,
          index.grammar_rules,
          index.lex_trie_root,
          index.colex_trie_root,
          0
        );
      }
    );
    tdc::StatPhase::wrap
    (
      "calculate_lex_trie_rank_ranges",
      [&] ()
      {
        uint64_t rank {1};
        calculate_trie_rank_ranges(index.lex_trie_root, rank);
      }
    );
    tdc::StatPhase::wrap
    (
      "calculate_colex_trie_rank_ranges",
      [&] ()
      {
        uint64_t rank {1};
        calculate_trie_rank_ranges(index.colex_trie_root, rank);
      }
    );
  }
  phases.to_json().str(json_output);
  return;
}

template <typename string_iterator_type>
void calculate_sl_factor
(
  string_iterator_type &rfirst,
  string_iterator_type &rlast,
  string_iterator_type rend
)
{
  uint8_t prev_sl_type {L};
  rfirst = rlast--;
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
  typename grammar_rules_iterator_type,
  typename sl_factor_iterator_type
>
void lookup_grammar_rule
(
  grammar_rules_iterator_type grammar_rules_begin,
  trie_node *node,
  sl_factor_iterator_type it,
  sl_factor_iterator_type end,
  int64_t const dist,
  uint64_t &leftmost_rank,
  uint64_t &rightmost_rank
)
{
  leftmost_rank = rightmost_rank = 0;
  while (node->branches.find(*it) != std::end(node->branches))
  {
    node = node->branches[*it];
    auto edge_it {std::next(grammar_rules_begin, node->edge_begin_dist)};
    auto edge_end {std::next(grammar_rules_begin, node->edge_end_dist)};
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
uint64_t backward_search_pattern_prefix
(
  gc_index const &index,
  uint64_t const begin_dist_L,
  uint64_t const end_dist_L,
  uint64_t const begin_dist_S,
  uint64_t const end_dist_S,
  pattern_iterator_type rfirst,
  pattern_iterator_type rlast
)
{
  uint64_t count {0};
  uint64_t leftmost_colex_rank {0};
  uint64_t rightmost_colex_rank {0};
  lookup_grammar_rule
  (
    std::begin(index.grammar_rules),
    index.colex_trie_root,
    rfirst,
    rlast,
    -1,
    leftmost_colex_rank,
    rightmost_colex_rank
  );
  if (leftmost_colex_rank != 0)
  {
#ifdef IS_LOGGED
    for (auto it {std::next(rlast)}; it != std::next(rfirst); ++it)
    {
      std::cout << *it;
    }
    std::cout << " -> (" << leftmost_colex_rank << ", " << rightmost_colex_rank << ")";
#endif
    if (begin_dist_L != end_dist_L)
    {
      uint64_t count_L
      {
        std::get<0>
        (
          index.colex_gc_bwt.range_search_2d
          (
            begin_dist_L,
            (end_dist_L - 1),
            leftmost_colex_rank,
            rightmost_colex_rank,
            false
          )
        )
      };
      count += count_L;
#ifdef IS_LOGGED
      std::cout << " -> L:(" << count_L << ")";
#endif
    }
    if (begin_dist_S != end_dist_S)
    {
      uint64_t count_S
      {
        std::get<0>
        (
          index.colex_gc_bwt.range_search_2d
          (
            begin_dist_S,
            (end_dist_S - 1),
            leftmost_colex_rank,
            rightmost_colex_rank,
            false
          )
        )
      };
      count += count_S;
#ifdef IS_LOGGED
      std::cout << " -> S:(" << count_S << ")";
#endif
    }
#ifdef IS_LOGGED
    std::cout << '\n';
#endif
  }
  return count;
}

template <typename pattern_iterator_type>
void backward_search_exact_sl_factor
(
  gc_index const &index,
  uint64_t &begin_dist_L,
  uint64_t &end_dist_L,
  uint64_t &begin_dist_S,
  uint64_t &end_dist_S,
  pattern_iterator_type rfirst,
  pattern_iterator_type rlast
)
{
  uint64_t colex_rank {0};
  lookup_grammar_rule
  (
    std::begin(index.grammar_rules),
    index.colex_trie_root,
    rfirst,
    rlast,
    -1,
    colex_rank,
    colex_rank
  );
  if (colex_rank != 0)
  {
    auto character_begin_dist
    {
      index.lex_gc_character_bucket_end_dists
      [
        index.colex_lex_permutation[colex_rank] - 1
      ]
    };
#ifdef IS_LOGGED
    for (auto it {std::next(rlast)}; it != std::next(rfirst); ++it)
    {
      std::cout << *it;
    }
    std::cout << " -> (" << colex_rank << ": " << index.colex_lex_permutation[colex_rank] << ")";
#endif
    if (begin_dist_L != end_dist_L)
    {
      begin_dist_L = character_begin_dist + index.colex_gc_bwt.rank(begin_dist_L, colex_rank);
      end_dist_L = character_begin_dist + index.colex_gc_bwt.rank(end_dist_L, colex_rank);
#ifdef IS_LOGGED
      std::cout << " -> L:(" << begin_dist_L << ", " << end_dist_L << ")";
#endif
    }
    if (begin_dist_S != end_dist_S)
    {
      begin_dist_S = character_begin_dist + index.colex_gc_bwt.rank(begin_dist_S, colex_rank);
      end_dist_S = character_begin_dist + index.colex_gc_bwt.rank(end_dist_S, colex_rank);
#ifdef IS_LOGGED
      std::cout << " -> S:(" << begin_dist_S << ", " << end_dist_S << ")";
#endif
    }
#ifdef IS_LOGGED
    std::cout << '\n';
#endif
  }
  else
  {
    begin_dist_L = end_dist_L;
    begin_dist_S = end_dist_S;
  }
  return;
}

template <typename pattern_iterator_type>
auto calculate_pattern_suffix_S_rlast
(
  pattern_iterator_type rit,
  pattern_iterator_type rend
)
{
  while
  (
    (std::prev(rit) != rend)
    &&
    (*std::prev(rit) <= *rit)
  )
  {
    --rit;
  }
  return std::prev(rit);
}

template <typename pattern_iterator_type>
auto backward_search_pattern_suffix
(
  gc_index const &index,
  uint64_t &begin_dist_L,
  uint64_t &end_dist_L,
  uint64_t &begin_dist_S,
  uint64_t &end_dist_S,
  pattern_iterator_type rbegin,
  pattern_iterator_type rend
)
{
  auto rfirst {rbegin};
  auto rlast {calculate_pattern_suffix_S_rlast(rbegin, rend)};
  uint64_t leftmost_lex_rank {0};
  uint64_t rightmost_lex_rank {0};
  if (rlast != rend)
  {
    lookup_grammar_rule
    (
      std::begin(index.grammar_rules),
      index.lex_trie_root,
      std::next(rlast),
      std::next(rfirst),
      1,
      leftmost_lex_rank,
      rightmost_lex_rank
    );
    if (leftmost_lex_rank != 0)
    {
      begin_dist_S = index.lex_gc_character_bucket_end_dists[leftmost_lex_rank - 1];
      end_dist_S = index.lex_gc_character_bucket_end_dists[rightmost_lex_rank];
#ifdef IS_LOGGED
      for (auto it {std::next(rlast)}; it != std::next(rfirst); ++it)
      {
        std::cout << *it;
      }
      std::cout << " -> (" << leftmost_lex_rank << ", " << rightmost_lex_rank << ")";
      std::cout << " -> (" << begin_dist_S << ", " << end_dist_S << ")\n";
#endif
    }
    calculate_sl_factor(rfirst, rlast, rend);
    if (begin_dist_S != end_dist_S)
    {
      uint64_t colex_rank {0};
      lookup_grammar_rule
      (
        std::begin(index.grammar_rules),
        index.colex_trie_root,
        rfirst,
        rlast,
        -1,
        colex_rank,
        colex_rank
      );
      if (colex_rank != 0)
      {
        auto character_begin_dist
        {
          index.lex_gc_character_bucket_end_dists
          [
            index.colex_lex_permutation[colex_rank] - 1
          ]
        };
        begin_dist_S = character_begin_dist + index.colex_gc_bwt.rank(begin_dist_S, colex_rank);
        end_dist_S = character_begin_dist + index.colex_gc_bwt.rank(end_dist_S, colex_rank);
#ifdef IS_LOGGED
        for (auto it {std::next(rlast)}; it != std::next(rfirst); ++it)
        {
          std::cout << *it;
        }
        std::cout
        << " -> (" << colex_rank << ": " << index.colex_lex_permutation[colex_rank] << ")"
        << " -> (" << begin_dist_S << ", " << end_dist_S << ")\n";
#endif
      }
    }
  }
  lookup_grammar_rule
  (
    std::begin(index.grammar_rules),
    index.lex_trie_root,
    std::next(rlast),
    std::next(rbegin),
    1,
    leftmost_lex_rank,
    rightmost_lex_rank
  );
  if (leftmost_lex_rank != 0)
  {
    begin_dist_L = index.lex_gc_character_bucket_end_dists[leftmost_lex_rank - 1];
    end_dist_L = index.lex_gc_character_bucket_end_dists[rightmost_lex_rank];
#ifdef IS_LOGGED
    for (auto it {std::next(rlast)}; it != std::next(rbegin); ++it)
    {
      std::cout << *it;
    }
    std::cout << " -> (" << leftmost_lex_rank << ", " << rightmost_lex_rank << ")";
    std::cout << " -> (" << begin_dist_L << ", " << end_dist_L << ")\n";
#endif
  }
  return rlast;
}

template <typename pattern_iterator_type>
uint64_t count
(
  gc_index const &index,
  pattern_iterator_type pattern_begin,
  pattern_iterator_type pattern_end
)
{
  uint64_t begin_dist_L {0};
  uint64_t end_dist_L {0};
  uint64_t begin_dist_S {0};
  uint64_t end_dist_S {0};
  auto rbegin {std::prev(pattern_end)};
  auto rend {std::prev(pattern_begin)};
  auto rfirst {rbegin};
  auto rlast {rbegin};
#ifdef IS_LOGGED
  for (auto it {pattern_begin}; it != pattern_end; ++it)
  {
    std::cout << *it;
  }
  std::cout << '\n';
#endif
  rlast = backward_search_pattern_suffix(index, begin_dist_L, end_dist_L, begin_dist_S, end_dist_S, rbegin, rend);
  if (rlast != rend)
  {
    while
    (
      (begin_dist_L != end_dist_L)
      ||
      (begin_dist_S != end_dist_S)
    )
    {
      calculate_sl_factor(rfirst, rlast, rend);
      if (rlast != rend)
      {
        backward_search_exact_sl_factor(index, begin_dist_L, end_dist_L, begin_dist_S, end_dist_S, rfirst, rlast);
      }
      else
      {
        return backward_search_pattern_prefix(index, begin_dist_L, end_dist_L, begin_dist_S, end_dist_S, rfirst, rlast);
      }
    }
  }
  return ((end_dist_L - begin_dist_L) + (end_dist_S - begin_dist_S));
}

template <typename text_type>
uint64_t calculate_max_sl_factor_size (text_type const &text)
{
  uint64_t max_size {0};
  auto rfirst {std::prev(std::end(text))};
  auto rlast {std::prev(std::end(text))};
  auto rend {std::prev(std::begin(text))};
  while (rlast != rend)
  {
    calculate_sl_factor(rfirst, rlast, rend);
    auto size {static_cast<uint64_t>(std::distance(rlast, rfirst))};
    if (max_size < size)
    {
      max_size = size;
    }
  }
  return max_size;
}

void calculate_sl_factor_size_number_pairs
(
  std::map<uint64_t, uint64_t> &size_number_pairs,
  std::string const text_path
)
{
  sdsl::int_vector<8> text;
  sdsl::load_vector_from_file(text, text_path);
  auto rfirst {std::prev(std::end(text))};
  auto rlast {std::prev(std::end(text))};
  auto rend {std::prev(std::begin(text))};
  while (rlast != rend)
  {
    calculate_sl_factor(rfirst, rlast, rend);
    auto size {static_cast<uint64_t>(std::distance(rlast, rfirst))};
    if (size_number_pairs.find(size) == std::end(size_number_pairs))
    {
      size_number_pairs[size] = 1;
    }
    else
    {
      ++size_number_pairs[size];
    }
  }
  return;
}
}

#endif
