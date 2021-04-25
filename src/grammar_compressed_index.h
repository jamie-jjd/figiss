#pragma once

#include <deque>
#include <map>
#include <memory>

#include "utility.h"

namespace project
{
constexpr uint64_t S {1};
constexpr uint64_t L {0};

template <typename SlTypesIterator>
constexpr bool IsRightmostLType (SlTypesIterator it)
{
  return
  (
    (*it == L)
    &&
    (*std::next(it) == S)
  );
}

template <typename SlTypesIterator>
constexpr bool IsLeftmostSType (SlTypesIterator it)
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
  typename Text,
  typename SlTypes
>
void CalculateSlTypes
(
  Text const &text,
  SlTypes &sl_types
)
{
  sl_types.resize(std::size(text));
  sdsl::util::set_to_value(sl_types, S);
  auto text_it {std::prev(std::end(text))};
  auto text_last_it {std::begin(text)};
  auto sl_types_it {std::prev(std::end(sl_types))};
  while (text_it != text_last_it)
  {
    if
    (
      (*std::prev(text_it) > *text_it)
      ||
      (
        (*std::prev(text_it) == *text_it)
        &&
        (*sl_types_it == L)
      )
    )
    {
      *std::prev(sl_types_it) = L;
    }
    --sl_types_it;
    --text_it;
  }
  return;
}

template
<
  typename Text,
  typename CharacterBucketOffsets
>
void CalculateCharacterBucketEndOffsets
(
  Text const &text,
  CharacterBucketOffsets &character_bucket_offsets
)
{
  auto begin {std::begin(character_bucket_offsets)};
  auto end {std::end(character_bucket_offsets)};
  sdsl::util::set_to_value(character_bucket_offsets, 0);
  for (auto const &character : text)
  {
    ++character_bucket_offsets[character];
  }
  std::partial_sum(begin, end, begin);
  return;
}

template
<
  typename Text,
  typename CharacterBucketOffsets
>
void CalculateCharacterBucketBeginOffsets
(
  Text const &text,
  CharacterBucketOffsets &character_bucket_offsets
)
{
  CalculateCharacterBucketEndOffsets(text, character_bucket_offsets);
  auto it {std::prev(std::end(character_bucket_offsets))};
  auto last_it {std::begin(character_bucket_offsets)};
  while (it != last_it)
  {
    *it = *std::prev(it);
    --it;
  }
  *last_it = 0;
  return;
}

template
<
  typename Text,
  typename SlTypes,
  typename CharacterBucketOffsets,
  typename TextOffsets
>
void BucketSortRightmostLTypeCharacters
(
  Text const &text,
  SlTypes const &sl_types,
  CharacterBucketOffsets &character_bucket_offsets,
  TextOffsets &text_offsets
)
{
  auto text_it {std::begin(text)};
  auto text_end {std::end(text)};
  while (text_it != text_end)
  {
    auto text_offset {std::distance(std::begin(text), text_it)};
    if (IsRightmostLType(std::next(std::begin(sl_types), text_offset)))
    {
      text_offsets[character_bucket_offsets[*text_it]++] = text_offset;
    }
    ++text_it;
  }
  return;
}

template
<
  typename Text,
  typename SlTypes,
  typename CharacterBucketOffsets,
  typename TextOffsets
>
void InduceSortLTypeCharacters
(
  Text const &text,
  SlTypes const &sl_types,
  CharacterBucketOffsets &character_bucket_offsets,
  TextOffsets &text_offsets
)
{
  auto text_offsets_it {std::next(std::begin(text_offsets))};
  auto text_offsets_last_it {std::end(text_offsets)};
  auto invalid_text_offset {std::size(text)};
  while (text_offsets_it != text_offsets_last_it)
  {
    auto text_offset {*text_offsets_it};
    if
    (
      (text_offset != invalid_text_offset)
      &&
      (text_offset != 0)
      &&
      (sl_types[text_offset - 1] == L)
    )
    {
      text_offsets[character_bucket_offsets[text[text_offset - 1]]++] = (text_offset - 1);
      *text_offsets_it = invalid_text_offset;
    }
    ++text_offsets_it;
  }
  return;
}

template
<
  typename Text,
  typename SlTypes,
  typename CharacterBucketOffsets,
  typename TextOffsets
>
void InduceSortSTypeCharacters
(
  Text const &text,
  SlTypes const &sl_types,
  CharacterBucketOffsets &character_bucket_offsets,
  TextOffsets &text_offsets
)
{
  auto text_offsets_it {std::prev(std::end(text_offsets))};
  auto text_offsets_last_it {std::begin(text_offsets)};
  auto invalid_text_offset {std::size(text)};
  while (text_offsets_it != text_offsets_last_it)
  {
    auto text_offset {*text_offsets_it};
    if
    (
      (text_offset != invalid_text_offset)
      &&
      (text_offset != 0)
      &&
      (sl_types[text_offset - 1] == S)
    )
    {
      text_offsets[--character_bucket_offsets[text[text_offset - 1]]] = (text_offset - 1);
      *text_offsets_it = invalid_text_offset;
    }
    --text_offsets_it;
  }
  return;
}

template <typename RandomAccessIterator>
auto MoveVaildEntriesToFront
(
  RandomAccessIterator first_it,
  RandomAccessIterator last_it,
  uint64_t const invalid_value
)
{
  auto it {first_it};
  auto new_last_it {first_it};
  while (it != last_it)
  {
    if (*it != invalid_value)
    {
      uint64_t temporary_value {*it};
      *it = *new_last_it;
      *new_last_it = temporary_value;
      ++new_last_it;
    }
    ++it;
  }
  return new_last_it;
}

template
<
  typename Text,
  typename SlTypes,
  typename GrammarRuleBeginOffsetsIterator,
  typename TemporaryLexTextIterator,
  typename GrammarRuleCounts
>
void CalculateGrammarRuleCountsBeginOffsetsAndTemporaryLexText
(
  Text const &text,
  SlTypes const &sl_types,
  GrammarRuleCounts &rule_counts,
  GrammarRuleBeginOffsetsIterator begin_offsets_begin,
  GrammarRuleBeginOffsetsIterator begin_offsets_end,
  TemporaryLexTextIterator temporary_compressed_text_begin
)
{
  std::deque<uint64_t> counts;
  auto invalid_text_offset {std::size(text)};
  uint64_t lex_rank {0};
  auto prev_rule_it {std::prev(std::end(text))};
  auto prev_sl_types_it {std::prev(std::end(sl_types))};
  auto begin_offsets_it {begin_offsets_begin};
  while (begin_offsets_it != begin_offsets_end)
  {
    uint64_t begin_offset {*begin_offsets_it};
    auto rule_it {std::next(std::begin(text), begin_offset)};
    auto sl_types_it {std::next(std::begin(sl_types), begin_offset)};
    if
    (
      (*prev_rule_it == *rule_it)
      &&
      (*prev_sl_types_it == *sl_types_it)
    )
    {
      do
      {
        ++prev_rule_it;
        ++rule_it;
        ++prev_sl_types_it;
        ++sl_types_it;
      }
      while
      (
        !IsLeftmostSType(prev_sl_types_it)
        &&
        !IsLeftmostSType(sl_types_it)
        &&
        (*prev_rule_it == *rule_it)
        &&
        (*prev_sl_types_it == *sl_types_it)
      );
      if
      (
        IsLeftmostSType(prev_sl_types_it)
        &&
        IsLeftmostSType(sl_types_it)
      )
      {
        --lex_rank;
        *begin_offsets_it = invalid_text_offset;
        ++counts.back();
      }
      else
      {
        counts.emplace_back(1);
      }
    }
    else
    {
      counts.emplace_back(1);
    }
    *std::next(temporary_compressed_text_begin, (begin_offset + 1) / 2) = ++lex_rank;
    prev_rule_it = std::next(std::begin(text), begin_offset);
    prev_sl_types_it = std::next(std::begin(sl_types), begin_offset);
    ++begin_offsets_it;
  }
  rule_counts.resize(std::size(counts));
  auto rule_counts_it {std::begin(rule_counts)};
  for (auto const &count : counts)
  {
    *rule_counts_it++ = count;
  }
  return;
}

template
<
  typename GrammarRuleSizes,
  typename SlTypes,
  typename GrammarRuleBeginOffsetsIterator
>
void CalculateGrammarRuleSizes
(
  GrammarRuleSizes &sizes,
  SlTypes const &sl_types,
  GrammarRuleBeginOffsetsIterator begin_offsets_begin,
  GrammarRuleBeginOffsetsIterator begin_offsets_end,
  uint64_t const invalid_value
)
{
  begin_offsets_end = MoveVaildEntriesToFront(begin_offsets_begin, begin_offsets_end, invalid_value);
  sizes.resize(std::distance(begin_offsets_begin, begin_offsets_end));
  auto sizes_it {std::begin(sizes)};
  auto begin_offsets_it {begin_offsets_begin};
  while (begin_offsets_it != begin_offsets_end)
  {
    auto sl_types_first_it {std::next(std::begin(sl_types), *begin_offsets_it)};
    auto sl_types_last_it {std::next(sl_types_first_it)};
    while (!IsLeftmostSType(sl_types_last_it))
    {
      ++sl_types_last_it;
    }
    *sizes_it = std::distance(sl_types_first_it, sl_types_last_it);
    ++sizes_it;
    ++begin_offsets_it;
  }
  return;
}

template
<
  typename Text,
  typename GrammarRuleSizes,
  typename GrammarRules,
  typename GrammarRuleBeginOffsetsIterator
>
void CalculateGrammarRules
(
  Text const &text,
  GrammarRuleSizes const &sizes,
  GrammarRules &rules,
  GrammarRuleBeginOffsetsIterator begin_offsets_begin
)
{
  rules.resize(std::accumulate(std::begin(sizes), std::end(sizes), 0));
  auto rules_it {std::begin(rules)};
  auto begin_offsets_it {begin_offsets_begin};
  auto sizes_it {std::begin(sizes)};
  auto sizes_end {std::end(sizes)};
  while (sizes_it != sizes_end)
  {
    auto rule_it {std::next(std::begin(text), *begin_offsets_it)};
    auto rule_end {std::next(std::begin(text), *begin_offsets_it + *sizes_it)};
    while (rule_it != rule_end)
    {
      *rules_it = *rule_it;
      ++rules_it;
      ++rule_it;
    }
    ++begin_offsets_it;
    ++sizes_it;
  }
  return;
}

template
<
  typename LexText,
  typename TemporaryLexTextIterator
>
void CalculateLexText
(
  LexText &text,
  uint64_t const width,
  TemporaryLexTextIterator begin,
  TemporaryLexTextIterator end,
  uint64_t const invalid_value
)
{
  end = MoveVaildEntriesToFront(begin, end, invalid_value);
  auto size {std::distance(begin, end) + 1};
  {
    text.width(width);
    text.resize(size);
    std::copy(begin, end, std::begin(text));
    *std::prev(std::end(text)) = 0;
  }
  return;
}

struct DynamicGrammarTrie
{
  using EdgeRange = std::pair<uint64_t, uint64_t>;
  using RankRange = std::pair<uint64_t, uint64_t>;

  struct Node
  {
    std::map<uint8_t, std::shared_ptr<Node>> branches;
    EdgeRange edge_range;
    RankRange rank_range;
    uint64_t count;

    Node () = default;

    Node
    (
      EdgeRange const edge_range_,
      RankRange const rank_range_
    )
    : edge_range {edge_range_},
      rank_range {rank_range_},
      count {}
    {
    }

    Node
    (
      EdgeRange const edge_range_,
      uint64_t const count_
    )
    : edge_range {edge_range_},
      rank_range {},
      count {count_}
    {
    }
  };

  using NodePointer = std::shared_ptr<Node>;

  NodePointer root;
  uint64_t size;
  int64_t step;

  DynamicGrammarTrie (int64_t const step_ = 1)
  : root {std::make_shared<Node>()},
    size {1},
    step {step_}
  {
  }
};

template
<
  typename File,
  typename Labels
>
void PrintDynamicGrammarTrie
(
  File &file,
  Labels const &labels,
  DynamicGrammarTrie const &trie
)
{
  file << trie.size << "\n";
  for (uint64_t i {}; i != std::size(labels); ++i)
  {
    std::cout << i << " ";
  }
  std::cout << "\n";
  Print(file, labels);
  std::deque<std::pair<DynamicGrammarTrie::NodePointer, uint64_t>> nodes;
  nodes.emplace_back(trie.root, 0);
  while (!nodes.empty())
  {
    auto node {std::get<0>(nodes.back())};
    auto depth {std::get<1>(nodes.back())};
    nodes.pop_back();
    file << depth << ":";
    file << "[" << std::get<0>(node->edge_range) << "," << std::get<1>(node->edge_range) << "]";
    file << "[" << std::get<0>(node->rank_range) << "," << std::get<1>(node->rank_range) << "]";
    file << "(" << node->count << ")";
    file << "\n";
    if (!node->branches.empty())
    {
      auto branches_it {std::rbegin(node->branches)};
      auto branches_end {std::rend(node->branches)};
      while (branches_it != branches_end)
      {
        nodes.emplace_back(std::get<1>(*branches_it), (depth + 1));
        ++branches_it;
      }
    }
  }
  return;
}

template <typename GrammarRulesIterator>
void InsertGrammarRuleSuffixAndCountIntoDynamicGrammarTrie
(
  DynamicGrammarTrie &trie,
  GrammarRulesIterator rules_begin,
  GrammarRulesIterator first,
  GrammarRulesIterator last,
  uint64_t const count
)
{
  auto node {trie.root};
  auto it {first};
  while (true)
  {
    auto character {*it};
    auto branches_it {node->branches.find(character)};
    if (branches_it == std::end(node->branches))
    {
      node->branches[character] = std::make_shared<DynamicGrammarTrie::Node>
      (
        DynamicGrammarTrie::EdgeRange
        {
          std::distance(rules_begin, it),
          (std::distance(rules_begin, last) - trie.step),
        },
        count
      );
      ++(trie.size);
      return;
    }
    else
    {
      auto child_node {std::get<1>(*branches_it)};
      auto edge_it {std::next(rules_begin, std::get<0>(child_node->edge_range))};
      auto edge_end {std::next(rules_begin, std::get<1>(child_node->edge_range) + trie.step)};
      while
      (
        (it != last)
        &&
        (edge_it != edge_end)
        &&
        (*it == *edge_it)
      )
      {
        it += trie.step;
        edge_it += trie.step;
      }
      if (edge_it == edge_end)
      {
        if (it != last)
        {
          node = child_node;
        }
        else
        {
          child_node->count += count;
          return;
        }
      }
      else
      {
        auto edge_it_offset {std::distance(rules_begin, edge_it)};
        auto internal_node
        {
          std::make_shared<DynamicGrammarTrie::Node>
          (
            DynamicGrammarTrie::EdgeRange
            {
              std::get<0>(child_node->edge_range),
              (edge_it_offset - trie.step)
            },
            0
          )
        };
        ++(trie.size);
        std::get<1>(*branches_it) = internal_node;
        internal_node->branches[*edge_it] = child_node;
        std::get<0>(child_node->edge_range) = edge_it_offset;
        if (it != last)
        {
          node = internal_node;
        }
        else
        {
          internal_node->count += count;
          return;
        }
      }
    }
  }
  return;
}

template
<
  typename GrammarRuleSizes,
  typename GrammarRules,
  typename GrammarRuleCounts
>
void InsertGrammarRuleSuffixesAndCountsIntoDynamicGrammarTrie
(
  GrammarRuleSizes const &rule_sizes,
  GrammarRules const &rules,
  GrammarRuleCounts const &rule_counts,
  DynamicGrammarTrie &lex_grammar_count_trie
)
{
  auto rule_sizes_it {std::begin(rule_sizes)};
  auto rules_it {std::begin(rules)};
  auto rule_counts_it {std::begin(rule_counts)};
  while (rules_it != std::end(rules))
  {
    auto rule_begin {rules_it};
    auto rule_end {std::next(rule_begin, *rule_sizes_it)};
    while (rule_begin != rule_end)
    {
      InsertGrammarRuleSuffixAndCountIntoDynamicGrammarTrie
      (
        lex_grammar_count_trie,
        std::begin(rules),
        rule_begin,
        rule_end,
        *rule_counts_it
      );
      ++rule_begin;
    }
    ++rule_sizes_it;
    rules_it = rule_end;
    ++rule_counts_it;
  }
  return;
}

void CalculateCumulativeGrammarCount (DynamicGrammarTrie &trie)
{
  std::deque<std::pair<DynamicGrammarTrie::NodePointer, bool>> nodes;
  nodes.emplace_back(trie.root, true);
  while (!nodes.empty())
  {
    auto node {std::get<0>(nodes.back())};
    auto &is_forward {std::get<1>(nodes.back())};
    if (is_forward)
    {
      is_forward = false;
      if (!node->branches.empty())
      {
        auto branches_it {std::begin(node->branches)};
        auto branches_end {std::end(node->branches)};
        while (branches_it != branches_end)
        {
          nodes.emplace_back(std::get<1>(*branches_it), true);
          ++branches_it;
        }
      }
    }
    else
    {
      if (!node->branches.empty())
      {
        auto branches_it {std::begin(node->branches)};
        auto branches_end {std::end(node->branches)};
        while (branches_it != branches_end)
        {
          node->count += std::get<1>(*branches_it)->count;
          ++branches_it;
        }
      }
      nodes.pop_back();
    }
  }
  return;
}

template <typename GrammarRulesIterator>
void InsertGrammarRuleAndRankIntoDynamicGrammarTrie
(
  DynamicGrammarTrie &trie,
  GrammarRulesIterator rules_begin,
  GrammarRulesIterator first,
  GrammarRulesIterator last,
  uint64_t const rank
)
{
  auto node {trie.root};
  auto it {first};
  while (true)
  {
    auto character {*it};
    auto branches_it {node->branches.find(character)};
    if (branches_it == std::end(node->branches))
    {
      node->branches[character] = std::make_shared<DynamicGrammarTrie::Node>
      (
        DynamicGrammarTrie::EdgeRange
        {
          std::distance(rules_begin, it),
          (std::distance(rules_begin, last) - trie.step)
        },
        DynamicGrammarTrie::RankRange{rank, rank}
      );
      ++(trie.size);
      return;
    }
    else
    {
      auto child_node {std::get<1>(*branches_it)};
      auto edge_it {std::next(rules_begin, std::get<0>(child_node->edge_range))};
      auto edge_end {std::next(rules_begin, std::get<1>(child_node->edge_range) + trie.step)};
      while
      (
        (it != last)
        &&
        (edge_it != edge_end)
        &&
        (*it == *edge_it)
      )
      {
        it += trie.step;
        edge_it += trie.step;
      }
      if (edge_it == edge_end)
      {
        if (it != last)
        {
          node = child_node;
        }
        else
        {
          child_node->rank_range = {rank, rank};
          return;
        }
      }
      else
      {
        auto edge_it_offset {std::distance(rules_begin, edge_it)};
        auto internal_node
        {
          std::make_shared<DynamicGrammarTrie::Node>
          (
            DynamicGrammarTrie::EdgeRange
            {
              std::get<0>(child_node->edge_range),
              (edge_it_offset - trie.step)
            },
            DynamicGrammarTrie::RankRange{}
          )
        };
        ++(trie.size);
        std::get<1>(*branches_it) = internal_node;
        internal_node->branches[*edge_it] = child_node;
        std::get<0>(child_node->edge_range) = edge_it_offset;
        if (it != last)
        {
          node = internal_node;
        }
        else
        {
          internal_node->rank_range = {rank, rank};
          return;
        }
      }
    }
  }
  return;
}

template
<
  typename GrammarRuleSizes,
  typename GrammarRules
>
void InsertGrammarRulesIntoDynamicGrammarTries
(
  GrammarRuleSizes const &rule_sizes,
  GrammarRules const &rules,
  DynamicGrammarTrie &lex_grammar_rank_trie,
  DynamicGrammarTrie &colex_grammar_rank_trie
)
{
  uint64_t lex_rank {1};
  auto rule_sizes_it {std::begin(rule_sizes)};
  auto rules_it {std::begin(rules)};
  while (rules_it != std::end(rules))
  {
    auto rule_begin {rules_it};
    auto rule_end {std::next(rule_begin, *rule_sizes_it)};
    InsertGrammarRuleAndRankIntoDynamicGrammarTrie
    (
      lex_grammar_rank_trie,
      std::begin(rules),
      rule_begin,
      rule_end,
      lex_rank
    );
    InsertGrammarRuleAndRankIntoDynamicGrammarTrie
    (
      colex_grammar_rank_trie,
      std::begin(rules),
      std::prev(rule_end),
      std::prev(rule_begin),
      lex_rank
    );
    ++lex_rank;
    ++rule_sizes_it;
    rules_it = rule_end;
  }
  return;
}

void CalculateCumulativeLexRankRanges (DynamicGrammarTrie &lex_grammar_rank_trie)
{
  std::deque<std::pair<DynamicGrammarTrie::NodePointer, bool>> nodes;
  nodes.emplace_back(lex_grammar_rank_trie.root, true);
  while (!nodes.empty())
  {
    auto node {std::get<0>(nodes.back())};
    auto &is_forward {std::get<1>(nodes.back())};
    if (is_forward)
    {
      is_forward = false;
      if (!node->branches.empty())
      {
        auto branches_it {std::rbegin(node->branches)};
        auto branches_end {std::rend(node->branches)};
        while (branches_it != branches_end)
        {
          nodes.emplace_back(std::get<1>(*branches_it), true);
          ++branches_it;
        }
      }
    }
    else
    {
      if (!node->branches.empty())
      {
        auto first_child_node {std::get<1>(*std::begin(node->branches))};
        auto last_child_node {std::get<1>(*std::rbegin(node->branches))};
        if (std::get<0>(node->rank_range) == 0)
        {
          std::get<0>(node->rank_range) = std::get<0>(first_child_node->rank_range);
        }
        std::get<1>(node->rank_range) = std::get<1>(last_child_node->rank_range);
      }
      nodes.pop_back();
    }
  }
  return;
}

template <typename LexToColex>
void CalculateCumulativeColexRankRangesAndLexToColex
(
  DynamicGrammarTrie colex_grammar_rank_trie,
  LexToColex &lex_to_colex
)
{
  uint64_t colex_rank {1};
  std::deque<std::pair<DynamicGrammarTrie::NodePointer, bool>> nodes;
  nodes.emplace_back(colex_grammar_rank_trie.root, true);
  while (!nodes.empty())
  {
    auto node {std::get<0>(nodes.back())};
    auto &is_forward {std::get<1>(nodes.back())};
    if (is_forward)
    {
      is_forward = false;
      auto &leftmost_rank {std::get<0>(node->rank_range)};
      auto &rightmost_rank {std::get<1>(node->rank_range)};
      if (leftmost_rank != 0)
      {
        lex_to_colex[leftmost_rank] = colex_rank;
        leftmost_rank = rightmost_rank = colex_rank++;
      }
      if (!node->branches.empty())
      {
        auto branches_it {std::rbegin(node->branches)};
        auto branches_end {std::rend(node->branches)};
        while (branches_it != branches_end)
        {
          nodes.emplace_back(std::get<1>(*branches_it), true);
          ++branches_it;
        }
      }
    }
    else
    {
      if (!node->branches.empty())
      {
        auto first_child_node {std::get<1>(*std::begin(node->branches))};
        auto last_child_node {std::get<1>(*std::rbegin(node->branches))};
        if (std::get<0>(node->rank_range) == 0)
        {
          std::get<0>(node->rank_range) = std::get<0>(first_child_node->rank_range);
        }
        std::get<1>(node->rank_range) = std::get<1>(last_child_node->rank_range);
      }
      nodes.pop_back();
    }
  }
  return;
}

template
<
  typename LexText,
  typename LexToColex,
  typename ColexBwt
>
void CalculateColexBwt
(
  LexText const &lex_text,
  LexToColex const &lex_to_colex,
  ColexBwt &colex_bwt
)
{
  sdsl::int_vector<> buffer;
  sdsl::qsufsort::construct_sa(buffer, lex_text);
  auto buffer_it {std::begin(buffer)};
  while (buffer_it != std::end(buffer))
  {
    if (*buffer_it != 0)
    {
      *buffer_it = lex_to_colex[lex_text[(*buffer_it - 1)]];
    }
    ++buffer_it;
  }
  sdsl::construct_im(colex_bwt, buffer);
  return;
}

struct StaticGrammarTrie
{
  sdsl::int_vector<> branch_end_offsets;
  sdsl::int_vector<> branch_characters;
  sdsl::int_vector<> edge_ranges;
  sdsl::int_vector<> rank_ranges;
  sdsl::int_vector<> counts;
};

template
<
  typename File,
  typename Labels
>
void PrintStaticGrammarTrie
(
  File &file,
  Labels const &labels,
  StaticGrammarTrie const &trie
)
{
  for (uint64_t i {}; i != std::size(labels); ++i)
  {
    std::cout << i << " ";
  }
  std::cout << "\n";
  Print(file, labels);
  Print(file, trie.branch_end_offsets);
  Print(file, trie.branch_characters);
  Print(file, trie.edge_ranges);
  Print(file, trie.rank_ranges);
  Print(file, trie.counts);
  return;
}

template
<
  typename FromVector,
  typename ToVector
>
void InitializeAndCopy
(
  FromVector const &from,
  ToVector &to
)
{
  to.width(sdsl::bits::hi(*std::max_element(std::begin(from), std::end(from))) + 1);
  to.resize(std::size(from));
  std::copy(std::begin(from), std::end(from), std::begin(to));
  return;
}

void ConvertDynamicToStaticGrammarTrie
(
  DynamicGrammarTrie const &dynamic_trie,
  StaticGrammarTrie &static_trie,
  bool const is_rank = true
)
{
  std::deque<uint64_t> branch_end_offsets;
  std::deque<uint8_t> branch_characters;
  std::deque<uint64_t> edge_ranges;
  std::deque<uint64_t> rank_ranges;
  std::deque<uint64_t> counts;
  std::deque<DynamicGrammarTrie::NodePointer> nodes;
  nodes.emplace_back(dynamic_trie.root);
  while (!nodes.empty())
  {
    auto node {nodes.front()}; nodes.pop_front();
    if (!node->branches.empty())
    {
      branch_end_offsets.emplace_back(std::size(node->branches));
      auto branches_it {std::begin(node->branches)};
      auto branches_end {std::end(node->branches)};
      while (branches_it != branches_end)
      {
        auto character {std::get<0>(*branches_it)};
        auto child_node {std::get<1>(*branches_it)};
        branch_characters.emplace_back(character);
        edge_ranges.emplace_back(std::get<0>(child_node->edge_range));
        edge_ranges.emplace_back(std::get<1>(child_node->edge_range));
        if (is_rank)
        {
          rank_ranges.emplace_back(std::get<0>(child_node->rank_range));
          rank_ranges.emplace_back(std::get<1>(child_node->rank_range));
        }
        else
        {
          counts.emplace_back(child_node->count);
        }
        nodes.emplace_back(child_node);
        ++branches_it;
      }
    }
  }
  std::partial_sum
  (
    std::begin(branch_end_offsets),
    std::end(branch_end_offsets),
    std::begin(branch_end_offsets)
  );
  InitializeAndCopy(branch_end_offsets, static_trie.branch_end_offsets);
  InitializeAndCopy(branch_characters, static_trie.branch_characters);
  InitializeAndCopy(edge_ranges, static_trie.edge_ranges);
  if (is_rank)
  {
    InitializeAndCopy(rank_ranges, static_trie.rank_ranges);
  }
  else
  {
    InitializeAndCopy(counts, static_trie.counts);
  }
  return;
}

struct Index
{
  sdsl::int_vector<8> grammar_rules;
  StaticGrammarTrie lex_grammar_count_trie;
  StaticGrammarTrie lex_grammar_rank_trie;
  StaticGrammarTrie colex_grammar_rank_trie;
  sdsl::int_vector<> colex_to_lex;
  sdsl::int_vector<> lex_rank_bucket_end_offsets;
  // RunLengthWaveletTree colex_bwt;
};

template <typename File>
void PrintIndex
(
  File &file,
  Index &index
)
{
  file << index.grammar_rules << "\n";
  PrintStaticGrammarTrie(file, index.grammar_rules, index.lex_grammar_count_trie);
  PrintStaticGrammarTrie(file, index.grammar_rules, index.lex_grammar_rank_trie);
  PrintStaticGrammarTrie(file, index.grammar_rules, index.colex_grammar_rank_trie);
  file << index.colex_to_lex << "\n";
  file << index.lex_rank_bucket_end_offsets << "\n";
  // file << index.colex_bwt << "\n";
  return;
}

void Construct
(
  Index &index,
  std::filesystem::path const &text_path
)
{
  sdsl::int_vector<8> text;
  {
    std::ifstream text_file {text_path};
    sdsl::load_vector_from_file(text, text_path);
    if (*std::prev(std::end(text)) != 0)
    {
      sdsl::append_zero_symbol(text);
      // std::cout << "zero symbol appended\n";
    }
    Print(std::cout, text);
  }
  sdsl::bit_vector sl_types;
  {
    CalculateSlTypes(text, sl_types);
    Print(std::cout, sl_types);
  }
  auto invalid_text_offset {std::size(text)};
  auto text_size_width {sdsl::bits::hi(std::size(text)) + 1};
  sdsl::int_vector<> text_offsets;
  {
    text_offsets.width(text_size_width);
    text_offsets.resize(std::size(text));
    sdsl::util::set_to_value(text_offsets, invalid_text_offset);
    // Print(std::cout, text_offsets);
  }
  sdsl::int_vector<> character_bucket_offsets;
  {
    character_bucket_offsets.width(text_size_width);
    character_bucket_offsets.resize(256);
  }
  {
    CalculateCharacterBucketBeginOffsets(text, character_bucket_offsets);
    // Print(std::cout, character_bucket_offsets);
    BucketSortRightmostLTypeCharacters(text, sl_types, character_bucket_offsets, text_offsets);
    // Print(std::cout, text_offsets);
    InduceSortLTypeCharacters(text, sl_types, character_bucket_offsets, text_offsets);
    // Print(std::cout, text_offsets);
    CalculateCharacterBucketEndOffsets(text, character_bucket_offsets);
    // Print(std::cout, character_bucket_offsets);
    InduceSortSTypeCharacters(text, sl_types, character_bucket_offsets, text_offsets);
    // Print(std::cout, text_offsets);
  }
  auto text_offsets_boundary {std::begin(text_offsets)};
  text_offsets_boundary = MoveVaildEntriesToFront
  (
    std::begin(text_offsets),
    std::end(text_offsets),
    invalid_text_offset
  );
  // Print(std::cout, text_offsets);
  auto grammar_rule_begin_offsets_begin {std::begin(text_offsets)};
  auto grammar_rule_begin_offsets_end {text_offsets_boundary};
  auto temporary_lex_text_begin {text_offsets_boundary};
  auto temporary_lex_text_end {std::end(text_offsets)};
  sdsl::int_vector<> grammar_rule_counts;
  CalculateGrammarRuleCountsBeginOffsetsAndTemporaryLexText
  (
    text,
    sl_types,
    grammar_rule_counts,
    grammar_rule_begin_offsets_begin,
    grammar_rule_begin_offsets_end,
    temporary_lex_text_begin
  );
  Print(std::cout, grammar_rule_counts);
  // Print(std::cout, grammar_rule_begin_offsets_begin, grammar_rule_begin_offsets_end);
  // Print(std::cout, temporary_lex_text_begin, temporary_lex_text_end);
  sdsl::int_vector<> grammar_rule_sizes;
  CalculateGrammarRuleSizes
  (
    grammar_rule_sizes,
    sl_types,
    grammar_rule_begin_offsets_begin,
    grammar_rule_begin_offsets_end,
    invalid_text_offset
  );
  Print(std::cout, grammar_rule_sizes);
  CalculateGrammarRules
  (
    text,
    grammar_rule_sizes,
    index.grammar_rules,
    grammar_rule_begin_offsets_begin
  );
  Print(std::cout, index.grammar_rules);
  auto grammar_ranks_size {std::size(grammar_rule_sizes) + 1};
  auto lex_text_width {sdsl::bits::hi(std::size(grammar_rule_sizes)) + 1};
  sdsl::int_vector<> lex_text;
  CalculateLexText
  (
    lex_text,
    lex_text_width,
    temporary_lex_text_begin,
    temporary_lex_text_end,
    invalid_text_offset
  );
  Print(std::cout, lex_text);
  {
    sdsl::util::clear(sl_types);
    sdsl::util::clear(text);
    sdsl::util::clear(text_offsets);
  }
  DynamicGrammarTrie lex_grammar_count_trie;
  DynamicGrammarTrie lex_grammar_rank_trie;
  DynamicGrammarTrie colex_grammar_rank_trie {-1};
  InsertGrammarRuleSuffixesAndCountsIntoDynamicGrammarTrie
  (
    grammar_rule_sizes,
    index.grammar_rules,
    grammar_rule_counts,
    lex_grammar_count_trie
  );
  // PrintDynamicGrammarTrie(std::cout, grammar_rules, lex_grammar_count_trie);
  CalculateCumulativeGrammarCount(lex_grammar_count_trie);
  PrintDynamicGrammarTrie(std::cout, index.grammar_rules, lex_grammar_count_trie);
  InsertGrammarRulesIntoDynamicGrammarTries
  (
    grammar_rule_sizes,
    index.grammar_rules,
    lex_grammar_rank_trie,
    colex_grammar_rank_trie
  );
  // PrintDynamicGrammarTrie(std::cout, index.grammar_rules, lex_grammar_rank_trie);
  // PrintDynamicGrammarTrie(std::cout, index.grammar_rules, colex_grammar_rank_trie);
  CalculateCumulativeLexRankRanges(lex_grammar_rank_trie);
  PrintDynamicGrammarTrie(std::cout, index.grammar_rules, lex_grammar_rank_trie);
  sdsl::int_vector<> lex_to_colex;
  lex_to_colex.width(lex_text_width);
  lex_to_colex.resize(grammar_ranks_size);
  lex_to_colex[0] = 0;
  CalculateCumulativeColexRankRangesAndLexToColex(colex_grammar_rank_trie, lex_to_colex);
  // Print(std::cout, lex_to_colex);
  PrintDynamicGrammarTrie(std::cout, index.grammar_rules, colex_grammar_rank_trie);
  ConvertDynamicToStaticGrammarTrie
  (
    lex_grammar_count_trie,
    index.lex_grammar_count_trie,
    false
  );
  PrintStaticGrammarTrie(std::cout, index.grammar_rules, index.lex_grammar_count_trie);
  ConvertDynamicToStaticGrammarTrie
  (
    lex_grammar_rank_trie,
    index.lex_grammar_rank_trie
  );
  PrintStaticGrammarTrie(std::cout, index.grammar_rules, index.lex_grammar_rank_trie);
  ConvertDynamicToStaticGrammarTrie
  (
    colex_grammar_rank_trie,
    index.colex_grammar_rank_trie
  );
  PrintStaticGrammarTrie(std::cout, index.grammar_rules, index.colex_grammar_rank_trie);
  index.colex_to_lex.width(lex_text_width);
  index.colex_to_lex.resize(grammar_ranks_size);
  for (uint64_t rank {}; rank != std::size(lex_to_colex); ++rank)
  {
    index.colex_to_lex[lex_to_colex[rank]] = rank;
  }
  // Print(std::cout, index.colex_to_lex);
  index.lex_rank_bucket_end_offsets.width(sdsl::bits::hi(std::size(lex_text)) + 1);
  index.lex_rank_bucket_end_offsets.resize(grammar_ranks_size);
  CalculateCharacterBucketEndOffsets
  (
    lex_text,
    index.lex_rank_bucket_end_offsets
  );
  Print(std::cout, index.lex_rank_bucket_end_offsets);
  // CalculateColexBwt(lex_text, lex_to_colex, index.colex_bwt);
  return;
}

// void SerializeIndex
// (
//   Index &index,
//   std::filesystem::path index_path
// )
// {
//   std::ofstream index_file {index_path};
//   sdsl::serialize(index.grammar_rules, index_file);
//   project::SerializeStaticGrammarTrie(index.lex_grammar_count_trie, index_file);
//   project::SerializeStaticGrammarTrie(index.lex_grammar_rank_trie, index_file);
//   project::SerializeStaticGrammarTrie(index.colex_grammar_rank_trie, index_file);
//   sdsl::serialize(index.colex_to_lex, index_file);
//   sdsl::serialize(index.lex_rank_bucket_end_offsets, index_file);
//   // sdsl::serialize(index.colex_bwt, index_file);
//   return;
// }
//
// void LoadIndex
// (
//   Index &index,
//   std::filesystem::path const &index_path
// )
// {
//   std::ifstream index_file {index_path};
//   LoadStaticGrammarTrie(index.lex_grammar_count_trie);
//   LoadStaticGrammarTrie(index.lex_grammar_rank_trie);
//   LoadStaticGrammarTrie(index.colex_grammar_rank_trie);
//   index.colex_to_lex.load(index_file);
//   index.lex_rank_bucket_end_offsets.load(index_file);
//   index.colex_bwt.load(index_file);
//   return;
// }
//
// template <typename TextIterator>
// void CalculateSlFactor
// (
//   TextIterator &reverse_first_it,
//   TextIterator &reverse_last_it,
//   TextIterator reverse_end
// )
// {
//   uint64_t prev_sl_type {L};
//   reverse_first_it = reverse_last_it--;
//   while
//   (
//     (reverse_last_it != reverse_end)
//     &&
//     ! (
//         (prev_sl_type == S)
//         &&
//         (*reverse_last_it > *std::next(reverse_last_it))
//       )
//   )
//   {
//     if
//     (
//       (prev_sl_type == L)
//       &&
//       (*reverse_last_it < *std::next(reverse_last_it))
//     )
//     {
//       prev_sl_type = S;
//     }
//     --reverse_last_it;
//   }
//   return;
// }
//
// template
// <
//   typename GrammarRulesIterator,
//   typename SlFactorIterator
// >
// void LookupGrammarRule
// (
//   GrammarRulesIterator grammar_rules_begin,
//   TrieNode *node,
//   SlFactorIterator sl_factor_it,
//   SlFactorIterator sl_factor_end,
//   int64_t const offset,
//   uint64_t &leftmost_rank,
//   uint64_t &rightmost_rank
// )
// {
//   leftmost_rank = rightmost_rank = 0;
//   auto branches_it {std::end(node->branches)};
//   while
//   (
//     (branches_it = node->branches.find(*sl_factor_it))
//     !=
//     std::end(node->branches)
//   )
//   {
//     node = std::get<1>(*branches_it);
//     auto edge_it {std::next(grammar_rules_begin, node->edge_begin_offset)};
//     auto edge_end {std::next(grammar_rules_begin, node->edge_end_offset)};
//     while
//     (
//       (sl_factor_it != sl_factor_end)
//       &&
//       (edge_it != edge_end)
//       &&
//       (*sl_factor_it == *edge_it)
//     )
//     {
//       sl_factor_it += offset;
//       edge_it += offset;
//     }
//     if (sl_factor_it != sl_factor_end)
//     {
//       if (edge_it != edge_end)
//       {
//         return;
//       }
//     }
//     else
//     {
//       rightmost_rank = node->rightmost_rank;
//       leftmost_rank = node->leftmost_rank;
//       return;
//     }
//   }
//   return;
// }
//
// template <typename PatternIterator>
// void BackwardSearchPatternPrefix
// (
//   Index const &index,
//   uint64_t &pattern_range_begin_offset_L,
//   uint64_t &pattern_range_end_offset_L,
//   uint64_t &pattern_range_begin_offset_S,
//   uint64_t &pattern_range_end_offset_S,
//   PatternIterator reverse_pattern_first_it,
//   PatternIterator reverse_pattern_last_it
// )
// {
//   uint64_t leftmost_colex_rank {0};
//   uint64_t rightmost_colex_rank {0};
//   LookupGrammarRule
//   (
//     std::begin(index.grammar_rules),
//     index.colex_grammar_rank_trie_root,
//     reverse_pattern_first_it,
//     reverse_pattern_last_it,
//     -1,
//     leftmost_colex_rank,
//     rightmost_colex_rank
//   );
//   if (leftmost_colex_rank != 0)
//   {
//     if (pattern_range_begin_offset_L != pattern_range_end_offset_L)
//     {
//       pattern_range_end_offset_L = std::get<0>
//       (
//         index.colex_bwt.range_search_2d
//         (
//           pattern_range_begin_offset_L,
//           (pattern_range_end_offset_L - 1),
//           leftmost_colex_rank,
//           rightmost_colex_rank,
//           false
//         )
//       );
//       pattern_range_begin_offset_L = 0;
//     }
//     if (pattern_range_begin_offset_S != pattern_range_end_offset_S)
//     {
//       pattern_range_end_offset_S = std::get<0>
//       (
//         index.colex_bwt.range_search_2d
//         (
//           pattern_range_begin_offset_S,
//           (pattern_range_end_offset_S - 1),
//           leftmost_colex_rank,
//           rightmost_colex_rank,
//           false
//         )
//       );
//       pattern_range_begin_offset_S = 0;
//     }
//   }
//   // {
//   //   auto it {std::next(reverse_pattern_last_it)};
//   //   auto end {std::next(reverse_pattern_first_it)};
//   //   while ( it != end)
//   //   {
//   //     std::cout << *it;
//   //     ++it;
//   //   }
//   //   std::cout
//   //   << "->(" << leftmost_colex_rank << "," << rightmost_colex_rank << ")"
//   //   << "->L:(" << pattern_range_begin_offset_L << "," << pattern_range_end_offset_L << ")"
//   //   << "->S:(" << pattern_range_begin_offset_S << "," << pattern_range_end_offset_S << ")\n";
//   // }
//   return;
// }
//
// template <typename PatternIterator>
// void BackwardSearchExactSlFactor
// (
//   Index const &index,
//   uint64_t &pattern_range_begin_offset_L,
//   uint64_t &pattern_range_end_offset_L,
//   uint64_t &pattern_range_begin_offset_S,
//   uint64_t &pattern_range_end_offset_S,
//   PatternIterator reverse_pattern_first_it,
//   PatternIterator reverse_pattern_last_it
// )
// {
//   uint64_t colex_rank {0};
//   LookupGrammarRule
//   (
//     std::begin(index.grammar_rules),
//     index.colex_grammar_rank_trie_root,
//     reverse_pattern_first_it,
//     reverse_pattern_last_it,
//     -1,
//     colex_rank,
//     colex_rank
//   );
//   if (colex_rank != 0)
//   {
//     auto character_begin_offset
//     {
//       index.lex_rank_bucket_end_offsets
//       [
//         index.colex_to_lex[colex_rank] - 1
//       ]
//     };
//     if (pattern_range_begin_offset_L != pattern_range_end_offset_L)
//     {
//       pattern_range_begin_offset_L =
//       (
//         character_begin_offset
//         + index.colex_bwt.rank
//         (
//           pattern_range_begin_offset_L,
//           colex_rank
//         )
//       );
//       pattern_range_end_offset_L =
//       (
//         character_begin_offset
//         + index.colex_bwt.rank
//         (
//           pattern_range_end_offset_L,
//           colex_rank
//         )
//       );
//     }
//     if (pattern_range_begin_offset_S != pattern_range_end_offset_S)
//     {
//       pattern_range_begin_offset_S =
//       (
//         character_begin_offset
//         + index.colex_bwt.rank
//         (
//           pattern_range_begin_offset_S,
//           colex_rank
//         )
//       );
//       pattern_range_end_offset_S =
//       (
//         character_begin_offset
//         + index.colex_bwt.rank
//         (
//           pattern_range_end_offset_S,
//           colex_rank
//         )
//       );
//     }
//   }
//   else
//   {
//     pattern_range_begin_offset_L = pattern_range_end_offset_L;
//     pattern_range_begin_offset_S = pattern_range_end_offset_S;
//   }
//   // {
//   //   auto it {std::next(reverse_pattern_last_it)};
//   //   auto end {std::next(reverse_pattern_first_it)};
//   //   while (it != end)
//   //   {
//   //     std::cout << *it;
//   //     ++it;
//   //   }
//   //   std::cout
//   //   << "->(" << colex_rank << ":" << index.colex_to_lex[colex_rank] << ")"
//   //   << "->L:(" << pattern_range_begin_offset_L << "," << pattern_range_end_offset_L << ")"
//   //   << "->S:(" << pattern_range_begin_offset_S << "," << pattern_range_end_offset_S << ")\n";
//   // }
//   return;
// }
//
// template <typename PatternIterator>
// auto CalculateReversePatternLastIteratorOfPatternSuffixS
// (
//   PatternIterator reverse_begin,
//   PatternIterator reverse_end
// )
// {
//   auto it {reverse_begin};
//   while
//   (
//     (std::prev(it) != reverse_end)
//     &&
//     (*std::prev(it) == *it)
//   )
//   {
//     --it;
//   }
//   if
//   (
//     (std::prev(it) != reverse_end)
//     &&
//     (*std::prev(it) < *it)
//   )
//   {
//     return reverse_begin;
//   }
//   return std::prev(it);
// }
//
// template <typename PatternIterator>
// auto BackwardSearchPatternSuffix
// (
//   Index const &index,
//   uint64_t &pattern_range_begin_offset_L,
//   uint64_t &pattern_range_end_offset_L,
//   uint64_t &pattern_range_begin_offset_S,
//   uint64_t &pattern_range_end_offset_S,
//   PatternIterator reverse_pattern_begin,
//   PatternIterator reverse_pattern_end
// )
// {
//   auto reverse_pattern_first_it {reverse_pattern_begin};
//   auto reverse_pattern_last_it
//   {
//     CalculateReversePatternLastIteratorOfPatternSuffixS
//     (
//       reverse_pattern_begin,
//       reverse_pattern_end
//     )
//   };
//   uint64_t leftmost_lex_rank {0};
//   uint64_t rightmost_lex_rank {0};
//   if
//   (
//     (reverse_pattern_last_it != reverse_pattern_begin)
//     &&
//     (reverse_pattern_last_it != reverse_pattern_end)
//   )
//   {
//     LookupGrammarRule
//     (
//       std::begin(index.grammar_rules),
//       index.lex_grammar_rank_trie_root,
//       std::next(reverse_pattern_last_it),
//       std::next(reverse_pattern_begin),
//       1,
//       leftmost_lex_rank,
//       rightmost_lex_rank
//     );
//     if (leftmost_lex_rank != 0)
//     {
//       pattern_range_begin_offset_S = index.lex_rank_bucket_end_offsets[leftmost_lex_rank - 1];
//       pattern_range_end_offset_S = index.lex_rank_bucket_end_offsets[rightmost_lex_rank];
//     }
//     // {
//     //   auto it {std::next(reverse_pattern_last_it)};
//     //   auto end {std::next(reverse_pattern_begin)};
//     //   while (it != end)
//     //   {
//     //     std::cout << *it;
//     //     ++it;
//     //   }
//     //   std::cout << "->(" << leftmost_lex_rank << "," << rightmost_lex_rank << ")";
//     //   std::cout << "->S:(" << pattern_range_begin_offset_S << "," << pattern_range_end_offset_S << ")\n";
//     // }
//     CalculateSlFactor(reverse_pattern_first_it, reverse_pattern_last_it, reverse_pattern_end);
//     if (pattern_range_begin_offset_S != pattern_range_end_offset_S)
//     {
//       uint64_t leftmost_colex_rank {0};
//       uint64_t rightmost_colex_rank {0};
//       LookupGrammarRule
//       (
//         std::begin(index.grammar_rules),
//         index.colex_grammar_rank_trie_root,
//         reverse_pattern_first_it,
//         reverse_pattern_last_it,
//         -1,
//         leftmost_colex_rank,
//         rightmost_colex_rank
//       );
//       if (leftmost_colex_rank != 0)
//       {
//         if (reverse_pattern_last_it != reverse_pattern_end)
//         {
//           auto character_begin_offset
//           {
//             index.lex_rank_bucket_end_offsets
//             [
//               index.colex_to_lex[leftmost_colex_rank] - 1
//             ]
//           };
//           pattern_range_begin_offset_S =
//           (
//             character_begin_offset
//             + index.colex_bwt.rank
//             (
//               pattern_range_begin_offset_S,
//               leftmost_colex_rank
//             )
//           );
//           pattern_range_end_offset_S =
//           (
//             character_begin_offset
//             + index.colex_bwt.rank
//             (
//               pattern_range_end_offset_S,
//               leftmost_colex_rank
//             )
//           );
//         }
//         else
//         {
//           pattern_range_end_offset_S = std::get<0>
//           (
//             index.colex_bwt.range_search_2d
//             (
//               pattern_range_begin_offset_S,
//               (pattern_range_end_offset_S - 1),
//               leftmost_colex_rank,
//               rightmost_colex_rank,
//               false
//             )
//           );
//           pattern_range_begin_offset_S = 0;
//         }
//       }
//       else
//       {
//         pattern_range_begin_offset_S = pattern_range_end_offset_S;
//       }
//       // {
//       //   auto it {std::next(reverse_pattern_last_it)};
//       //   auto end {std::next(reverse_pattern_first_it)};
//       //   while (it != end)
//       //   {
//       //     std::cout << *it;
//       //     ++it;
//       //   }
//       //   std::cout
//       //   << "->(" << leftmost_colex_rank << ":" << index.colex_to_lex[leftmost_colex_rank]
//       //   << "," << rightmost_colex_rank << ")"
//       //   << "->S:(" << pattern_range_begin_offset_S << "," << pattern_range_end_offset_S << ")\n";
//       // }
//     }
//   }
//   else if (reverse_pattern_last_it == reverse_pattern_begin)
//   {
//     CalculateSlFactor
//     (
//       reverse_pattern_first_it,
//       reverse_pattern_last_it,
//       reverse_pattern_end
//     );
//   }
//   LookupGrammarRule
//   (
//     std::begin(index.grammar_rules),
//     index.lex_grammar_rank_trie_root,
//     std::next(reverse_pattern_last_it),
//     std::next(reverse_pattern_begin),
//     1,
//     leftmost_lex_rank,
//     rightmost_lex_rank
//   );
//   if (leftmost_lex_rank != 0)
//   {
//     pattern_range_begin_offset_L = index.lex_rank_bucket_end_offsets[leftmost_lex_rank - 1];
//     pattern_range_end_offset_L = index.lex_rank_bucket_end_offsets[rightmost_lex_rank];
//   }
//   // {
//   //   auto it {std::next(reverse_pattern_last_it)};
//   //   auto end {std::next(reverse_pattern_begin)};
//   //   while (it != end)
//   //   {
//   //     std::cout << *it;
//   //     ++it;
//   //   }
//   //   std::cout << "->(" << leftmost_lex_rank << "," << rightmost_lex_rank << ")";
//   //   std::cout << "->L:(" << pattern_range_begin_offset_L << "," << pattern_range_end_offset_L << ")\n";
//   // }
//   return reverse_pattern_last_it;
// }
//
// template <typename PatternIterator>
// uint64_t Count
// (
//   Index const &index,
//   PatternIterator pattern_begin,
//   PatternIterator pattern_end
// )
// {
//   uint64_t pattern_range_begin_offset_L {0};
//   uint64_t pattern_range_end_offset_L {0};
//   uint64_t pattern_range_begin_offset_S {0};
//   uint64_t pattern_range_end_offset_S {0};
//   auto reverse_pattern_begin {std::prev(pattern_end)};
//   auto reverse_pattern_end {std::prev(pattern_begin)};
//   auto reverse_pattern_first_it {reverse_pattern_begin};
//   auto reverse_pattern_last_it {reverse_pattern_begin};
//   // {
//   //   for (auto it {pattern_begin}; it != pattern_end; ++it)
//   //   {
//   //     std::cout << *it;
//   //   }
//   //   std::cout << "\n";
//   // }
//   reverse_pattern_last_it = BackwardSearchPatternSuffix
//   (
//     index,
//     pattern_range_begin_offset_L,
//     pattern_range_end_offset_L,
//     pattern_range_begin_offset_S,
//     pattern_range_end_offset_S,
//     reverse_pattern_begin,
//     reverse_pattern_end
//   );
//   if (reverse_pattern_last_it != reverse_pattern_end)
//   {
//     while
//     (
//       (pattern_range_begin_offset_L != pattern_range_end_offset_L)
//       ||
//       (pattern_range_begin_offset_S != pattern_range_end_offset_S)
//     )
//     {
//       CalculateSlFactor(reverse_pattern_first_it, reverse_pattern_last_it, reverse_pattern_end);
//       if (reverse_pattern_last_it!= reverse_pattern_end)
//       {
//         BackwardSearchExactSlFactor
//         (
//           index,
//           pattern_range_begin_offset_L,
//           pattern_range_end_offset_L,
//           pattern_range_begin_offset_S,
//           pattern_range_end_offset_S,
//           reverse_pattern_first_it,
//           reverse_pattern_last_it
//         );
//       }
//       else
//       {
//         BackwardSearchPatternPrefix
//         (
//           index,
//           pattern_range_begin_offset_L,
//           pattern_range_end_offset_L,
//           pattern_range_begin_offset_S,
//           pattern_range_end_offset_S,
//           reverse_pattern_first_it,
//           reverse_pattern_last_it
//         );
//         break;
//       }
//     }
//   }
//   return
//   (
//     (pattern_range_end_offset_L - pattern_range_begin_offset_L)
//     +
//     (pattern_range_end_offset_S - pattern_range_begin_offset_S)
//   );
// }
//
// template <typename Text>
// uint64_t CalculateMaxSlFactorSize (Text const &text)
// {
//   uint64_t max_size {0};
//   auto reverse_first_it {std::prev(std::end(text))};
//   auto reverse_last_it {std::prev(std::end(text))};
//   auto reverse_end {std::prev(std::begin(text))};
//   while (reverse_last_it != reverse_end)
//   {
//     CalculateSlFactor(reverse_first_it, reverse_last_it, reverse_end);
//     auto size {static_cast<uint64_t>(std::distance(reverse_last_it, reverse_first_it))};
//     if (max_size < size)
//     {
//       max_size = size;
//     }
//   }
//   return max_size;
// }
}
