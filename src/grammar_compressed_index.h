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
  return ((*it == L) && (*std::next(it) == S));
}

template <typename SlTypesIterator>
constexpr bool IsLeftmostSType (SlTypesIterator it)
{
  return ((*it == S) && (*std::prev(it) == L));
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
    if ((*std::prev(text_it) > *text_it) || ((*std::prev(text_it) == *text_it) && (*sl_types_it == L)))
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
    *it-- = *std::prev(it);
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
    if((text_offset != invalid_text_offset) && (text_offset != 0) && (sl_types[text_offset - 1] == L))
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
    if ((text_offset != invalid_text_offset) && (text_offset != 0) && (sl_types[text_offset - 1] == S))
    {
      text_offsets[--character_bucket_offsets[text[text_offset - 1]]] = (text_offset - 1);
      *text_offsets_it = invalid_text_offset;
    }
    --text_offsets_it;
  }
  return;
}

template <typename Iterator>
auto MoveVaildEntriesToFront
(
  Iterator begin,
  Iterator end,
  uint64_t const invalid_value
)
{
  auto it {begin};
  auto last {begin};
  while (it != end)
  {
    if (*it != invalid_value)
    {
      uint64_t value {*it}; *it = *last; *last = value;
      ++last;
    }
    ++it;
  }
  return last;
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
  uint64_t lex_rank {};
  auto prev_rule_it {std::prev(std::end(text))};
  auto prev_sl_types_it {std::prev(std::end(sl_types))};
  auto begin_offsets_it {begin_offsets_begin};
  while (begin_offsets_it != begin_offsets_end)
  {
    uint64_t begin_offset {*begin_offsets_it};
    auto rule_it {std::next(std::begin(text), begin_offset)};
    auto sl_types_it {std::next(std::begin(sl_types), begin_offset)};
    if ((*prev_rule_it == *rule_it) && (*prev_sl_types_it == *sl_types_it))
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
      if (IsLeftmostSType(prev_sl_types_it) && IsLeftmostSType(sl_types_it))
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
  int8_t step;

  DynamicGrammarTrie (int8_t const step_ = 1)
  : root {std::make_shared<Node>()},
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
  DynamicGrammarTrie const &trie,
  bool const is_rank = true
)
{
  uint64_t size {};
  for (uint64_t i {}; i != std::size(labels); ++i)
  {
    file << i << " ";
  }
  file << "\n";
  Print(file, labels);
  std::deque<std::pair<DynamicGrammarTrie::NodePointer, uint64_t>> nodes;
  nodes.emplace_back(trie.root, 0);
  while (!nodes.empty())
  {
    ++size;
    auto node {std::get<0>(nodes.back())};
    auto depth {std::get<1>(nodes.back())};
    nodes.pop_back();
    if (depth != 0)
    {
      file << depth << ":" << labels[std::get<0>(node->edge_range)];
      file << "[" << std::get<0>(node->edge_range) << "," << std::get<1>(node->edge_range) << "]";
      if (is_rank)
      {
        file << "[" << std::get<0>(node->rank_range) << "," << std::get<1>(node->rank_range) << "]";
      }
      else
      {
        file << "(" << node->count << ")";
      }
      file << "\n";
    }
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
  file << (--size) << "\n";
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
      return;
    }
    else
    {
      auto child_node {std::get<1>(*branches_it)};
      auto edge_it {std::next(rules_begin, std::get<0>(child_node->edge_range))};
      auto edge_end {std::next(rules_begin, std::get<1>(child_node->edge_range) + trie.step)};
      while ((it != last) && (edge_it != edge_end) && (*it == *edge_it))
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
      return;
    }
    else
    {
      auto child_node {std::get<1>(*branches_it)};
      auto edge_it {std::next(rules_begin, std::get<0>(child_node->edge_range))};
      auto edge_end {std::next(rules_begin, std::get<1>(child_node->edge_range) + trie.step)};
      while ((it != last) && (edge_it != edge_end) && (*it == *edge_it))
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
  sdsl::bit_vector level_order;
  sdsl::bit_vector::select_1_type level_order_select;
  sdsl::int_vector<> edge_begin_offsets;
  sdsl::int_vector<> edge_prev_end_offsets;
  sdsl::int_vector<> leftmost_ranks;
  sdsl::int_vector<> rightmost_ranks;
  sdsl::int_vector<> counts;
  int8_t step;
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
  StaticGrammarTrie const &trie,
  bool const is_rank = true
)
{
  for (uint64_t i {}; i != std::size(labels); ++i)
  {
    file << i << " ";
  }
  file << "\n";
  Print(file, labels);
  std::deque<std::pair<uint64_t, uint64_t>> offsets;
  uint64_t begin_offset {};
  uint64_t end_offset {trie.level_order_select(1)};
  for (auto offset {end_offset - 1}; offset != (begin_offset - 1); --offset)
  {
    offsets.emplace_back(offset, 1);
  }
  uint64_t size {};
  while (!offsets.empty())
  {
    ++size;
    auto offset {std::get<0>(offsets.back())};
    auto depth {std::get<1>(offsets.back())};
    offsets.pop_back();
    file << depth << ":" << labels[trie.edge_begin_offsets[offset]];
    file << "[" << trie.edge_begin_offsets[offset] << "," << trie.edge_prev_end_offsets[offset] << "]";
    if (is_rank)
    {
      file << "[" << trie.leftmost_ranks[offset] << "," << trie.rightmost_ranks[offset] << "]";
    }
    else
    {
      file << "(" << trie.counts[offset] << ")";
    }
    file << "\n";
    begin_offset = trie.level_order_select(offset + 1) - (offset + 1) + 1;
    end_offset = trie.level_order_select(offset + 2) - (offset + 2) + 1;
    for (auto offset {end_offset - 1}; offset != (begin_offset - 1); --offset)
    {
      offsets.emplace_back(offset, depth + 1);
    }
  }
  file << size << "\n";
  return;
}

template <typename File>
void SerializeStaticGrammarTrie
(
  StaticGrammarTrie const &trie,
  File &file
)
{
  sdsl::serialize(trie.level_order, file);
  sdsl::serialize(trie.level_order_select, file);
  sdsl::serialize(trie.edge_begin_offsets, file);
  sdsl::serialize(trie.edge_prev_end_offsets, file);
  sdsl::serialize(trie.leftmost_ranks, file);
  sdsl::serialize(trie.rightmost_ranks, file);
  sdsl::serialize(trie.counts, file);
  sdsl::write_member(trie.step, file);
  // {
  //   std::cout << "level_order: " << sdsl::size_in_bytes(trie.level_order) << "\n";
  //   std::cout << "level_order_select: " << sdsl::size_in_bytes(trie.level_order_select) << "\n";
  //   std::cout << "edge_begin_offsets: " << sdsl::size_in_bytes(trie.edge_begin_offsets) << "\n";
  //   std::cout << "edge_prev_end_offsets: " << sdsl::size_in_bytes(trie.edge_prev_end_offsets) << "\n";
  //   std::cout << "leftmost_ranks: " << sdsl::size_in_bytes(trie.leftmost_ranks) << "\n";
  //   std::cout << "rightmost_ranks: " << sdsl::size_in_bytes(trie.rightmost_ranks) << "\n";
  //   std::cout << "counts: " << sdsl::size_in_bytes(trie.counts) << "\n";
  // }
  return;
}

template <typename File>
void LoadStaticGrammarTrie
(
  StaticGrammarTrie &trie,
  File &file
)
{
  trie.level_order.load(file);
  trie.level_order_select.load(file);
  trie.level_order_select.set_vector(&trie.level_order);
  trie.edge_begin_offsets.load(file);
  trie.edge_prev_end_offsets.load(file);
  trie.leftmost_ranks.load(file);
  trie.rightmost_ranks.load(file);
  trie.counts.load(file);
  sdsl::read_member(trie.step, file);
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

void ConstructStaticGrammarTrie
(
  DynamicGrammarTrie const &dynamic_trie,
  StaticGrammarTrie &static_trie,
  bool const is_rank = true
)
{
  std::deque<uint8_t> level_order;
  std::deque<uint64_t> edge_begin_offsets;
  std::deque<uint64_t> edge_prev_end_offsets;
  std::deque<uint64_t> leftmost_ranks;
  std::deque<uint64_t> rightmost_ranks;
  std::deque<uint64_t> counts;
  std::deque<DynamicGrammarTrie::NodePointer> nodes;
  nodes.emplace_back(dynamic_trie.root);
  while (!nodes.empty())
  {
    auto node {nodes.front()}; nodes.pop_front();
    if (!node->branches.empty())
    {
      auto branches_it {std::begin(node->branches)};
      auto branches_end {std::end(node->branches)};
      while (branches_it != branches_end)
      {
        auto child_node {std::get<1>(*branches_it)};
        level_order.emplace_back(0);
        edge_begin_offsets.emplace_back(std::get<0>(child_node->edge_range));
        edge_prev_end_offsets.emplace_back(std::get<1>(child_node->edge_range));
        if (is_rank)
        {
          leftmost_ranks.emplace_back(std::get<0>(child_node->rank_range));
          rightmost_ranks.emplace_back(std::get<1>(child_node->rank_range));
        }
        else
        {
          counts.emplace_back(child_node->count);
        }
        nodes.emplace_back(child_node);
        ++branches_it;
      }
    }
    level_order.emplace_back(1);
  }
  {
    static_trie.level_order.resize(std::size(level_order));
    std::copy(std::begin(level_order), std::end(level_order), std::begin(static_trie.level_order));
    static_trie.level_order_select = decltype(static_trie.level_order_select){&static_trie.level_order};
  }
  InitializeAndCopy(edge_begin_offsets, static_trie.edge_begin_offsets);
  InitializeAndCopy(edge_prev_end_offsets, static_trie.edge_prev_end_offsets);
  if (is_rank)
  {
    InitializeAndCopy(leftmost_ranks, static_trie.leftmost_ranks);
    InitializeAndCopy(rightmost_ranks, static_trie.rightmost_ranks);
  }
  else
  {
    InitializeAndCopy(counts, static_trie.counts);
  }
  static_trie.step = dynamic_trie.step;
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
  sdsl::wt_int<> colex_bwt;
  // RunLengthWaveletTree colex_bwt;
};

template <typename File>
void PrintIndex
(
  File &file,
  Index &index
)
{
  Print(file, index.grammar_rules);
  PrintStaticGrammarTrie(file, index.grammar_rules, index.lex_grammar_count_trie, false);
  PrintStaticGrammarTrie(file, index.grammar_rules, index.lex_grammar_rank_trie);
  PrintStaticGrammarTrie(file, index.grammar_rules, index.colex_grammar_rank_trie);
  Print(file, index.colex_to_lex);
  Print(file, index.lex_rank_bucket_end_offsets);
  Print(file, index.colex_bwt);
  return;
}

void ConstructIndex
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
    // Print(std::cout, text);
  }
  sdsl::bit_vector sl_types;
  {
    CalculateSlTypes(text, sl_types);
    // Print(std::cout, sl_types);
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
  {
    text_offsets_boundary = MoveVaildEntriesToFront
    (
      std::begin(text_offsets),
      std::end(text_offsets),
      invalid_text_offset
    );
    // Print(std::cout, text_offsets);
  }
  auto grammar_rule_begin_offsets_begin {std::begin(text_offsets)};
  auto grammar_rule_begin_offsets_end {text_offsets_boundary};
  auto temporary_lex_text_begin {text_offsets_boundary};
  auto temporary_lex_text_end {std::end(text_offsets)};
  sdsl::int_vector<> grammar_rule_counts;
  {
    CalculateGrammarRuleCountsBeginOffsetsAndTemporaryLexText
    (
      text,
      sl_types,
      grammar_rule_counts,
      grammar_rule_begin_offsets_begin,
      grammar_rule_begin_offsets_end,
      temporary_lex_text_begin
    );
    // Print(std::cout, grammar_rule_counts);
    // Print(std::cout, grammar_rule_begin_offsets_begin, grammar_rule_begin_offsets_end);
    // Print(std::cout, temporary_lex_text_begin, temporary_lex_text_end);
  }
  sdsl::int_vector<> grammar_rule_sizes;
  {
    CalculateGrammarRuleSizes
    (
      grammar_rule_sizes,
      sl_types,
      grammar_rule_begin_offsets_begin,
      grammar_rule_begin_offsets_end,
      invalid_text_offset
    );
    // Print(std::cout, grammar_rule_sizes);
  }
  {
    CalculateGrammarRules
    (
      text,
      grammar_rule_sizes,
      index.grammar_rules,
      grammar_rule_begin_offsets_begin
    );
    // Print(std::cout, index.grammar_rules);
  }
  sdsl::int_vector<> lex_text;
  sdsl::int_vector<> lex_to_colex;
  auto grammar_ranks_size {std::size(grammar_rule_sizes) + 1};
  auto lex_text_width {sdsl::bits::hi(std::size(grammar_rule_sizes)) + 1};
  {
    CalculateLexText
    (
      lex_text,
      lex_text_width,
      temporary_lex_text_begin,
      temporary_lex_text_end,
      invalid_text_offset
    );
    // Print(std::cout, lex_text);    
  }
  {
    sdsl::util::clear(sl_types);
    sdsl::util::clear(text);
    sdsl::util::clear(text_offsets);
  }
  {
    DynamicGrammarTrie lex_grammar_count_trie;
    InsertGrammarRuleSuffixesAndCountsIntoDynamicGrammarTrie
    (
      grammar_rule_sizes,
      index.grammar_rules,
      grammar_rule_counts,
      lex_grammar_count_trie
    );
    // PrintDynamicGrammarTrie(std::cout, grammar_rules, lex_grammar_count_trie, false);
    CalculateCumulativeGrammarCount(lex_grammar_count_trie);
    // PrintDynamicGrammarTrie(std::cout, index.grammar_rules, lex_grammar_count_trie, false);
    ConstructStaticGrammarTrie
    (
      lex_grammar_count_trie,
      index.lex_grammar_count_trie,
      false
    );
    // PrintStaticGrammarTrie(std::cout, index.grammar_rules, index.lex_grammar_count_trie, false);
  }
  {
    DynamicGrammarTrie lex_grammar_rank_trie;
    DynamicGrammarTrie colex_grammar_rank_trie {-1};
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
    // PrintDynamicGrammarTrie(std::cout, index.grammar_rules, lex_grammar_rank_trie);
    ConstructStaticGrammarTrie
    (
      lex_grammar_rank_trie,
      index.lex_grammar_rank_trie
    );
    // PrintStaticGrammarTrie(std::cout, index.grammar_rules, index.lex_grammar_rank_trie);
    lex_to_colex.width(lex_text_width);
    lex_to_colex.resize(grammar_ranks_size);
    lex_to_colex[0] = 0;
    CalculateCumulativeColexRankRangesAndLexToColex(colex_grammar_rank_trie, lex_to_colex);
    // Print(std::cout, lex_to_colex);
    // PrintDynamicGrammarTrie(std::cout, index.grammar_rules, colex_grammar_rank_trie);
    ConstructStaticGrammarTrie
    (
      colex_grammar_rank_trie,
      index.colex_grammar_rank_trie
    );
    // PrintStaticGrammarTrie(std::cout, index.grammar_rules, index.colex_grammar_rank_trie);
  }
  {
    index.colex_to_lex.width(lex_text_width);
    index.colex_to_lex.resize(grammar_ranks_size);
    for (uint64_t rank {}; rank != std::size(lex_to_colex); ++rank)
    {
      index.colex_to_lex[lex_to_colex[rank]] = rank;
    }
    // Print(std::cout, index.colex_to_lex);
  }
  {
    index.lex_rank_bucket_end_offsets.width(sdsl::bits::hi(std::size(lex_text)) + 1);
    index.lex_rank_bucket_end_offsets.resize(grammar_ranks_size);
    CalculateCharacterBucketEndOffsets
    (
      lex_text,
      index.lex_rank_bucket_end_offsets
    );
    // Print(std::cout, index.lex_rank_bucket_end_offsets);
  }
  {
    CalculateColexBwt(lex_text, lex_to_colex, index.colex_bwt);
    // Print(std::cout, index.colex_bwt);
  }
  return;
}

void SerializeIndex
(
  Index &index,
  std::filesystem::path index_path
)
{
  std::fstream index_file(index_path, std::ios_base::out | std::ios_base::trunc);
  sdsl::serialize(index.grammar_rules, index_file);
  // std::cout << "grammar_rules: " << sdsl::size_in_bytes(index.grammar_rules) << "\n";
  SerializeStaticGrammarTrie(index.lex_grammar_count_trie, index_file);
  SerializeStaticGrammarTrie(index.lex_grammar_rank_trie, index_file);
  SerializeStaticGrammarTrie(index.colex_grammar_rank_trie, index_file);
  sdsl::serialize(index.colex_to_lex, index_file);
  sdsl::serialize(index.lex_rank_bucket_end_offsets, index_file);
  sdsl::serialize(index.colex_bwt, index_file);
  return;
}

void LoadIndex
(
  Index &index,
  std::filesystem::path const &index_path
)
{
  std::ifstream index_file {index_path};
  index.grammar_rules.load(index_file);
  LoadStaticGrammarTrie(index.lex_grammar_count_trie, index_file);
  LoadStaticGrammarTrie(index.lex_grammar_rank_trie, index_file);
  LoadStaticGrammarTrie(index.colex_grammar_rank_trie, index_file);
  index.colex_to_lex.load(index_file);
  index.lex_rank_bucket_end_offsets.load(index_file);
  index.colex_bwt.load(index_file);
  return;
}

template <typename Range>
constexpr bool IsEmptyRange (Range const &range)
{
  return (std::get<0>(range) > std::get<1>(range));
}

template <typename Range>
uint64_t CalculateRangeSize (Range const &range)
{
  if (std::get<0>(range) <= std::get<1>(range))
  {
    return (std::get<1>(range) - std::get<0>(range) + 1);
  }
  return 0;
}

template <typename TextIterator>
void CalculateSlFactor
(
  TextIterator const rend,
  TextIterator &rfirst,
  TextIterator &rlast
)
{
  uint64_t prev_sl_type {L};
  rfirst = rlast--;
  while ((rlast != rend) && !((prev_sl_type == S) && (*rlast > *std::next(rlast))))
  {
    if((prev_sl_type == L) && (*rlast < *std::next(rlast)))
    {
      prev_sl_type = S;
    }
    --rlast;
  }
  return;
}

template
<
  typename Labels,
  typename SlFactorIterator
>
uint64_t LookUpSlFactorCountInStaticGrammarTrie
(
  Labels const &labels,
  StaticGrammarTrie const &trie,
  SlFactorIterator it,
  SlFactorIterator last
)
{
  uint64_t begin_offset {};
  uint64_t end_offset {trie.level_order_select(1)};
  uint64_t offset {};
  while (begin_offset != end_offset)
  {
    while (begin_offset != end_offset)
    {
      offset = begin_offset + (end_offset - begin_offset) / 2;
      auto branch_character {labels[trie.edge_begin_offsets[offset]]};
      if (*it == branch_character)
      {
        break;
      }
      else if (*it < branch_character)
      {
        end_offset = offset;
      }
      else
      {
        begin_offset = offset + 1;
      }
    }
    if (begin_offset != end_offset)
    {
      auto edge_it {std::next(std::begin(labels), trie.edge_begin_offsets[offset])};
      auto edge_end {std::next(std::begin(labels), trie.edge_prev_end_offsets[offset] + trie.step)};
      while ((it != last) && (edge_it != edge_end) && (*it == *edge_it))
      {
        it += trie.step;
        edge_it += trie.step;
      }
      if (it != last)
      {
        if (edge_it != edge_end)
        {
          break;
        }
        else
        {
          begin_offset = trie.level_order_select(offset + 1) - (offset + 1) + 1;
          end_offset = trie.level_order_select(offset + 2) - (offset + 2) + 1;
        }
      }
      else
      {
        return trie.counts[offset];
      }
    }
  }
  return 0;
}

template
<
  typename Labels,
  typename SlFactorIterator
>
std::pair<uint64_t, uint64_t> LookUpSlFactorRankRangeInStaticGrammarTrie
(
  Labels const &labels,
  StaticGrammarTrie const &trie,
  SlFactorIterator it,
  SlFactorIterator last
)
{
  uint64_t begin_offset {};
  uint64_t end_offset {trie.level_order_select(1)};
  uint64_t offset {};
  while (begin_offset != end_offset)
  {
    while (begin_offset != end_offset)
    {
      offset = begin_offset + (end_offset - begin_offset) / 2;
      auto branch_character {labels[trie.edge_begin_offsets[offset]]};
      if (*it == branch_character)
      {
        break;
      }
      else if (*it < branch_character)
      {
        end_offset = offset;
      }
      else
      {
        begin_offset = offset + 1;
      }
    }
    if (begin_offset != end_offset)
    {
      auto edge_it {std::next(std::begin(labels), trie.edge_begin_offsets[offset])};
      auto edge_end {std::next(std::begin(labels), trie.edge_prev_end_offsets[offset] + trie.step)};
      while ((it != last) && (edge_it != edge_end) && (*it == *edge_it))
      {
        it += trie.step;
        edge_it += trie.step;
      }
      if (it != last)
      {
        if (edge_it != edge_end)
        {
          break;
        }
        else
        {
          begin_offset = trie.level_order_select(offset + 1) - (offset + 1) + 1;
          end_offset = trie.level_order_select(offset + 2) - (offset + 2) + 1;
        }
      }
      else
      {
        return
        {
          trie.leftmost_ranks[offset],
          trie.rightmost_ranks[offset]
        };
      }
    }
  }
  return {1, 0};
}

template
<
  typename PatternRange,
  typename PatternIterator
>
void BackwardSearchPatternPrefix
(
  Index const &index,
  PatternRange &pattern_range_l,
  PatternRange &pattern_range_s,
  PatternIterator rfirst,
  PatternIterator rlast
)
{
  auto colex_rank_range
  {
    LookUpSlFactorRankRangeInStaticGrammarTrie
    (
      index.grammar_rules,
      index.colex_grammar_rank_trie,
      rfirst,
      rlast
    )
  };
  if (!IsEmptyRange(colex_rank_range))
  {
    if (!IsEmptyRange(pattern_range_l))
    {
      pattern_range_l =
      {
        1,
        std::get<0>
        (
          index.colex_bwt.range_search_2d
          (
            std::get<0>(pattern_range_l),
            std::get<1>(pattern_range_l),
            std::get<0>(colex_rank_range),
            std::get<1>(colex_rank_range),
            false
          )
        )
      };
    }
    if (!IsEmptyRange(pattern_range_s))
    {
      pattern_range_s =
      {
        1,
        std::get<0>
        (
          index.colex_bwt.range_search_2d
          (
            std::get<0>(pattern_range_s),
            std::get<1>(pattern_range_s),
            std::get<0>(colex_rank_range),
            std::get<1>(colex_rank_range),
            false
          )
        )
      };
    }
  }
  // {
  //   Print(std::cout, rfirst, rlast, -1);
  //   std::cout
  //   << "->[" << std::get<0>(colex_rank_range) << "," << std::get<1>(colex_rank_range) << "]"
  //   << "->L:[" << std::get<0>(pattern_range_l) << "," << std::get<1>(pattern_range_l) << "]"
  //   << "->S:[" << std::get<0>(pattern_range_s) << "," << std::get<1>(pattern_range_s) << "]\n";
  // }
  return;
}

template
<
  typename PatternRange,
  typename PatternIterator
>
auto BackwardSearchPatternProperSubstring
(
  Index const &index,
  PatternRange &pattern_range_l,
  PatternRange &pattern_range_s,
  PatternIterator rfirst,
  PatternIterator rlast
)
{
  auto colex_rank_range
  {
    LookUpSlFactorRankRangeInStaticGrammarTrie
    (
      index.grammar_rules,
      index.colex_grammar_rank_trie,
      rfirst,
      rlast
    )
  };
  auto colex_rank {std::get<0>(colex_rank_range)};
  if (!IsEmptyRange(colex_rank_range))
  {
    auto begin_offset {index.lex_rank_bucket_end_offsets[index.colex_to_lex[colex_rank] - 1]};
    if (!IsEmptyRange(pattern_range_l))
    {
      pattern_range_l =
      {
        begin_offset + index.colex_bwt.rank(std::get<0>(pattern_range_l), colex_rank),
        begin_offset + index.colex_bwt.rank(std::get<1>(pattern_range_l) + 1, colex_rank) - 1
      };
    }
    if (!IsEmptyRange(pattern_range_s))
    {
      pattern_range_s =
      {
        begin_offset + index.colex_bwt.rank(std::get<0>(pattern_range_s), colex_rank),
        begin_offset + index.colex_bwt.rank(std::get<1>(pattern_range_s) + 1, colex_rank) - 1
      };
    }
  }
  else
  {
    pattern_range_l = pattern_range_s = {1, 0};
  }
  // {
  //   Print(std::cout, rfirst, rlast, -1);
  //   std::cout
  //   << "->[" << colex_rank << ":" << index.colex_to_lex[colex_rank] << "]"
  //   << "->L:[" << std::get<0>(pattern_range_l) << "," << std::get<1>(pattern_range_l) << "]"
  //   << "->S:[" << std::get<0>(pattern_range_s) << "," << std::get<1>(pattern_range_s) << "]\n";
  // }
  return rlast;
}

template <typename PatternIterator>
auto CalculatePatternSuffixS
(
  PatternIterator rbegin,
  PatternIterator rend
)
{
  auto it {rbegin};
  while ((std::prev(it) != rend) && (*std::prev(it) == *it))
  {
    --it;
  }
  if ((std::prev(it) != rend) && (*std::prev(it) < *it))
  {
    return rbegin;
  }
  return std::prev(it);
}

template
<
  typename PatternRange,
  typename PatternIterator
>
auto BackwardSearchPatternSuffix
(
  Index const &index,
  PatternRange &pattern_range_l,
  PatternRange &pattern_range_s,
  PatternIterator rbegin,
  PatternIterator rend
)
{
  auto rfirst {rbegin};
  auto rlast {CalculatePatternSuffixS(rbegin, rend)};
  if (rlast == rend)
  {
    std::get<1>(pattern_range_l) = LookUpSlFactorCountInStaticGrammarTrie
    (
      index.grammar_rules,
      index.lex_grammar_count_trie,
      std::next(rend),
      std::next(rbegin)
    );
    // {
    //   Print(std::cout, std::next(rend), std::next(rbegin));
    //   std::cout << "->L:[" << std::get<0>(pattern_range_l) << "," << std::get<1>(pattern_range_l) << "]\n";
    // }
  }
  else
  {
    if (rlast != rbegin)
    {
      auto lex_rank_range
      {
        LookUpSlFactorRankRangeInStaticGrammarTrie
        (
          index.grammar_rules,
          index.lex_grammar_rank_trie,
          std::next(rlast),
          std::next(rbegin)
        )
      };
      if (!IsEmptyRange(lex_rank_range))
      {
        pattern_range_s =
        {
          index.lex_rank_bucket_end_offsets[std::get<0>(lex_rank_range) - 1],
          (index.lex_rank_bucket_end_offsets[std::get<1>(lex_rank_range)] - 1)
        };
      }
      // {
      //   Print(std::cout, std::next(rlast), std::next(rbegin));
      //   std::cout << "->[" << std::get<0>(lex_rank_range) << "," << std::get<1>(lex_rank_range) << "]";
      //   std::cout << "->S:[" << std::get<0>(pattern_range_s) << "," << std::get<1>(pattern_range_s) << "]\n";
      // }
      CalculateSlFactor(rend, rfirst, rlast);
      if (!IsEmptyRange(pattern_range_s))
      {
        auto colex_rank_range
        {
          LookUpSlFactorRankRangeInStaticGrammarTrie
          (
            index.grammar_rules,
            index.colex_grammar_rank_trie,
            rfirst,
            rlast
          )
        };
        if (!IsEmptyRange(colex_rank_range))
        {
          if (rlast == rend)
          {
            pattern_range_s =
            {
              1,
              std::get<0>
              (
                index.colex_bwt.range_search_2d
                (
                  std::get<0>(pattern_range_s),
                  std::get<1>(pattern_range_s),
                  std::get<0>(colex_rank_range),
                  std::get<1>(colex_rank_range),
                  false
                )
              )
            };
          }
          else
          {
            auto colex_rank {std::get<0>(colex_rank_range)};
            auto begin_offset {index.lex_rank_bucket_end_offsets[index.colex_to_lex[colex_rank] - 1]};
            pattern_range_s =
            {
              begin_offset + index.colex_bwt.rank(std::get<0>(pattern_range_s), colex_rank),
              begin_offset + index.colex_bwt.rank(std::get<1>(pattern_range_s) + 1, colex_rank) - 1
            };
          }
        }
        else
        {
          pattern_range_s = {1, 0};
        }
        // {
        //   Print(std::cout, rfirst, rlast, -1);
        //   std::cout
        //   << "->[" << std::get<0>(colex_rank_range)
        //   << ":" << index.colex_to_lex[std::get<0>(colex_rank_range)]
        //   << "," << std::get<1>(colex_rank_range) << "]"
        //   << "->S:[" << std::get<0>(pattern_range_s) << "," << std::get<1>(pattern_range_s) << "]\n";
        // }
      }
    }
    else
    {
      CalculateSlFactor(rend, rfirst, rlast);
    }
    if (rlast == rend)
    {
      std::get<1>(pattern_range_l) = LookUpSlFactorCountInStaticGrammarTrie
      (
        index.grammar_rules,
        index.lex_grammar_count_trie,
        std::next(rend),
        std::next(rbegin)
      );
      // {
      //   Print(std::cout, std::next(rend), std::next(rbegin));
      //   std::cout << "->L:[" << std::get<0>(pattern_range_l) << "," << std::get<1>(pattern_range_l) << "]\n";
      // }
    }
    else
    {
      auto lex_rank_range
      {
        LookUpSlFactorRankRangeInStaticGrammarTrie
        (
          index.grammar_rules,
          index.lex_grammar_rank_trie,
          std::next(rlast),
          std::next(rbegin)
        )
      };
      if (!IsEmptyRange(lex_rank_range))
      {
        pattern_range_l =
        {
          index.lex_rank_bucket_end_offsets[std::get<0>(lex_rank_range) - 1],
          (index.lex_rank_bucket_end_offsets[std::get<1>(lex_rank_range)] - 1)
        };
      }
      // {
      //   Print(std::cout, std::next(rlast), std::next(rbegin));
      //   std::cout << "->[" << std::get<0>(lex_rank_range) << "," << std::get<1>(lex_rank_range) << "]";
      //   std::cout << "->L:[" << std::get<0>(pattern_range_l) << "," << std::get<1>(pattern_range_l) << "]\n";
      // }
    }
  }
  return rlast;
}

template <typename PatternIterator>
uint64_t Count
(
  Index const &index,
  PatternIterator begin,
  PatternIterator end
)
{
  std::pair<uint64_t, uint64_t> pattern_range_l {1, 0};
  std::pair<uint64_t, uint64_t> pattern_range_s {1, 0};
  auto rbegin {std::prev(end)};
  auto rend {std::prev(begin)};
  auto rfirst {rbegin};
  auto rlast {rbegin};
  // Print(std::cout, begin, end);
  rlast = BackwardSearchPatternSuffix(index, pattern_range_l, pattern_range_s, rbegin, rend);
  if (rlast != rend)
  {
    while (!IsEmptyRange(pattern_range_l) || !IsEmptyRange(pattern_range_s))
    {
      CalculateSlFactor(rend, rfirst, rlast);
      if (rlast != rend)
      {
        rlast = BackwardSearchPatternProperSubstring(index, pattern_range_l, pattern_range_s, rfirst, rlast);
      }
      else
      {
        BackwardSearchPatternPrefix(index, pattern_range_l, pattern_range_s, rfirst, rlast);
        break;
      }
    }
  }
  return (CalculateRangeSize(pattern_range_l) + CalculateRangeSize(pattern_range_s));
}
}
