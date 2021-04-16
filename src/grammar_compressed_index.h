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
constexpr bool IsRightmostLType (SlTypesIterator iterator)
{
  return
  (
    (*iterator == L)
    &&
    (*std::next(iterator) == S)
  );
}

template <typename SlTypesIterator>
constexpr bool IsLeftmostSType (SlTypesIterator iterator)
{
  return
  (
    (*iterator == S)
    &&
    (*std::prev(iterator) == L)
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
  auto text_iterator {std::prev(std::end(text))};
  auto text_last_iterator {std::begin(text)};
  auto sl_types_iterator {std::prev(std::end(sl_types))};
  while (text_iterator != text_last_iterator)
  {
    if
    (
      (*std::prev(text_iterator) > *text_iterator)
      ||
      (
        (*std::prev(text_iterator) == *text_iterator)
        &&
        (*sl_types_iterator == L)
      )
    )
    {
      *std::prev(sl_types_iterator) = L;
    }
    --sl_types_iterator;
    --text_iterator;
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
  auto iterator {std::prev(std::end(character_bucket_offsets))};
  auto last_iterator {std::begin(character_bucket_offsets)};
  while (iterator != last_iterator)
  {
    *iterator = *std::prev(iterator);
    --iterator;
  }
  *last_iterator = 0;
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
  auto text_iterator {std::begin(text)};
  auto text_end {std::end(text)};
  while (text_iterator != text_end)
  {
    auto text_offset {std::distance(std::begin(text), text_iterator)};
    if (IsRightmostLType(std::next(std::begin(sl_types), text_offset)))
    {
      text_offsets[character_bucket_offsets[*text_iterator]++] = text_offset;
    }
    ++text_iterator;
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
  auto text_offsets_iterator {std::next(std::begin(text_offsets))};
  auto text_offsets_last_iterator {std::end(text_offsets)};
  auto invalid_text_offset {std::size(text)};
  while (text_offsets_iterator != text_offsets_last_iterator)
  {
    auto text_offset {*text_offsets_iterator};
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
      *text_offsets_iterator = invalid_text_offset;
    }
    ++text_offsets_iterator;
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
  auto text_offsets_iterator {std::prev(std::end(text_offsets))};
  auto text_offsets_last_iterator {std::begin(text_offsets)};
  auto invalid_text_offset {std::size(text)};
  while (text_offsets_iterator != text_offsets_last_iterator)
  {
    auto text_offset {*text_offsets_iterator};
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
      *text_offsets_iterator = invalid_text_offset;
    }
    --text_offsets_iterator;
  }
  return;
}

template <typename RandomAccessIterator>
auto MoveVaildEntriesToFront
(
  RandomAccessIterator first_iterator,
  RandomAccessIterator last_iterator,
  uint64_t const invalid_value
)
{
  auto iterator {first_iterator};
  auto new_last_iterator {first_iterator};
  while (iterator != last_iterator)
  {
    if (*iterator != invalid_value)
    {
      uint64_t temporary_value {*iterator};
      *iterator = *new_last_iterator;
      *new_last_iterator = temporary_value;
      ++new_last_iterator;
    }
    ++iterator;
  }
  return new_last_iterator;
}

template
<
  typename Text,
  typename SlTypes,
  typename GrammarRuleBeginOffsetsIterator,
  typename TemporaryGrammarCompressedTextIterator,
  typename GrammarRuleCounts
>
void CalculateGrammarRuleCountsBeginOffsetsAndTemporaryGrammarCompressedText
(
  Text const &text,
  SlTypes const &sl_types,
  GrammarRuleCounts &rule_counts,
  GrammarRuleBeginOffsetsIterator begin_offsets_begin,
  GrammarRuleBeginOffsetsIterator begin_offsets_end,
  TemporaryGrammarCompressedTextIterator temporary_compressed_text_begin
)
{
  std::deque<uint64_t> counts;
  auto invalid_text_offset {std::size(text)};
  uint64_t lex_rank {0};
  auto previous_rule_iterator {std::prev(std::end(text))};
  auto previous_sl_types_iterator {std::prev(std::end(sl_types))};
  auto begin_offsets_iterator {begin_offsets_begin};
  while (begin_offsets_iterator != begin_offsets_end)
  {
    uint64_t begin_offset {*begin_offsets_iterator};
    auto rule_iterator {std::next(std::begin(text), begin_offset)};
    auto sl_types_iterator {std::next(std::begin(sl_types), begin_offset)};
    if
    (
      (*previous_rule_iterator == *rule_iterator)
      &&
      (*previous_sl_types_iterator == *sl_types_iterator)
    )
    {
      do
      {
        ++previous_rule_iterator;
        ++rule_iterator;
        ++previous_sl_types_iterator;
        ++sl_types_iterator;
      }
      while
      (
        !IsLeftmostSType(previous_sl_types_iterator)
        &&
        !IsLeftmostSType(sl_types_iterator)
        &&
        (*previous_rule_iterator == *rule_iterator)
        &&
        (*previous_sl_types_iterator == *sl_types_iterator)
      );
      if
      (
        IsLeftmostSType(previous_sl_types_iterator)
        &&
        IsLeftmostSType(sl_types_iterator)
      )
      {
        --lex_rank;
        *begin_offsets_iterator = invalid_text_offset;
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
    previous_rule_iterator = std::next(std::begin(text), begin_offset);
    previous_sl_types_iterator = std::next(std::begin(sl_types), begin_offset);
    ++begin_offsets_iterator;
  }
  rule_counts.resize(std::size(counts));
  auto rule_counts_iterator {std::begin(rule_counts)};
  for (auto const &count : counts)
  {
    *rule_counts_iterator++ = count;
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
  auto sizes_iterator {std::begin(sizes)};
  auto begin_offsets_iterator {begin_offsets_begin};
  while (begin_offsets_iterator != begin_offsets_end)
  {
    auto sl_types_first_iterator {std::next(std::begin(sl_types), *begin_offsets_iterator)};
    auto sl_types_last_iterator {std::next(sl_types_first_iterator)};
    while (!IsLeftmostSType(sl_types_last_iterator))
    {
      ++sl_types_last_iterator;
    }
    *sizes_iterator = std::distance(sl_types_first_iterator, sl_types_last_iterator);
    ++sizes_iterator;
    ++begin_offsets_iterator;
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
  auto rules_iterator {std::begin(rules)};
  auto begin_offsets_iterator {begin_offsets_begin};
  auto sizes_iterator {std::begin(sizes)};
  auto sizes_end {std::end(sizes)};
  while (sizes_iterator != sizes_end)
  {
    auto rule_iterator {std::next(std::begin(text), *begin_offsets_iterator)};
    auto rule_end {std::next(std::begin(text), *begin_offsets_iterator + *sizes_iterator)};
    while (rule_iterator != rule_end)
    {
      *rules_iterator = *rule_iterator;
      ++rules_iterator;
      ++rule_iterator;
    }
    ++begin_offsets_iterator;
    ++sizes_iterator;
  }
  return;
}

template
<
  typename GrammarCompressedText,
  typename TemporaryGrammarCompressedTextIterator
>
void CalculateGrammarCompressedText
(
  GrammarCompressedText &text,
  uint64_t const width,
  TemporaryGrammarCompressedTextIterator begin,
  TemporaryGrammarCompressedTextIterator end,
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

struct GrammarCountTrie
{
  using EdgeRange = std::pair<uint64_t, uint64_t>;

  struct Node
  {
    std::map<uint64_t, std::shared_ptr<Node>> branches;
    EdgeRange edge_range;
    uint64_t count;

    Node
    (
      EdgeRange const edge_range_ = EdgeRange{},
      uint64_t const count_ = 0
    )
    : edge_range {edge_range_},
      count {count_}
    {
    }
  };

  using NodePointer = std::shared_ptr<Node>;

  NodePointer root;
  int64_t edge_offset;

  GrammarCountTrie (int64_t const edge_offset_)
  : root {std::make_shared<Node>()},
    edge_offset {edge_offset_}
  {
  }
};

template
<
  typename File,
  typename Labels
>
void PrintGrammarCountTrie
(
  File &file,
  Labels const &labels,
  GrammarCountTrie const &trie
)
{
  uint64_t size {};
  std::deque<std::pair<GrammarCountTrie::NodePointer, uint64_t>> nodes_depth;
  nodes_depth.emplace_back(trie.root, 0);
  while (!nodes_depth.empty())
  {
    auto current_node {std::get<0>(nodes_depth.back())};
    auto depth {std::get<1>(nodes_depth.back())};
    nodes_depth.pop_back();
    ++size;
    file << depth << ":";
    auto edge_iterator {std::next(std::begin(labels), std::get<0>(current_node->edge_range))};
    auto edge_end {std::next(std::begin(labels), std::get<1>(current_node->edge_range))};
    while (edge_iterator != edge_end)
    {
      file << *edge_iterator;
      edge_iterator += trie.edge_offset;
    }
    file << "(" << current_node->count << ")\n";
    if (!current_node->branches.empty())
    {
      auto reverse_branches_iterator {std::rbegin(current_node->branches)};
      auto reverse_branches_end {std::rend(current_node->branches)};
      while (reverse_branches_iterator != reverse_branches_end)
      {
        nodes_depth.emplace_back(std::get<1>(*reverse_branches_iterator), (depth + 1));
        ++reverse_branches_iterator;
      }
    }
  }
  file << size << "\n";
  return;
}

template <typename GrammarRulesIterator>
void InsertGrammarRuleSuffixAndCountIntoGrammarCountTrie
(
  GrammarCountTrie &trie,
  GrammarRulesIterator rules_begin,
  GrammarRulesIterator rule_iterator,
  GrammarRulesIterator rule_end,
  uint64_t const count
)
{
  auto current_node {trie.root};
  auto offset {trie.edge_offset};
  while (true)
  {
    auto character {*rule_iterator};
    auto branches_iterator {current_node->branches.find(character)};
    if (branches_iterator == std::end(current_node->branches))
    {
      current_node->branches[character] = std::make_shared<GrammarCountTrie::Node>
      (
        GrammarCountTrie::EdgeRange
        {
          std::distance(rules_begin, rule_iterator),
          std::distance(rules_begin, rule_end)
        },
        count
      );
      return;
    }
    else
    {
      auto child_node {std::get<1>(*branches_iterator)};
      auto edge_iterator {std::next(rules_begin, std::get<0>(child_node->edge_range))};
      auto edge_end {std::next(rules_begin, std::get<1>(child_node->edge_range))};
      while
      (
        (rule_iterator != rule_end)
        &&
        (edge_iterator != edge_end)
        &&
        (*rule_iterator == *edge_iterator)
      )
      {
        rule_iterator += offset;
        edge_iterator += offset;
      }
      if (edge_iterator == edge_end)
      {
        if (rule_iterator != rule_end)
        {
          current_node = child_node;
        }
        else
        {
          child_node->count += count;
          return;
        }
      }
      else
      {
        auto edge_iterator_offset {std::distance(rules_begin, edge_iterator)};
        auto internal_node
        {
          std::make_shared<GrammarCountTrie::Node>
          (
            GrammarCountTrie::EdgeRange
            {
              std::get<0>(child_node->edge_range),
              edge_iterator_offset
            }
          )
        };
        std::get<1>(*branches_iterator) = internal_node;
        internal_node->branches[*edge_iterator] = child_node;
        std::get<0>(child_node->edge_range) = edge_iterator_offset;
        if (rule_iterator != rule_end)
        {
          current_node = internal_node;
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
void InsertGrammarRuleSuffixesAndCountsIntoGrammarCountTrie
(
  GrammarRuleSizes const &rule_sizes,
  GrammarRules const &rules,
  GrammarRuleCounts const &rule_counts,
  GrammarCountTrie &lex_grammar_count_trie
)
{
  auto rule_sizes_iterator {std::begin(rule_sizes)};
  auto rules_iterator {std::begin(rules)};
  auto rule_counts_iterator {std::begin(rule_counts)};
  while (rules_iterator != std::end(rules))
  {
    auto rule_begin {rules_iterator};
    auto rule_end {std::next(rule_begin, *rule_sizes_iterator)};
    while (rule_begin != rule_end)
    {
      // Print(std::cout, rule_begin, rule_end);
      InsertGrammarRuleSuffixAndCountIntoGrammarCountTrie
      (
        lex_grammar_count_trie,
        std::begin(rules),
        rule_begin,
        rule_end,
        *rule_counts_iterator
      );
      // PrintGrammarCountTrie(std::cout, rules, lex_grammar_count_trie);
      ++rule_begin;
    }
    ++rule_sizes_iterator;
    rules_iterator = rule_end;
    ++rule_counts_iterator;
  }
  return;
}

void CalculateCumulativeGrammarCount (GrammarCountTrie &trie)
{
  std::deque<std::pair<GrammarCountTrie::NodePointer, bool>> nodes_is_forward;
  nodes_is_forward.emplace_back(trie.root, true);
  while (!nodes_is_forward.empty())
  {
    auto current_node {std::get<0>(nodes_is_forward.back())};
    auto &is_forward {std::get<1>(nodes_is_forward.back())};
    if (is_forward)
    {
      is_forward = false;
      if (!current_node->branches.empty())
      {
        auto branches_iterator {std::begin(current_node->branches)};
        auto branches_end {std::end(current_node->branches)};
        while (branches_iterator != branches_end)
        {
          nodes_is_forward.emplace_back(std::get<1>(*branches_iterator), true);
          ++branches_iterator;
        }
      }
    }
    else
    {
      if (!current_node->branches.empty())
      {
        auto branches_iterator {std::begin(current_node->branches)};
        auto branches_end {std::end(current_node->branches)};
        while (branches_iterator != branches_end)
        {
          current_node->count += std::get<1>(*branches_iterator)->count;
          ++branches_iterator;
        }
      }
      nodes_is_forward.pop_back();
    }
  }
  return;
}

struct GrammarRankTrie
{
  using EdgeRange = std::pair<uint64_t, uint64_t>;
  using RankRange = std::pair<uint64_t, uint64_t>;

  struct Node
  {
    std::map<uint64_t, std::shared_ptr<Node>> branches;
    EdgeRange edge_range;
    RankRange rank_range;

    Node
    (
      EdgeRange const edge_range_ = EdgeRange{},
      RankRange const rank_range_ = RankRange{}
    )
    : edge_range {edge_range_},
      rank_range {rank_range_}
    {
    }
  };

  using NodePointer = std::shared_ptr<Node>;

  NodePointer root;
  int64_t edge_offset;

  GrammarRankTrie (int64_t const edge_offset_)
  : root {std::make_shared<Node>()},
    edge_offset {edge_offset_}
  {
  }
};

template
<
  typename File,
  typename Labels
>
void PrintGrammarRankTrie
(
  File &file,
  Labels const &labels,
  GrammarRankTrie const &trie
)
{
  uint64_t size {};
  std::deque<std::pair<GrammarRankTrie::NodePointer, uint64_t>> nodes_depth;
  nodes_depth.emplace_back(trie.root, 0);
  while (!nodes_depth.empty())
  {
    auto current_node {std::get<0>(nodes_depth.back())};
    auto depth {std::get<1>(nodes_depth.back())};
    nodes_depth.pop_back();
    ++size;
    file << depth << ":";
    auto edge_iterator {std::next(std::begin(labels), std::get<0>(current_node->edge_range))};
    auto edge_end {std::next(std::begin(labels), std::get<1>(current_node->edge_range))};
    while (edge_iterator != edge_end)
    {
      file << *edge_iterator;
      edge_iterator += trie.edge_offset;
    }
    PrintPair(file, current_node->rank_range, "\n");
    if (!current_node->branches.empty())
    {
      auto reverse_branches_iterator {std::rbegin(current_node->branches)};
      auto reverse_branches_end {std::rend(current_node->branches)};
      while (reverse_branches_iterator != reverse_branches_end)
      {
        nodes_depth.emplace_back(std::get<1>(*reverse_branches_iterator), (depth + 1));
        ++reverse_branches_iterator;
      }
    }
  }
  file << size << "\n";
  return;
}

template <typename GrammarRulesIterator>
void InsertGrammarRuleAndRankIntoGrammarRankTrie
(
  GrammarRankTrie &trie,
  GrammarRulesIterator rules_begin,
  GrammarRulesIterator rule_iterator,
  GrammarRulesIterator rule_end,
  uint64_t const rank
)
{
  auto current_node {trie.root};
  auto offset {trie.edge_offset};
  while (true)
  {
    auto character {*rule_iterator};
    auto branches_iterator {current_node->branches.find(character)};
    if (branches_iterator == std::end(current_node->branches))
    {
      current_node->branches[character] = std::make_shared<GrammarRankTrie::Node>
      (
        GrammarRankTrie::EdgeRange
        {
          std::distance(rules_begin, rule_iterator),
          std::distance(rules_begin, rule_end)
        },
        GrammarRankTrie::RankRange{rank, rank}
      );
      return;
    }
    else
    {
      auto child_node {std::get<1>(*branches_iterator)};
      auto edge_iterator {std::next(rules_begin, std::get<0>(child_node->edge_range))};
      auto edge_end {std::next(rules_begin, std::get<1>(child_node->edge_range))};
      while
      (
        (rule_iterator != rule_end)
        &&
        (edge_iterator != edge_end)
        &&
        (*rule_iterator == *edge_iterator)
      )
      {
        rule_iterator += offset;
        edge_iterator += offset;
      }
      if (edge_iterator == edge_end)
      {
        if (rule_iterator != rule_end)
        {
          current_node = child_node;
        }
        else
        {
          child_node->rank_range = {rank, rank};
          return;
        }
      }
      else
      {
        auto edge_iterator_offset {std::distance(rules_begin, edge_iterator)};
        auto internal_node
        {
          std::make_shared<GrammarRankTrie::Node>
          (
            GrammarRankTrie::EdgeRange
            {
              std::get<0>(child_node->edge_range),
              edge_iterator_offset
            }
          )
        };
        std::get<1>(*branches_iterator) = internal_node;
        internal_node->branches[*edge_iterator] = child_node;
        std::get<0>(child_node->edge_range) = edge_iterator_offset;
        if (rule_iterator != rule_end)
        {
          current_node = internal_node;
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
void InsertGrammarRulesIntoGrammarRankTries
(
  GrammarRuleSizes const &rule_sizes,
  GrammarRules const &rules,
  GrammarRankTrie &lex_grammar_rank_trie,
  GrammarRankTrie &colex_grammar_rank_trie
)
{
  uint64_t lex_rank {1};
  auto rule_sizes_iterator {std::begin(rule_sizes)};
  auto rules_iterator {std::begin(rules)};
  while (rules_iterator != std::end(rules))
  {
    auto rule_begin {rules_iterator};
    auto rule_end {std::next(rule_begin, *rule_sizes_iterator)};
    InsertGrammarRuleAndRankIntoGrammarRankTrie
    (
      lex_grammar_rank_trie,
      std::begin(rules),
      rule_begin,
      rule_end,
      lex_rank
    );
    InsertGrammarRuleAndRankIntoGrammarRankTrie
    (
      colex_grammar_rank_trie,
      std::begin(rules),
      std::prev(rule_end),
      std::prev(rule_begin),
      lex_rank
    );
    ++lex_rank;
    ++rule_sizes_iterator;
    rules_iterator = rule_end;
  }
  return;
}

void CalculateTemporaryLexGrammarTrieRankRanges (GrammarRankTrie &lex_grammar_rank_trie)
{
  std::deque<std::pair<GrammarRankTrie::NodePointer, bool>> nodes_is_forward;
  nodes_is_forward.emplace_back(lex_grammar_rank_trie.root, true);
  while (!nodes_is_forward.empty())
  {
    auto current_node {std::get<0>(nodes_is_forward.back())};
    auto &is_forward {std::get<1>(nodes_is_forward.back())};
    if (is_forward)
    {
      is_forward = false;
      if (!current_node->branches.empty())
      {
        auto reverse_branches_iterator {std::rbegin(current_node->branches)};
        auto reverse_branches_end {std::rend(current_node->branches)};
        while (reverse_branches_iterator != reverse_branches_end)
        {
          nodes_is_forward.emplace_back(std::get<1>(*reverse_branches_iterator), true);
          ++reverse_branches_iterator;
        }
      }
    }
    else
    {
      if (!current_node->branches.empty())
      {
        auto first_child_node {std::get<1>(*std::begin(current_node->branches))};
        auto last_child_node {std::get<1>(*std::rbegin(current_node->branches))};
        if (std::get<0>(current_node->rank_range) == 0)
        {
          std::get<0>(current_node->rank_range) = std::get<0>(first_child_node->rank_range);
        }
        std::get<1>(current_node->rank_range) = std::get<1>(last_child_node->rank_range);
      }
      nodes_is_forward.pop_back();
    }
  }
  return;
}

template <typename LexToColex>
void CalculateTemporaryColexGrammarTrieRankRangesAndLexToColex
(
  GrammarRankTrie colex_grammar_rank_trie,
  LexToColex &lex_to_colex
)
{
  uint64_t colex_rank {1};
  std::deque<std::pair<GrammarRankTrie::NodePointer, bool>> nodes_is_forward;
  nodes_is_forward.emplace_back(colex_grammar_rank_trie.root, true);
  while (!nodes_is_forward.empty())
  {
    auto current_node {std::get<0>(nodes_is_forward.back())};
    auto &is_forward {std::get<1>(nodes_is_forward.back())};
    if (is_forward)
    {
      is_forward = false;
      auto &leftmost_rank {std::get<0>(current_node->rank_range)};
      auto &rightmost_rank {std::get<1>(current_node->rank_range)};
      if (leftmost_rank != 0)
      {
        lex_to_colex[leftmost_rank] = colex_rank;
        leftmost_rank = rightmost_rank = colex_rank++;
      }
      if (!current_node->branches.empty())
      {
        auto reverse_branches_iterator {std::rbegin(current_node->branches)};
        auto reverse_branches_end {std::rend(current_node->branches)};
        while (reverse_branches_iterator != reverse_branches_end)
        {
          nodes_is_forward.emplace_back(std::get<1>(*reverse_branches_iterator), true);
          ++reverse_branches_iterator;
        }
      }
    }
    else
    {
      if (!current_node->branches.empty())
      {
        auto first_child_node {std::get<1>(*std::begin(current_node->branches))};
        auto last_child_node {std::get<1>(*std::rbegin(current_node->branches))};
        if (std::get<0>(current_node->rank_range) == 0)
        {
          std::get<0>(current_node->rank_range) = std::get<0>(first_child_node->rank_range);
        }
        std::get<1>(current_node->rank_range) = std::get<1>(last_child_node->rank_range);
      }
      nodes_is_forward.pop_back();
    }
  }
  return;
}

// template
// <
//   typename GrammarCompressedText,
//   typename LexToColexOrderMapping,
//   typename ColexGrammarCompressedBwt
// >
// void CalculateColexGrammarCompressedBwt
// (
//   GrammarCompressedText const &grammar_compressed_text,
//   LexToColexOrderMapping const &lex_to_colex,
//   ColexGrammarCompressedBwt &colex_grammar_compressed_bwt
// )
// {
//   sdsl::int_vector<> buffer;
//   sdsl::qsufsort::construct_sa(buffer, grammar_compressed_text);
//   auto buffer_iterator {std::begin(buffer)};
//   auto buffer_end {std::end(buffer)};
//   while (buffer_iterator != buffer_end)
//   {
//     if (*buffer_iterator != 0)
//     {
//       *buffer_iterator = lex_to_colex[grammar_compressed_text[(*buffer_iterator - 1)]];
//     }
//     ++buffer_iterator;
//   }
//   sdsl::construct_im(colex_grammar_compressed_bwt, buffer);
//   return;
// }
//
// void CalculateTrieRankRanges
// (
//   TrieNode *current_node,
//   uint64_t &rank
// )
// {
//   if (current_node->leftmost_rank != 0)
//   {
//     current_node->leftmost_rank = current_node->rightmost_rank = ++rank;
//   }
//   if (!current_node->branches.empty())
//   {
//     auto branches_begin {std::begin(current_node->branches)};
//     auto branches_iterator {branches_begin};
//     auto branches_end {std::end(current_node->branches)};
//     while (branches_iterator != branches_end)
//     {
//       CalculateTrieRankRanges(std::get<1>(*branches_iterator), rank);
//       ++branches_iterator;
//     }
//     auto first_child_node {std::get<1>(*branches_begin)};
//     auto last_child_node {std::get<1>(*std::prev(branches_end))};
//     if (current_node->leftmost_rank == 0)
//     {
//       current_node->leftmost_rank = first_child_node->leftmost_rank;
//     }
//     current_node->rightmost_rank = last_child_node->rightmost_rank;
//   }
//   return;
// }

struct Index
{
  // StaticGrammarTrie lex_grammar_rank_trie;
  // StaticGrammarTrie colex_grammar_rank_trie;
  sdsl::int_vector<> colex_to_lex;
  // sdsl::int_vector<> lex_grammar_compressed_character_bucket_end_offsets;
  // sdsl::wm_int<> colex_grammar_compressed_bwt;
};

// template <typename File>
// void PrintIndex
// (
//   File &file,
//   Index &index
// )
// {
//   file << index.grammar_rule_sizes << "\n";
//   file << index.grammar_rules << "\n";
//   PrintGrammarTrie(file, std::begin(index.grammar_rules), index.lex_grammar_rank_trie_root, 1);
//   PrintGrammarTrie(file, std::begin(index.grammar_rules), index.colex_grammar_rank_trie_root, -1);
//   file << index.colex_to_lex << "\n";
//   file << index.lex_grammar_compressed_character_bucket_end_offsets << "\n";
//   file << index.colex_grammar_compressed_bwt << "\n";
//   return;
// }

void Construct
(
  // Index &index,
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
    Print(std::cout, std::begin(text), std::end(text));
  }
  sdsl::bit_vector sl_types;
  {
    CalculateSlTypes(text, sl_types);
    Print(std::cout, std::begin(sl_types), std::end(sl_types));
  }
  auto invalid_text_offset {std::size(text)};
  auto text_size_width {sdsl::bits::hi(std::size(text)) + 1};
  sdsl::int_vector<> text_offsets;
  {
    text_offsets.width(text_size_width);
    text_offsets.resize(std::size(text));
    sdsl::util::set_to_value(text_offsets, invalid_text_offset);
    // Print(std::cout, std::begin(text_offsets), std::end(text_offsets));
  }
  sdsl::int_vector<> character_bucket_offsets;
  {
    character_bucket_offsets.width(text_size_width);
    character_bucket_offsets.resize(256);
  }
  {
    CalculateCharacterBucketBeginOffsets(text, character_bucket_offsets);
    // Print(std::cout, std::begin(character_bucket_offsets), std::end(character_bucket_offsets));
    BucketSortRightmostLTypeCharacters(text, sl_types, character_bucket_offsets, text_offsets);
    // Print(std::cout, std::begin(text_offsets), std::end(text_offsets));
    InduceSortLTypeCharacters(text, sl_types, character_bucket_offsets, text_offsets);
    // Print(std::cout, std::begin(text_offsets), std::end(text_offsets));
    CalculateCharacterBucketEndOffsets(text, character_bucket_offsets);
    // Print(std::cout, std::begin(character_bucket_offsets), std::end(character_bucket_offsets));
    InduceSortSTypeCharacters(text, sl_types, character_bucket_offsets, text_offsets);
    // Print(std::cout, std::begin(text_offsets), std::end(text_offsets));
  }
  auto text_offsets_boundary {std::begin(text_offsets)};
  text_offsets_boundary = MoveVaildEntriesToFront
  (
    std::begin(text_offsets),
    std::end(text_offsets),
    invalid_text_offset
  );
  // Print(std::cout, std::begin(text_offsets), std::end(text_offsets));
  auto grammar_rule_begin_offsets_begin {std::begin(text_offsets)};
  auto grammar_rule_begin_offsets_end {text_offsets_boundary};
  auto temporary_grammar_compressed_text_begin {text_offsets_boundary};
  auto temporary_grammar_compressed_text_end {std::end(text_offsets)};
  sdsl::int_vector<> grammar_rule_counts;
  CalculateGrammarRuleCountsBeginOffsetsAndTemporaryGrammarCompressedText
  (
    text,
    sl_types,
    grammar_rule_counts,
    grammar_rule_begin_offsets_begin,
    grammar_rule_begin_offsets_end,
    temporary_grammar_compressed_text_begin
  );
  Print(std::cout, std::begin(grammar_rule_counts), std::end(grammar_rule_counts));
  // Print(std::cout, grammar_rule_begin_offsets_begin, grammar_rule_begin_offsets_end);
  // Print(std::cout, temporary_grammar_compressed_text_begin, temporary_grammar_compressed_text_end);
  sdsl::int_vector<> grammar_rule_sizes;
  CalculateGrammarRuleSizes
  (
    grammar_rule_sizes,
    sl_types,
    grammar_rule_begin_offsets_begin,
    grammar_rule_begin_offsets_end,
    invalid_text_offset
  );
  Print(std::cout, std::begin(grammar_rule_sizes), std::end(grammar_rule_sizes));
  sdsl::int_vector<8> grammar_rules;
  CalculateGrammarRules
  (
    text,
    grammar_rule_sizes,
    grammar_rules,
    grammar_rule_begin_offsets_begin
  );
  Print(std::cout, std::begin(grammar_rules), std::end(grammar_rules));
  // auto grammar_compressed_alphabet_size {std::size(grammar_rule_sizes) + 1};
  auto grammar_compressed_text_width {sdsl::bits::hi(std::size(grammar_rule_sizes)) + 1};
  sdsl::int_vector<> grammar_compressed_text;
  CalculateGrammarCompressedText
  (
    grammar_compressed_text,
    grammar_compressed_text_width,
    temporary_grammar_compressed_text_begin,
    temporary_grammar_compressed_text_end,
    invalid_text_offset
  );
  Print(std::cout, std::begin(grammar_compressed_text), std::end(grammar_compressed_text));
  {
    sdsl::util::clear(sl_types);
    sdsl::util::clear(text);
    sdsl::util::clear(text_offsets);
  }
  GrammarCountTrie lex_grammar_count_trie {1};
  GrammarRankTrie lex_grammar_rank_trie {1};
  GrammarRankTrie colex_grammar_rank_trie {-1};
  {
    InsertGrammarRuleSuffixesAndCountsIntoGrammarCountTrie
    (
      grammar_rule_sizes,
      grammar_rules,
      grammar_rule_counts,
      lex_grammar_count_trie
    );
    // PrintGrammarCountTrie(std::cout, grammar_rules, lex_grammar_count_trie);
    CalculateCumulativeGrammarCount(lex_grammar_count_trie);
    InsertGrammarRulesIntoGrammarRankTries
    // PrintGrammarCountTrie(std::cout, grammar_rules, lex_grammar_count_trie);
    (
      grammar_rule_sizes,
      grammar_rules,
      lex_grammar_rank_trie,
      colex_grammar_rank_trie
    );
    PrintGrammarRankTrie(std::cout, grammar_rules, lex_grammar_rank_trie);
    PrintGrammarRankTrie(std::cout, grammar_rules, colex_grammar_rank_trie);
    // CalculateTemporaryLexGrammarTrieRankRanges(lex_grammar_rank_trie);
    // PrintGrammarRankTrie(std::cout, grammar_rules, lex_grammar_rank_trie);
    // sdsl::int_vector<> lex_to_colex;
    // lex_to_colex.width(grammar_compressed_text_width);
    // lex_to_colex.resize(grammar_compressed_alphabet_size);
    // lex_to_colex[0] = 0;
    // CalculateTemporaryColexGrammarTrieRankRangesAndLexToColex(colex_grammar_rank_trie, lex_to_colex);
    // // Print(std::cout, std::begin(lex_to_colex), std::end(lex_to_colex));
    // PrintGrammarRankTrie(std::cout, grammar_rules, colex_grammar_rank_trie);
    // index.colex_to_lex.width(grammar_compressed_text_width);
    // index.colex_to_lex.resize(grammar_compressed_alphabet_size);
    // for (uint64_t rank {}; rank != std::size(lex_to_colex); ++rank)
    // {
    //   index.colex_to_lex[lex_to_colex[rank]] = rank;
    // }
    // Print(std::cout, std::begin(index.colex_to_lex), std::end(index.colex_to_lex));
    // InsertGrammarRulesSuffixesIntoTemporaryLexGrammarTrie
    // (
    //   grammar_rule_counts,
    //   grammar_rule_sizes,
    //   grammar_rules,
    //   lex_grammar_rank_trie
    // );
    // // PrintGrammarRankTrie(std::cout, grammar_rules, lex_grammar_rank_trie);
    // CalculateGrammarRankTrieCount(lex_grammar_rank_trie);
    // PrintGrammarRankTrie(std::cout, grammar_rules, lex_grammar_rank_trie);
  }
  // index.lex_grammar_compressed_character_bucket_end_offsets.width(sdsl::bits::hi(grammar_compressed_text_size) + 1);
  // index.lex_grammar_compressed_character_bucket_end_offsets.resize(grammar_compressed_alphabet_size);
  // CalculateCharacterBucketEndOffsets
  // (
  //   grammar_compressed_text,
  //   index.lex_grammar_compressed_character_bucket_end_offsets
  // );
  // CalculateColexGrammarCompressedBwt
  // (
  //   grammar_compressed_text,
  //   lex_to_colex,
  //   index.colex_grammar_compressed_bwt
  // );
  return;
}

// void Serialize
// (
//   Index &index,
//   std::filesystem::path index_path
// )
// {
//   std::ofstream index_file {index_path};
//   // todo: serialize grammar trie
//   sdsl::serialize(index.colex_to_lex, index_file);
//   sdsl::serialize(index.lex_grammar_compressed_character_bucket_end_offsets, index_file);
//   sdsl::serialize(index.colex_grammar_compressed_bwt, index_file);
//   return;
// }
//
// void Load
// (
//   Index &index,
//   std::filesystem::path const &index_path
// )
// {
//   sdsl::util::clear(index.colex_to_lex);
//   sdsl::util::clear(index.lex_grammar_compressed_character_bucket_end_offsets);
//   sdsl::util::clear(index.colex_grammar_compressed_bwt);
//   std::ifstream index_file {index_path};
//   // todo: load grammar trie
//   index.colex_to_lex.load(index_file);
//   index.lex_grammar_compressed_character_bucket_end_offsets.load(index_file);
//   index.colex_grammar_compressed_bwt.load(index_file);
//   return;
// }
//
// template <typename TextIterator>
// void CalculateSlFactor
// (
//   TextIterator &reverse_first_iterator,
//   TextIterator &reverse_last_iterator,
//   TextIterator reverse_end
// )
// {
//   uint64_t previous_sl_type {L};
//   reverse_first_iterator = reverse_last_iterator--;
//   while
//   (
//     (reverse_last_iterator != reverse_end)
//     &&
//     ! (
//         (previous_sl_type == S)
//         &&
//         (*reverse_last_iterator > *std::next(reverse_last_iterator))
//       )
//   )
//   {
//     if
//     (
//       (previous_sl_type == L)
//       &&
//       (*reverse_last_iterator < *std::next(reverse_last_iterator))
//     )
//     {
//       previous_sl_type = S;
//     }
//     --reverse_last_iterator;
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
//   TrieNode *current_node,
//   SlFactorIterator sl_factor_iterator,
//   SlFactorIterator sl_factor_end,
//   int64_t const offset,
//   uint64_t &leftmost_rank,
//   uint64_t &rightmost_rank
// )
// {
//   leftmost_rank = rightmost_rank = 0;
//   auto branches_iterator {std::end(current_node->branches)};
//   while
//   (
//     (branches_iterator = current_node->branches.find(*sl_factor_iterator))
//     !=
//     std::end(current_node->branches)
//   )
//   {
//     current_node = std::get<1>(*branches_iterator);
//     auto edge_iterator {std::next(grammar_rules_begin, current_node->edge_begin_offset)};
//     auto edge_end {std::next(grammar_rules_begin, current_node->edge_end_offset)};
//     while
//     (
//       (sl_factor_iterator != sl_factor_end)
//       &&
//       (edge_iterator != edge_end)
//       &&
//       (*sl_factor_iterator == *edge_iterator)
//     )
//     {
//       sl_factor_iterator += offset;
//       edge_iterator += offset;
//     }
//     if (sl_factor_iterator != sl_factor_end)
//     {
//       if (edge_iterator != edge_end)
//       {
//         return;
//       }
//     }
//     else
//     {
//       rightmost_rank = current_node->rightmost_rank;
//       leftmost_rank = current_node->leftmost_rank;
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
//   PatternIterator reverse_pattern_first_iterator,
//   PatternIterator reverse_pattern_last_iterator
// )
// {
//   uint64_t leftmost_colex_rank {0};
//   uint64_t rightmost_colex_rank {0};
//   LookupGrammarRule
//   (
//     std::begin(index.grammar_rules),
//     index.colex_grammar_rank_trie_root,
//     reverse_pattern_first_iterator,
//     reverse_pattern_last_iterator,
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
//         index.colex_grammar_compressed_bwt.range_search_2d
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
//         index.colex_grammar_compressed_bwt.range_search_2d
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
//   //   auto iterator {std::next(reverse_pattern_last_iterator)};
//   //   auto end {std::next(reverse_pattern_first_iterator)};
//   //   while ( iterator != end)
//   //   {
//   //     std::cout << *iterator;
//   //     ++iterator;
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
//   PatternIterator reverse_pattern_first_iterator,
//   PatternIterator reverse_pattern_last_iterator
// )
// {
//   uint64_t colex_rank {0};
//   LookupGrammarRule
//   (
//     std::begin(index.grammar_rules),
//     index.colex_grammar_rank_trie_root,
//     reverse_pattern_first_iterator,
//     reverse_pattern_last_iterator,
//     -1,
//     colex_rank,
//     colex_rank
//   );
//   if (colex_rank != 0)
//   {
//     auto character_begin_offset
//     {
//       index.lex_grammar_compressed_character_bucket_end_offsets
//       [
//         index.colex_to_lex[colex_rank] - 1
//       ]
//     };
//     if (pattern_range_begin_offset_L != pattern_range_end_offset_L)
//     {
//       pattern_range_begin_offset_L =
//       (
//         character_begin_offset
//         + index.colex_grammar_compressed_bwt.rank
//         (
//           pattern_range_begin_offset_L,
//           colex_rank
//         )
//       );
//       pattern_range_end_offset_L =
//       (
//         character_begin_offset
//         + index.colex_grammar_compressed_bwt.rank
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
//         + index.colex_grammar_compressed_bwt.rank
//         (
//           pattern_range_begin_offset_S,
//           colex_rank
//         )
//       );
//       pattern_range_end_offset_S =
//       (
//         character_begin_offset
//         + index.colex_grammar_compressed_bwt.rank
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
//   //   auto iterator {std::next(reverse_pattern_last_iterator)};
//   //   auto end {std::next(reverse_pattern_first_iterator)};
//   //   while (iterator != end)
//   //   {
//   //     std::cout << *iterator;
//   //     ++iterator;
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
//   auto iterator {reverse_begin};
//   while
//   (
//     (std::prev(iterator) != reverse_end)
//     &&
//     (*std::prev(iterator) == *iterator)
//   )
//   {
//     --iterator;
//   }
//   if
//   (
//     (std::prev(iterator) != reverse_end)
//     &&
//     (*std::prev(iterator) < *iterator)
//   )
//   {
//     return reverse_begin;
//   }
//   return std::prev(iterator);
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
//   auto reverse_pattern_first_iterator {reverse_pattern_begin};
//   auto reverse_pattern_last_iterator
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
//     (reverse_pattern_last_iterator != reverse_pattern_begin)
//     &&
//     (reverse_pattern_last_iterator != reverse_pattern_end)
//   )
//   {
//     LookupGrammarRule
//     (
//       std::begin(index.grammar_rules),
//       index.lex_grammar_rank_trie_root,
//       std::next(reverse_pattern_last_iterator),
//       std::next(reverse_pattern_begin),
//       1,
//       leftmost_lex_rank,
//       rightmost_lex_rank
//     );
//     if (leftmost_lex_rank != 0)
//     {
//       pattern_range_begin_offset_S = index.lex_grammar_compressed_character_bucket_end_offsets[leftmost_lex_rank - 1];
//       pattern_range_end_offset_S = index.lex_grammar_compressed_character_bucket_end_offsets[rightmost_lex_rank];
//     }
//     // {
//     //   auto iterator {std::next(reverse_pattern_last_iterator)};
//     //   auto end {std::next(reverse_pattern_begin)};
//     //   while (iterator != end)
//     //   {
//     //     std::cout << *iterator;
//     //     ++iterator;
//     //   }
//     //   std::cout << "->(" << leftmost_lex_rank << "," << rightmost_lex_rank << ")";
//     //   std::cout << "->S:(" << pattern_range_begin_offset_S << "," << pattern_range_end_offset_S << ")\n";
//     // }
//     CalculateSlFactor(reverse_pattern_first_iterator, reverse_pattern_last_iterator, reverse_pattern_end);
//     if (pattern_range_begin_offset_S != pattern_range_end_offset_S)
//     {
//       uint64_t leftmost_colex_rank {0};
//       uint64_t rightmost_colex_rank {0};
//       LookupGrammarRule
//       (
//         std::begin(index.grammar_rules),
//         index.colex_grammar_rank_trie_root,
//         reverse_pattern_first_iterator,
//         reverse_pattern_last_iterator,
//         -1,
//         leftmost_colex_rank,
//         rightmost_colex_rank
//       );
//       if (leftmost_colex_rank != 0)
//       {
//         if (reverse_pattern_last_iterator != reverse_pattern_end)
//         {
//           auto character_begin_offset
//           {
//             index.lex_grammar_compressed_character_bucket_end_offsets
//             [
//               index.colex_to_lex[leftmost_colex_rank] - 1
//             ]
//           };
//           pattern_range_begin_offset_S =
//           (
//             character_begin_offset
//             + index.colex_grammar_compressed_bwt.rank
//             (
//               pattern_range_begin_offset_S,
//               leftmost_colex_rank
//             )
//           );
//           pattern_range_end_offset_S =
//           (
//             character_begin_offset
//             + index.colex_grammar_compressed_bwt.rank
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
//             index.colex_grammar_compressed_bwt.range_search_2d
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
//       //   auto iterator {std::next(reverse_pattern_last_iterator)};
//       //   auto end {std::next(reverse_pattern_first_iterator)};
//       //   while (iterator != end)
//       //   {
//       //     std::cout << *iterator;
//       //     ++iterator;
//       //   }
//       //   std::cout
//       //   << "->(" << leftmost_colex_rank << ":" << index.colex_to_lex[leftmost_colex_rank]
//       //   << "," << rightmost_colex_rank << ")"
//       //   << "->S:(" << pattern_range_begin_offset_S << "," << pattern_range_end_offset_S << ")\n";
//       // }
//     }
//   }
//   else if (reverse_pattern_last_iterator == reverse_pattern_begin)
//   {
//     CalculateSlFactor
//     (
//       reverse_pattern_first_iterator,
//       reverse_pattern_last_iterator,
//       reverse_pattern_end
//     );
//   }
//   LookupGrammarRule
//   (
//     std::begin(index.grammar_rules),
//     index.lex_grammar_rank_trie_root,
//     std::next(reverse_pattern_last_iterator),
//     std::next(reverse_pattern_begin),
//     1,
//     leftmost_lex_rank,
//     rightmost_lex_rank
//   );
//   if (leftmost_lex_rank != 0)
//   {
//     pattern_range_begin_offset_L = index.lex_grammar_compressed_character_bucket_end_offsets[leftmost_lex_rank - 1];
//     pattern_range_end_offset_L = index.lex_grammar_compressed_character_bucket_end_offsets[rightmost_lex_rank];
//   }
//   // {
//   //   auto iterator {std::next(reverse_pattern_last_iterator)};
//   //   auto end {std::next(reverse_pattern_begin)};
//   //   while (iterator != end)
//   //   {
//   //     std::cout << *iterator;
//   //     ++iterator;
//   //   }
//   //   std::cout << "->(" << leftmost_lex_rank << "," << rightmost_lex_rank << ")";
//   //   std::cout << "->L:(" << pattern_range_begin_offset_L << "," << pattern_range_end_offset_L << ")\n";
//   // }
//   return reverse_pattern_last_iterator;
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
//   auto reverse_pattern_first_iterator {reverse_pattern_begin};
//   auto reverse_pattern_last_iterator {reverse_pattern_begin};
//   // {
//   //   for (auto iterator {pattern_begin}; iterator != pattern_end; ++iterator)
//   //   {
//   //     std::cout << *iterator;
//   //   }
//   //   std::cout << "\n";
//   // }
//   reverse_pattern_last_iterator = BackwardSearchPatternSuffix
//   (
//     index,
//     pattern_range_begin_offset_L,
//     pattern_range_end_offset_L,
//     pattern_range_begin_offset_S,
//     pattern_range_end_offset_S,
//     reverse_pattern_begin,
//     reverse_pattern_end
//   );
//   if (reverse_pattern_last_iterator != reverse_pattern_end)
//   {
//     while
//     (
//       (pattern_range_begin_offset_L != pattern_range_end_offset_L)
//       ||
//       (pattern_range_begin_offset_S != pattern_range_end_offset_S)
//     )
//     {
//       CalculateSlFactor(reverse_pattern_first_iterator, reverse_pattern_last_iterator, reverse_pattern_end);
//       if (reverse_pattern_last_iterator!= reverse_pattern_end)
//       {
//         BackwardSearchExactSlFactor
//         (
//           index,
//           pattern_range_begin_offset_L,
//           pattern_range_end_offset_L,
//           pattern_range_begin_offset_S,
//           pattern_range_end_offset_S,
//           reverse_pattern_first_iterator,
//           reverse_pattern_last_iterator
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
//           reverse_pattern_first_iterator,
//           reverse_pattern_last_iterator
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
//   auto reverse_first_iterator {std::prev(std::end(text))};
//   auto reverse_last_iterator {std::prev(std::end(text))};
//   auto reverse_end {std::prev(std::begin(text))};
//   while (reverse_last_iterator != reverse_end)
//   {
//     CalculateSlFactor(reverse_first_iterator, reverse_last_iterator, reverse_end);
//     auto size {static_cast<uint64_t>(std::distance(reverse_last_iterator, reverse_first_iterator))};
//     if (max_size < size)
//     {
//       max_size = size;
//     }
//   }
//   return max_size;
// }
}
