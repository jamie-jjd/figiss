#pragma once

#include "utility.h"

namespace project
{
template <typename Labels>
struct DynamicGrammarTrie
{
  using EdgeRange = std::pair<uint64_t, uint64_t>;
  using RankRange = std::pair<uint64_t, uint64_t>;

  struct Node
  {
    std::map<typename Labels::value_type, std::shared_ptr<Node>> children;
    EdgeRange edge_range;
    RankRange rank_range;
    uint64_t count;

    Node
    (
      EdgeRange const edge_range_ = {},
      RankRange const rank_range_ = {},
      uint64_t const count_ = {}
    )
    : edge_range {edge_range_},
      rank_range {rank_range_},
      count {count_}
    {
    }
  };

  using NodePointer = std::shared_ptr<Node>;
  using LabelsIterator = typename Labels::iterator;
  using LabelsRange = std::pair<LabelsIterator, LabelsIterator>;

  NodePointer root;
  LabelsRange labels_range;
  int8_t step;

  DynamicGrammarTrie
  (
    LabelsRange const labels_range_,
    int8_t const step_ = 1
  )
  : root {std::make_shared<Node>()},
    labels_range {labels_range_},
    step {step_}
  {
  }
};

template <typename File, typename DynamicGrammarTrie>
void PrintDynamicGrammarTrie
(
  File &file,
  DynamicGrammarTrie const &trie,
  bool const is_rank = true, // false: count
  bool const is_preorder = true // false: level order
)
{
  auto labels_begin {std::get<0>(trie.labels_range)};
  auto labels_end {std::get<1>(trie.labels_range)};
  for (auto it {labels_begin}; it != labels_end; ++it)
  {
    file << std::distance(labels_begin, it) % 10;
    file << (it != std::prev(labels_end) ? " " : "\n");
  }
  Print(file, labels_begin, labels_end);
  std::deque<std::pair<typename DynamicGrammarTrie::NodePointer, uint64_t>> nodes;
  nodes.emplace_back(trie.root, 0);
  uint64_t size {};
  while (!nodes.empty())
  {
    ++size;
    auto node {std::get<0>(nodes.back())};
    auto depth {std::get<1>(nodes.back())};
    if (is_preorder)
    {
      nodes.pop_back();
    }
    else
    {
      node = std::get<0>(nodes.front());
      depth = std::get<1>(nodes.front());
      nodes.pop_front();
    }
    if (depth != 0)
    {
      file << depth << ":" << *std::next(labels_begin, std::get<0>(node->edge_range));
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
    if (is_preorder)
    {
      for (auto it {std::rbegin(node->children)}; it != std::rend(node->children); ++it)
      {
        nodes.emplace_back(std::get<1>(*it), (depth + 1));
      }
    }
    else
    {
      for (auto it {std::begin(node->children)}; it != std::end(node->children); ++it)
      {
        nodes.emplace_back(std::get<1>(*it), (depth + 1));
      }
    }
  }
  file << (--size) << "\n";
  return;
}

template <typename DynamicGrammarTrie, typename GrammarRulesIterator>
void InsertGrammarRuleInformation
(
  DynamicGrammarTrie &trie,
  GrammarRulesIterator first,
  GrammarRulesIterator last,
  uint64_t const value,
  bool const is_rank = true // false: count
)
{
  auto node {trie.root};
  auto labels_begin {std::get<0>(trie.labels_range)};
  auto it {first};
  while (true)
  {
    auto character {*it};
    auto children_it {node->children.find(character)};
    if (children_it == std::end(node->children))
    {
      node->children[character] = std::make_shared<typename DynamicGrammarTrie::Node>
      (
        typename DynamicGrammarTrie::EdgeRange
        {
          std::distance(labels_begin, it),
          (std::distance(labels_begin, last) - trie.step),
        },
        typename DynamicGrammarTrie::RankRange{value, value},
        value
      );
      return;
    }
    else
    {
      auto child {std::get<1>(*children_it)};
      auto edge_it {std::next(labels_begin, std::get<0>(child->edge_range))};
      auto edge_end {std::next(labels_begin, std::get<1>(child->edge_range) + trie.step)};
      while ((it != last) && (edge_it != edge_end) && (*it == *edge_it))
      {
        it += trie.step;
        edge_it += trie.step;
      }
      if (edge_it == edge_end)
      {
        if (it != last)
        {
          node = child;
        }
        else
        {
          if (is_rank)
          {
            child->rank_range = {value, value};
          }
          else
          {
            child->count += value;
          }
          return;
        }
      }
      else
      {
        auto edge_it_offset {std::distance(labels_begin, edge_it)};
        auto internal_node
        {
          std::make_shared<typename DynamicGrammarTrie::Node>
          (
            typename DynamicGrammarTrie::EdgeRange
            {
              std::get<0>(child->edge_range),
              (edge_it_offset - trie.step)
            }
          )
        };
        std::get<1>(*children_it) = internal_node;
        internal_node->children[*edge_it] = child;
        std::get<0>(child->edge_range) = edge_it_offset;
        if (it != last)
        {
          node = internal_node;
        }
        else
        {
          if (is_rank)
          {
            internal_node->rank_range = {value, value};
          }
          else
          {
            internal_node->count += value;
          }
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
  typename GrammarRuleCounts,
  typename DynamicGrammarTrie
>
void InsertGrammarRuleSuffixesAndCounts
(
  GrammarRuleSizes const &rule_sizes,
  GrammarRuleCounts const &rule_counts,
  DynamicGrammarTrie &lex_grammar_count_trie
)
{
  auto rules_it {std::get<0>(lex_grammar_count_trie.labels_range)};
  auto rule_counts_it {std::begin(rule_counts)};
  for (auto const &rule_size : rule_sizes)
  {
    auto rule_end {std::next(rules_it, rule_size)};
    while (rules_it != rule_end)
    {
      InsertGrammarRuleInformation
      (
        lex_grammar_count_trie,
        rules_it,
        rule_end,
        *rule_counts_it,
        false
      );
      ++rules_it;
    }
    ++rule_counts_it;
  }
  return;
}

template
<
  typename GrammarRuleSizes,
  typename DynamicGrammarTrie
>
void InsertGrammarRules
(
  GrammarRuleSizes const &rule_sizes,
  DynamicGrammarTrie &lex_grammar_rank_trie,
  DynamicGrammarTrie &colex_grammar_rank_trie
)
{
  uint64_t lex_rank {1};
  auto rules_it {std::get<0>(lex_grammar_rank_trie.labels_range)};
  for (auto const &rule_size : rule_sizes)
  {
    auto rule_end {std::next(rules_it, rule_size)};
    InsertGrammarRuleInformation
    (
      lex_grammar_rank_trie,
      rules_it,
      rule_end,
      lex_rank
    );
    InsertGrammarRuleInformation
    (
      colex_grammar_rank_trie,
      std::prev(rule_end),
      std::prev(rules_it),
      lex_rank
    );
    ++lex_rank;
    rules_it = rule_end;
  }
  return;
}

template <typename DynamicGrammarTrie>
void CalculateCumulativeGrammarCount (DynamicGrammarTrie &trie)
{
  std::deque<std::pair<typename DynamicGrammarTrie::NodePointer, bool>> nodes;
  nodes.emplace_back(trie.root, true);
  while (!nodes.empty())
  {
    auto node {std::get<0>(nodes.back())};
    auto &is_forward {std::get<1>(nodes.back())};
    if (is_forward)
    {
      is_forward = false;
      for (auto const &pair : node->children)
      {
        nodes.emplace_back(std::get<1>(pair), true);
      }
    }
    else
    {
      for (auto const &pair : node->children)
      {
        node->count += std::get<1>(pair)->count;
      }
      nodes.pop_back();
    }
  }
  return;
}

template <typename DynamicGrammarTrie>
void CalculateCumulativeLexRankRanges (DynamicGrammarTrie &trie)
{
  std::deque<std::pair<typename DynamicGrammarTrie::NodePointer, bool>> nodes;
  nodes.emplace_back(trie.root, true);
  while (!nodes.empty())
  {
    auto node {std::get<0>(nodes.back())};
    auto &is_forward {std::get<1>(nodes.back())};
    if (is_forward)
    {
      is_forward = false;
      for (auto it {std::rbegin(node->children)}; it != std::rend(node->children); ++it)
      {
        nodes.emplace_back(std::get<1>(*it), true);
      }
    }
    else
    {
      if (!node->children.empty())
      {
        auto first_child {std::get<1>(*std::begin(node->children))};
        auto last_child {std::get<1>(*std::rbegin(node->children))};
        if (std::get<0>(node->rank_range) == 0)
        {
          std::get<0>(node->rank_range) = std::get<0>(first_child->rank_range);
        }
        std::get<1>(node->rank_range) = std::get<1>(last_child->rank_range);
      }
      nodes.pop_back();
    }
  }
  return;
}

template <typename DynamicGrammarTrie, typename LexToColex>
void CalculateCumulativeColexRankRangesAndLexToColex (DynamicGrammarTrie trie, LexToColex &lex_to_colex)
{
  uint64_t colex_rank {1};
  std::deque<std::pair<typename DynamicGrammarTrie::NodePointer, bool>> nodes;
  nodes.emplace_back(trie.root, true);
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
      for (auto it {std::rbegin(node->children)}; it != std::rend(node->children); ++it)
      {
        nodes.emplace_back(std::get<1>(*it), true);
      }
    }
    else
    {
      if (!node->children.empty())
      {
        auto first_child {std::get<1>(*std::begin(node->children))};
        auto last_child {std::get<1>(*std::rbegin(node->children))};
        if (std::get<0>(node->rank_range) == 0)
        {
          std::get<0>(node->rank_range) = std::get<0>(first_child->rank_range);
        }
        std::get<1>(node->rank_range) = std::get<1>(last_child->rank_range);
      }
      nodes.pop_back();
    }
  }
  return;
}
}
