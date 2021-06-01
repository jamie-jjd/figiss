#pragma once

#include "utility.h"

namespace project
{
struct DynamicGrammarTrie
{
  using EdgeRange = std::pair<uint64_t, uint64_t>;
  using RankRange = std::pair<uint64_t, uint64_t>;

  struct Node
  {
    std::map<uint64_t, std::shared_ptr<Node>> children;
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

  NodePointer root;
  int8_t step;

  DynamicGrammarTrie (int8_t const step_ = 1)
  : root {std::make_shared<Node>()},
    step {step_}
  {
  }
};

template <typename File, typename Labels>
void PrintDynamicGrammarTrie
(
  File &file,
  Labels const &labels,
  DynamicGrammarTrie const &trie,
  bool const is_rank = true, // false: count
  bool const is_preorder = true // false: level order
)
{
  for (uint64_t i {}; i != std::size(labels); ++i)
  {
    file << (i % 10);
    file << ((i != std::size(labels) - 1) ? " " : "\n");
  }
  Print(file, labels);
  std::deque<std::pair<DynamicGrammarTrie::NodePointer, uint64_t>> nodes;
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

template <typename GrammarRulesIterator>
void InsertGrammarRuleInformation
(
  DynamicGrammarTrie &trie,
  GrammarRulesIterator rules_begin,
  GrammarRulesIterator first,
  GrammarRulesIterator last,
  uint64_t const value,
  bool const is_rank = true // false: count
)
{
  auto node {trie.root};
  auto it {first};
  while (true)
  {
    auto character {*it};
    auto children_it {node->children.find(character)};
    if (children_it == std::end(node->children))
    {
      node->children[character] = std::make_shared<DynamicGrammarTrie::Node>
      (
        DynamicGrammarTrie::EdgeRange
        {
          std::distance(rules_begin, it),
          (std::distance(rules_begin, last) - trie.step),
        },
        DynamicGrammarTrie::RankRange{value, value},
        value
      );
      return;
    }
    else
    {
      auto child {std::get<1>(*children_it)};
      auto edge_it {std::next(rules_begin, std::get<0>(child->edge_range))};
      auto edge_end {std::next(rules_begin, std::get<1>(child->edge_range) + trie.step)};
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
        auto edge_it_offset {std::distance(rules_begin, edge_it)};
        auto internal_node
        {
          std::make_shared<DynamicGrammarTrie::Node>
          (
            DynamicGrammarTrie::EdgeRange
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
  typename GrammarRules,
  typename GrammarRuleCounts
>
void InsertGrammarRuleSuffixesAndCounts
(
  GrammarRuleSizes const &rule_sizes,
  GrammarRules const &rules,
  GrammarRuleCounts const &rule_counts,
  DynamicGrammarTrie &lex_grammar_count_trie
)
{
  auto rules_it {std::begin(rules)};
  auto rule_counts_it {std::begin(rule_counts)};
  for (auto const &rule_size : rule_sizes)
  {
    auto rule_end {std::next(rules_it, rule_size)};
    while (rules_it != rule_end)
    {
      InsertGrammarRuleInformation
      (
        lex_grammar_count_trie,
        std::begin(rules),
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
  typename GrammarRules
>
void InsertGrammarRules
(
  GrammarRuleSizes const &rule_sizes,
  GrammarRules const &rules,
  DynamicGrammarTrie &lex_grammar_rank_trie,
  DynamicGrammarTrie &colex_grammar_rank_trie
)
{
  uint64_t lex_rank {1};
  auto rules_it {std::begin(rules)};
  for (auto const &rule_size : rule_sizes)
  {
    auto rule_end {std::next(rules_it, rule_size)};
    InsertGrammarRuleInformation
    (
      lex_grammar_rank_trie,
      std::begin(rules),
      rules_it,
      rule_end,
      lex_rank
    );
    InsertGrammarRuleInformation
    (
      colex_grammar_rank_trie,
      std::begin(rules),
      std::prev(rule_end),
      std::prev(rules_it),
      lex_rank
    );
    ++lex_rank;
    rules_it = rule_end;
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
