#pragma once

#include <sdsl/int_vector.hpp>

#include "utility.h"

namespace project
{
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

template <typename File, typename Labels>
void PrintStaticGrammarTrie
(
  File &file,
  Labels const &labels,
  StaticGrammarTrie const &trie,
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
  std::deque<std::pair<uint64_t, uint64_t>> offsets;
  uint64_t begin_offset {};
  uint64_t end_offset {trie.level_order_select(1)};
  if (is_preorder)
  {
    for (auto offset {end_offset - 1}; offset != (begin_offset - 1); --offset)
    {
      offsets.emplace_back(offset, 1);
    }
  }
  else
  {
    for (auto offset {begin_offset}; offset != end_offset; ++offset)
    {
      offsets.emplace_back(offset, 1);
    }
  }
  uint64_t size {};
  while (!offsets.empty())
  {
    ++size;
    auto offset {std::get<0>(offsets.back())};
    auto depth {std::get<1>(offsets.back())};
    if (is_preorder)
    {
      offsets.pop_back();
    }
    else
    {
      offset = std::get<0>(offsets.front());
      depth = std::get<1>(offsets.front());
      offsets.pop_front();
    }
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
    if (is_preorder)
    {
      for (auto offset {end_offset - 1}; offset != (begin_offset - 1); --offset)
      {
        offsets.emplace_back(offset, depth + 1);
      }
    }
    else
    {
      for (auto offset {begin_offset}; offset != end_offset; ++offset)
      {
        offsets.emplace_back(offset, depth + 1);
      }
    }
  }
  file << size << "\n";
  return;
}

template <typename FromVector, typename ToVector>
void ResizeAndCopy (FromVector const &from, ToVector &to)
{
  to.width(sdsl::bits::hi(*std::max_element(std::begin(from), std::end(from))) + 1);
  to.resize(std::size(from));
  std::copy(std::begin(from), std::end(from), std::begin(to));
  return;
}

template <typename DynamicGrammarTrie>
void ConstructStaticGrammarTrie
(
  DynamicGrammarTrie const &dynamic_trie,
  StaticGrammarTrie &static_trie,
  bool const is_rank = true // false: count
)
{
  std::deque<uint8_t> level_order;
  std::deque<uint64_t> edge_begin_offsets;
  std::deque<uint64_t> edge_prev_end_offsets;
  std::deque<uint64_t> leftmost_ranks;
  std::deque<uint64_t> rightmost_ranks;
  std::deque<uint64_t> counts;
  std::deque<typename DynamicGrammarTrie::NodePointer> nodes;
  nodes.emplace_back(dynamic_trie.root);
  while (!nodes.empty())
  {
    auto node {nodes.front()}; nodes.pop_front();
    for (auto const &pair : node->children)
    {
      auto child {std::get<1>(pair)};
      edge_begin_offsets.emplace_back(std::get<0>(child->edge_range));
      edge_prev_end_offsets.emplace_back(std::get<1>(child->edge_range));
      if (is_rank)
      {
        leftmost_ranks.emplace_back(std::get<0>(child->rank_range));
        rightmost_ranks.emplace_back(std::get<1>(child->rank_range));
      }
      else
      {
        counts.emplace_back(child->count);
      }
      nodes.emplace_back(child);
      level_order.emplace_back(0);
    }
    level_order.emplace_back(1);
  }
  {
    static_trie.level_order.resize(std::size(level_order));
    std::copy(std::begin(level_order), std::end(level_order), std::begin(static_trie.level_order));
    static_trie.level_order_select = decltype(static_trie.level_order_select)(&static_trie.level_order);
  }
  ResizeAndCopy(edge_begin_offsets, static_trie.edge_begin_offsets);
  ResizeAndCopy(edge_prev_end_offsets, static_trie.edge_prev_end_offsets);
  if (is_rank)
  {
    ResizeAndCopy(leftmost_ranks, static_trie.leftmost_ranks);
    ResizeAndCopy(rightmost_ranks, static_trie.rightmost_ranks);
  }
  else
  {
    ResizeAndCopy(counts, static_trie.counts);
  }
  static_trie.step = dynamic_trie.step;
  return;
}

template <typename File, typename Node = InformationNode<std::string, uint64_t>>
uint64_t SerializeStaticGrammarTrie
(
  StaticGrammarTrie const &trie,
  File &file,
  std::shared_ptr<Node> root = nullptr
)
{
  if (root == nullptr)
  {
    sdsl::serialize(trie.level_order, file);
    sdsl::serialize(trie.level_order_select, file);
    sdsl::serialize(trie.edge_begin_offsets, file);
    sdsl::serialize(trie.edge_prev_end_offsets, file);
    sdsl::serialize(trie.leftmost_ranks, file);
    sdsl::serialize(trie.rightmost_ranks, file);
    sdsl::serialize(trie.counts, file);
    sdsl::write_member(trie.step, file);
  }
  else
  {
    {
      auto node {std::make_shared<Node>("level_order")};
      node->value = sdsl::serialize(trie.level_order, file);
      root->value += node->value;
      root->children.emplace_back(node);
    }
    {
      auto node {std::make_shared<Node>("level_order_select")};
      node->value = sdsl::serialize(trie.level_order_select, file);
      root->value += node->value;
      root->children.emplace_back(node);
    }
    {
      auto node {std::make_shared<Node>("edge_begin_offsets")};
      node->value = sdsl::serialize(trie.edge_begin_offsets, file);
      root->value += node->value;
      root->children.emplace_back(node);
    }
    {
      auto node {std::make_shared<Node>("edge_prev_end_offsets")};
      node->value = sdsl::serialize(trie.edge_prev_end_offsets, file);
      root->value += node->value;
      root->children.emplace_back(node);
    }
    {
      auto node {std::make_shared<Node>("leftmost_ranks")};
      node->value = sdsl::serialize(trie.leftmost_ranks, file);
      root->value += node->value;
      root->children.emplace_back(node);
    }
    {
      auto node {std::make_shared<Node>("rightmost_ranks")};
      node->value = sdsl::serialize(trie.rightmost_ranks, file);
      root->value += node->value;
      root->children.emplace_back(node);
    }
    {
      auto node {std::make_shared<Node>("counts")};
      node->value = sdsl::serialize(trie.counts, file);
      root->value += node->value;
      root->children.emplace_back(node);
    }
    {
      auto node {std::make_shared<Node>("step")};
      node->value = sdsl::write_member(trie.step, file);
      root->value += node->value;
      root->children.emplace_back(node);
    }
    return root->value;
  }
  return 0;
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
}
