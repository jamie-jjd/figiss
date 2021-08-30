#pragma once

#include <deque>
#include <map>
#include <memory>

namespace figiss
{
class GrammarTrie
{
public:

  struct Node
  {
    std::map<uint64_t, std::shared_ptr<Node>> children;
    std::pair<uint64_t, uint64_t> rank_range;
    uint64_t count;
  };

  GrammarTrie (uint64_t const length_width)
  : root_ {std::make_shared<Node>()},
    length_width_ {length_width}
  {
  }

  template <typename Iterator>
  void Insert (Iterator it, Iterator last);
  void CalculateCumulativeCountAndRankRange ();

  inline auto GetRoot () const
  {
    return root_;
  }

  friend std::ostream& operator<< (std::ostream &out, GrammarTrie const &grammar_trie);

private:

  std::shared_ptr<Node> root_;
  uint64_t length_width_;

};

template <typename Iterator>
void GrammarTrie::Insert (Iterator it, Iterator end)
{
  auto node {root_};
  while (it != end)
  {
    if (node->children.find(*it) == node->children.end())
    {
      node->children[*it] = std::make_shared<Node>();
    }
    node = node->children[*it];
    ++it;
  }
  ++(node->count);
  std::get<0>(node->rank_range) = 1;
  return;
}

void GrammarTrie::CalculateCumulativeCountAndRankRange ()
{
  uint64_t rank {};
  std::deque<std::pair<std::shared_ptr<GrammarTrie::Node>, bool>> nodes;
  nodes.emplace_back(root_, true);
  while (!nodes.empty())
  {
    auto node {std::get<0>(nodes.back())};
    auto &is_forward {std::get<1>(nodes.back())};
    if (is_forward)
    {
      is_forward = false;
      if (std::get<0>(node->rank_range))
      {
        std::get<0>(node->rank_range) = std::get<1>(node->rank_range) = ++rank;
      }
      for (auto it {std::rbegin(node->children)}; it != std::rend(node->children); ++it)
      {
        nodes.emplace_back(std::get<1>(*it), true);
      }
    }
    else
    {
      for (auto const &pair : node->children)
      {
        node->count += std::get<1>(pair)->count;
      }
      auto first_child {std::get<1>(*std::begin(node->children))};
      auto last_child {std::get<1>(*std::rend(node->children))};
      if (std::get<0>(node->rank_range) == 0)
      {
        std::get<0>(node->rank_range) = std::get<0>(first_child->rank_range);
      }
      std::get<1>(node->rank_range) = std::get<1>(last_child->rank_range);
      nodes.pop_back();
    }
  }
  return;
}

std::ostream& operator<< (std::ostream &out, GrammarTrie const &grammar_trie)
{
  out << "depth:(byte_rank,length)[leftmost_rank:rightmost_rank](count)\n";
  std::deque<std::tuple<uint64_t, std::shared_ptr<GrammarTrie::Node>, uint64_t>> nodes;
  nodes.emplace_back(0, grammar_trie.root_, 0);
  uint64_t size {};
  uint64_t divisor {1ULL << grammar_trie.length_width_};
  while (!nodes.empty())
  {
    ++size;
    auto label {std::get<0>(nodes.back())};
    auto node {std::get<1>(nodes.back())};
    auto depth {std::get<2>(nodes.back())};
    nodes.pop_back();
    if (depth != 0)
    {
      out << depth << ":"
      << "(" << (label / divisor) << "," << (label % divisor) << ")"
      << "[" << std::get<0>(node->rank_range) << "," << std::get<1>(node->rank_range) << "]"
      << "(" << node->count << ")\n";
    }
    for (auto it {std::rbegin(node->children)}; it != std::rend(node->children); ++it)
    {
      nodes.emplace_back(std::get<0>(*it), std::get<1>(*it), depth + 1);
    }
  }
  out << "|nodes|:" << (size - 1) << "\n";
  return out;
}
}
