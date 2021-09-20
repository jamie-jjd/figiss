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
    uint64_t label;
    uint64_t count;
    uint64_t lex_rank;
    std::map<uint64_t, std::shared_ptr<Node>> children;

    Node (uint64_t const label = 0)
    : label {label},
      count {},
      lex_rank {},
      children {}
    {
    }
  };

  GrammarTrie (uint8_t const length_width)
  : root_ {std::make_shared<Node>()},
    length_width_ {length_width}
  {
  }

  template <typename Iterator>
  void Insert (Iterator it, Iterator end);
  void CalculateCumulativeCountAndRank ();

  inline auto GetRoot () const
  {
    return root_;
  }

  friend std::ostream& operator<< (std::ostream& out, GrammarTrie const& grammar_trie);

private:

  std::shared_ptr<Node> root_;
  uint8_t length_width_;

};

template <typename Iterator>
void GrammarTrie::Insert (Iterator it, Iterator end)
{
  auto node {root_};
  while (it != end)
  {
    if (node->children.find(*it) == node->children.end())
    {
      node->children[*it] = std::make_shared<Node>(*it);
    }
    node = node->children[*it];
    ++(node->count);
    ++it;
  }
  node->lex_rank = 1;
  return;
}

std::ostream& operator<< (std::ostream& out, GrammarTrie const& grammar_trie)
{
  out << "depth:(byte_rank,length)(count)(lex_rank)\n";
  std::deque<std::pair<std::shared_ptr<GrammarTrie::Node>, uint64_t>> nodes;
  nodes.emplace_back(grammar_trie.root_, 0);
  uint64_t size {};
  uint64_t divisor {1ULL << grammar_trie.length_width_};
  while (!nodes.empty())
  {
    ++size;
    auto node {std::get<0>(nodes.back())};
    auto depth {std::get<1>(nodes.back())};
    nodes.pop_back();
    if (depth != 0)
    {
      out << depth << ":"
      << "(" << (node->label / divisor) << "," << (node->label % divisor) << ")"
      << "(" << node->count << ")"
      << "(" << node->lex_rank << ")\n";
    }
    for (auto it {std::rbegin(node->children)}; it != std::rend(node->children); ++it)
    {
      nodes.emplace_back(std::get<1>(*it), depth + 1);
    }
  }
  out << "|nodes|:" << (size - 1) << "\n";
  out << "length_width:" << static_cast<uint64_t>(grammar_trie.length_width_) << "\n";
  return out;
}
}
