#pragma once

#include <deque>
#include <map>
#include <memory>

#include <sdsl/wavelet_trees.hpp>

#include "integer_alphabet.h"
#include "sparse_prefix_sum.h"

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

  GrammarTrie (ByteAlphabet const& byte_alphabet, IntegerAlphabet const& run_length_alphabet)
  : root_ {std::make_shared<Node>()},
    length_width_ {static_cast<uint8_t>(run_length_alphabet.GetAlphabetWidth() - byte_alphabet.GetEffectiveAlphabetWidth())},
    run_length_alphabet_ {run_length_alphabet}
  {
  }

  template <typename Iterator>
  void Insert (Iterator it, Iterator end);
  void CalculateCumulativeCountAndRank ();

  inline auto GetRoot () const
  {
    return root_;
  }

  inline uint8_t GetLengthWidth () const
  {
    return length_width_;
  }

  inline uint8_t GetAlphabetWidth () const
  {
    return run_length_alphabet_.GetAlphabetWidth();
  }

  inline uint8_t GetEffectiveAlphabetWidth () const
  {
    return run_length_alphabet_.GetEffectiveAlphabetWidth();
  }

  friend std::ostream& operator<< (std::ostream& out, GrammarTrie const& grammar_trie);

private:

  std::shared_ptr<Node> root_;
  uint8_t length_width_;
  IntegerAlphabet run_length_alphabet_;

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
    ++it;
  }
  ++(node->count);
  node->lex_rank = 1;
  return;
}

void GrammarTrie::CalculateCumulativeCountAndRank ()
{
  uint64_t lex_rank {};
  std::deque<std::pair<std::shared_ptr<GrammarTrie::Node>, bool>> nodes;
  nodes.emplace_back(root_, true);
  while (!nodes.empty())
  {
    auto node {std::get<0>(nodes.back())};
    auto& is_forward {std::get<1>(nodes.back())};
    if (is_forward)
    {
      is_forward = false;
      if (node->lex_rank)
      {
        node->lex_rank = ++lex_rank;
      }
      for (auto it {std::rbegin(node->children)}; it != std::rend(node->children); ++it)
      {
        nodes.emplace_back(std::get<1>(*it), true);
      }
    }
    else
    {
      for (auto const& pair : node->children)
      {
        node->count += std::get<1>(pair)->count;
      }
      nodes.pop_back();
    }
  }
  return;
}

std::ostream& operator<< (std::ostream& out, GrammarTrie const& grammar_trie)
{
  out << "depth:(byte_rank,length)(count)(lex_rank)\n";
  std::deque<std::pair<std::shared_ptr<GrammarTrie::Node>, uint64_t>> nodes;
  nodes.emplace_back(grammar_trie.root_, 0);
  uint64_t size {};
  uint64_t divisor {1ULL << grammar_trie.GetLengthWidth()};
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
  out << "length_width_:" << static_cast<uint64_t>(grammar_trie.GetLengthWidth()) << "\n";
  out << "alphabet_width_:" << static_cast<uint64_t>(grammar_trie.GetAlphabetWidth()) << "\n";
  out << "effective_alphabet_width_:" << static_cast<uint64_t>(grammar_trie.GetEffectiveAlphabetWidth()) << "\n";
  return out;
}
}
