#pragma once

#include <deque>
#include <map>
#include <memory>

#include <sdsl/wavelet_trees.hpp>

#include "symbol_bucket_offsets.h"

namespace figiss
{
class GrammarTrie
{
public:

  struct Node
  {
    uint64_t label;
    uint64_t count;
    std::pair<uint64_t, uint64_t> rank_range;
    std::map<uint64_t, std::shared_ptr<Node>> children;

    Node (uint64_t const label = 0): label {label} {}
  };

  GrammarTrie (uint8_t const effective_alphabet_width, uint64_t const length_width)
  : root_ {std::make_shared<Node>()},
    effective_alphabet_width_ {effective_alphabet_width},
    length_width_ {length_width}
  {
  }

  template <typename Iterator>
  void Insert (Iterator it, Iterator end);
  void CalculateCumulativeCountAndRankRange ();

  inline auto GetRoot () const
  {
    return root_;
  }

  inline uint8_t GetEffectiveAlphabetWidth () const
  {
    return effective_alphabet_width_;
  }

  inline uint64_t GetLengthWidth () const
  {
    return length_width_;
  }

  friend std::ostream& operator<< (std::ostream& out, GrammarTrie const& grammar_trie);

private:

  std::shared_ptr<Node> root_;
  uint8_t effective_alphabet_width_;
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
      node->children[*it] = std::make_shared<Node>(*it);
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
    auto& is_forward {std::get<1>(nodes.back())};
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
      if (!node->children.empty())
      {
        for (auto const& pair : node->children)
        {
          node->count += std::get<1>(pair)->count;
        }
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

std::ostream& operator<< (std::ostream& out, GrammarTrie const& grammar_trie)
{
  out << "depth:(byte_rank,length)(count)[leftmost_rank:rightmost_rank]\n";
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
      << "[" << std::get<0>(node->rank_range) << "," << std::get<1>(node->rank_range) << "]\n";
    }
    for (auto it {std::rbegin(node->children)}; it != std::rend(node->children); ++it)
    {
      nodes.emplace_back(std::get<1>(*it), depth + 1);
    }
  }
  out << "|nodes|:" << (size - 1) << "\n";
  return out;
}

class GrammarXbwtTrie
{
public:

  static constexpr uint64_t kSpecialSymbol {1};

  GrammarXbwtTrie () = default;
  GrammarXbwtTrie (GrammarTrie const& grammar_trie);
  GrammarXbwtTrie& operator= (GrammarXbwtTrie&&);

  // uint64_t Serialize
  // (
  //   std::ostream& out,
  //   std::shared_ptr<SpaceNode> parent = nullptr,
  //   std::string const name = ""
  // );
  // void Load (std::istream& in);

  // template <typename Iterator>
  // uint64_t Count (Iterator it, Iterator end);

  // friend std::ostream& operator<< (std::ostream& out, GrammarXbwtTrie const& grammar_xbwt_trie);

private:

  void CalculateReverseGrammarRulesAndGrammarTrieNodePointers
  (
    GrammarTrie const& grammar_trie,
    sdsl::int_vector<>& reverse_grammar_rules,
    std::deque<std::shared_ptr<GrammarTrie::Node>>& grammar_trie_node_pointers
  );

  sdsl::bit_vector trie_bits_;
  sdsl::bit_vector::rank_1_type trie_rank_1_;
  sdsl::bit_vector::select_1_type trie_select_1_;
  sdsl::wt_rlmn
  <
    sdsl::sd_vector<>,
    typename sdsl::sd_vector<>::rank_1_type,
    typename sdsl::sd_vector<>::select_1_type,
    sdsl::wt_ap<>
  >
  bwt_;
  SymbolBucketOffsets symbol_bucket_offsets_;
  sdsl::int_vector<> counts_;
  sdsl::int_vector<> lex_rank_range_;

};

GrammarXbwtTrie::GrammarXbwtTrie (GrammarTrie const& grammar_trie)
{
  grammar_trie.GetLengthWidth();
  sdsl::int_vector<> reverse_grammar_rules;
  std::deque<std::shared_ptr<GrammarTrie::Node>> grammar_trie_node_pointers;
  CalculateReverseGrammarRulesAndGrammarTrieNodePointers
  (
    grammar_trie,
    reverse_grammar_rules,
    grammar_trie_node_pointers
  );
}

GrammarXbwtTrie& GrammarXbwtTrie::operator= (GrammarXbwtTrie&& grammar_xbwt_trie)
{
  if (this != &grammar_xbwt_trie)
  {
    trie_bits_ = std::move(grammar_xbwt_trie.trie_bits_);
    trie_rank_1_ = std::move(grammar_xbwt_trie.trie_rank_1_);
    trie_rank_1_.set_vector(&trie_bits_);
    trie_select_1_ = std::move(grammar_xbwt_trie.trie_select_1_);
    trie_select_1_.set_vector(&trie_bits_);
    bwt_ = std::move(grammar_xbwt_trie.bwt_);
    symbol_bucket_offsets_ = std::move(grammar_xbwt_trie.symbol_bucket_offsets_);
    counts_ = std::move(grammar_xbwt_trie.counts_);
    lex_rank_range_ = std::move(grammar_xbwt_trie.lex_rank_range_);
  }
  return *this;
}

// uint64_t GrammarXbwtTrie::Serialize
// (
//   std::ostream& out,
//   std::shared_ptr<SpaceNode> parent,
//   std::string const name
// )
// {
//   uint64_t size_in_bytes {};
//   if (!parent)
//   {
//     sdsl::serialize(level_order_bits_, out);
//     sdsl::serialize(level_order_select_1_, out);
//     sdsl::serialize(labels_, out);
//     sdsl::serialize(counts_, out);
//   }
//   else
//   {
//     auto node {std::make_shared<SpaceNode>(name)};
//     node->AddLeaf("level_order_bits_", sdsl::serialize(level_order_bits_, out));
//     node->AddLeaf("level_order_select_1_", sdsl::serialize(level_order_select_1_, out));
//     node->AddLeaf("labels_", sdsl::serialize(labels_, out));
//     node->AddLeaf("counts_", sdsl::serialize(counts_, out));
//     parent->AddChild(node);
//     size_in_bytes = node->GetSizeInBytes();
//   }
//   return size_in_bytes;
// }

// void GrammarXbwtTrie::Load (std::istream& in)
// {
//   level_order_bits_.load(in);
//   level_order_select_1_.load(in);
//   level_order_select_1_.set_vector(&level_order_bits_);
//   labels_.load(in);
//   counts_.load(in);
//   return;
// }

// template <typename Iterator>
// uint64_t GrammarXbwtTrie::Count (Iterator it, Iterator end)
// {
// }

// std::ostream& operator<< (std::ostream& out, GrammarTrie const& grammar_trie)
// {
//   return out;
// }

void GrammarXbwtTrie::CalculateReverseGrammarRulesAndGrammarTrieNodePointers
(
  GrammarTrie const& grammar_trie,
  sdsl::int_vector<>& reverse_grammar_rules,
  std::deque<std::shared_ptr<GrammarTrie::Node>>& grammar_trie_node_pointers
)
{
  std::deque<std::pair<std::shared_ptr<GrammarTrie::Node>, bool>> nodes;
  std::deque<std::shared_ptr<GrammarTrie::Node>> ancestors;
  std::deque<uint64_t> temp_reverse_grammar_rules;
  auto root {grammar_trie.GetRoot()};
  for (auto it {std::rbegin(root->children)}; it != std::rend(root->children); ++it)
  {
    nodes.emplace_back(std::get<1>(*it), true);
  }
  while (!nodes.empty())
  {
    auto node {std::get<0>(nodes.back())};
    auto& is_forward {std::get<1>(nodes.back())};
    if (is_forward)
    {
      is_forward = false;
      ancestors.push_back(node);
      if
      (
        node->children.empty()
        ||
        std::get<0>(node->rank_range) < std::get<0>(std::get<1>(*std::begin(node->children))->rank_range)
      )
      {
        for (auto it {std::rbegin(ancestors)}; it != std::rend(ancestors); ++it)
        {
          temp_reverse_grammar_rules.push_back((*it)->label);
          grammar_trie_node_pointers.push_back(*it);
        }
        temp_reverse_grammar_rules.push_back(GrammarXbwtTrie::kSpecialSymbol);
        grammar_trie_node_pointers.push_back(node);
      }
      for (auto it {std::rbegin(node->children)}; it != std::rend(node->children); ++it)
      {
        nodes.emplace_back(std::get<1>(*it), true);
      }
    }
    else
    {
      nodes.pop_back();
      ancestors.pop_back();
    }
  }
  reverse_grammar_rules.width(grammar_trie.GetEffectiveAlphabetWidth() + grammar_trie.GetLengthWidth());
  reverse_grammar_rules.resize(std::size(temp_reverse_grammar_rules) + 1);
  std::copy
  (
    std::begin(temp_reverse_grammar_rules),
    std::end(temp_reverse_grammar_rules),
    std::begin(reverse_grammar_rules)
  );
  {
    uint64_t divisor {1ULL << grammar_trie.GetLengthWidth()};
    for (auto const& symbol : reverse_grammar_rules)
    {
      std::cout << "(" << (symbol / divisor) << "," << (symbol % divisor) << ")";
    }
    std::cout << "\n";
  }
  return;
}
}
