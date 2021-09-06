#pragma once

#include <deque>
#include <map>
#include <memory>

#include <sdsl/wavelet_trees.hpp>

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
    uint64_t rank;
    std::map<uint64_t, std::shared_ptr<Node>> children;

    Node (uint64_t const label = 0)
    : label {label},
      count {},
      rank {},
      children {}
    {
    }
  };

  GrammarTrie (uint8_t const effective_alphabet_width, uint64_t const length_width)
  : root_ {std::make_shared<Node>()},
    effective_alphabet_width_ {effective_alphabet_width},
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

  inline uint8_t GetEffectiveAlphabetWidth () const
  {
    return effective_alphabet_width_;
  }

  inline uint64_t GetLengthWidth () const
  {
    return length_width_;
  }

  inline uint64_t GetLabelWidth () const
  {
    return effective_alphabet_width_ + length_width_;
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
  node->rank = 1;
  return;
}

void GrammarTrie::CalculateCumulativeCountAndRank ()
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
      if (node->rank)
      {
        node->rank = ++rank;
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
  out << "depth:(byte_rank,length)(count)(rank)\n";
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
      << "(" << node->rank << ")\n";
    }
    for (auto it {std::rbegin(node->children)}; it != std::rend(node->children); ++it)
    {
      nodes.emplace_back(std::get<1>(*it), depth + 1);
    }
  }
  out << "|nodes|:" << (size - 1) << "\n";
  out << "label_width:" << grammar_trie.GetLabelWidth() << "\n";
  out << "length_width:" << grammar_trie.GetLengthWidth() << "\n";
  return out;
}

class GrammarXbwtTrie
{
public:

  static constexpr uint64_t kSpecialSymbol {1};

  GrammarXbwtTrie () = default;
  GrammarXbwtTrie (GrammarTrie const& grammar_trie);
  GrammarXbwtTrie& operator= (GrammarXbwtTrie&&);

  uint64_t Serialize
  (
    std::ostream& out,
    std::shared_ptr<SpaceNode> parent = nullptr,
    std::string const name = ""
  );
  void Load (std::istream& in);

  // template <typename Iterator>
  // uint64_t Count (Iterator it, Iterator end);

  friend std::ostream& operator<< (std::ostream& out, GrammarXbwtTrie const& grammar_xbwt_trie);

private:

  void CalculateReverseGrammarRulesAndGrammarTrieNodes
  (
    GrammarTrie const& grammar_trie,
    sdsl::int_vector<>& reverse_grammar_rules,
    std::deque<std::shared_ptr<GrammarTrie::Node>>& grammar_trie_nodes
  );
  void CalculateReducedCommomPrefixLengths
  (
    sdsl::int_vector<> const& reverse_grammar_rules,
    sdsl::int_vector<> const& suffix_array,
    sdsl::int_vector<>& reduced_commom_prefix_lengths
  );

  sdsl::bit_vector trie_bits_;
  sdsl::bit_vector::rank_1_type trie_rank_1_;
  sdsl::bit_vector::select_1_type trie_select_1_;
  sdsl::wm_int<> child_labels_;
  SparsePrefixSum label_bucket_offsets_;
  SparsePrefixSum cumulative_counts_;
  sdsl::int_vector<> colex_to_lex_rank_;

};

GrammarXbwtTrie::GrammarXbwtTrie (GrammarTrie const& grammar_trie)
{
  grammar_trie.GetLengthWidth();
  sdsl::int_vector<> reverse_grammar_rules;
  std::deque<std::shared_ptr<GrammarTrie::Node>> grammar_trie_nodes;
  CalculateReverseGrammarRulesAndGrammarTrieNodes
  (
    grammar_trie,
    reverse_grammar_rules,
    grammar_trie_nodes
  );
  // {
  //   uint64_t divisor {1ULL << grammar_trie.GetLengthWidth()};
  //   for (auto const& label : reverse_grammar_rules)
  //   {
  //     std::cout << "(" << (label / divisor) << "," << (label % divisor) << ")";
  //   }
  //   std::cout << "\n";
  //   std::cout << "(byte_rank,length)(count)(rank)\n";
  //   for (auto node : grammar_trie_nodes)
  //   {
  //     std::cout
  //     << "(" << (node->label / divisor) << "," << (node->label % divisor) << ")"
  //     << "(" << node->count << ")"
  //     << "(" << node->rank << ")\n";
  //   }
  // }
  sdsl::int_vector<> suffix_array;
  sdsl::qsufsort::construct_sa(suffix_array, reverse_grammar_rules);
  sdsl::int_vector<> reduced_commom_prefix_lengths;
  CalculateReducedCommomPrefixLengths(reverse_grammar_rules, suffix_array, reduced_commom_prefix_lengths);
  // {
  //   uint64_t divisor {1ULL << grammar_trie.GetLengthWidth()};
  //   std::cout << "i\tRCP\tBWT\tSA\treverse_grammar_rules[i:]\n";
  //   for (uint64_t i {}; i != std::size(suffix_array); ++i)
  //   {
  //     auto it {std::next(std::begin(reverse_grammar_rules), suffix_array[i])};
  //     std::cout << i << "\t";
  //     std::cout << reduced_commom_prefix_lengths[i] << "\t";
  //     if (suffix_array[i])
  //     {
  //       auto label {reverse_grammar_rules[suffix_array[i] - 1]};
  //       std::cout << "(" << (label / divisor) << "," << (label % divisor) << ")\t";
  //     }
  //     else
  //     {
  //       std::cout
  //       << "(" << (GrammarXbwtTrie::kSpecialSymbol / divisor)
  //       << "," << (GrammarXbwtTrie::kSpecialSymbol % divisor)
  //       << ")\t";
  //     }
  //     std::cout << suffix_array[i] << "\t";
  //     while (it != std::end(reverse_grammar_rules))
  //     {
  //       std::cout << "(" << (*it / divisor) << "," << (*it % divisor) << ")";
  //       ++it;
  //     }
  //     std::cout << "\n";
  //   }
  // }
  {
    std::deque<uint8_t> trie_bits;
    std::set<uint64_t> label_set;
    std::deque<uint64_t> child_labels;
    std::deque<uint64_t> label_bucket_offsets;
    std::deque<uint64_t> counts;
    std::deque<uint64_t> colex_to_lex_rank;
    {
      label_bucket_offsets.push_back(0);
      label_bucket_offsets.push_back(0);
      counts.push_back(0);
    }
    for (uint64_t i {1}, j {suffix_array[i]}; i != std::size(suffix_array); ++i, j = suffix_array[i])
    {
      uint64_t label {};
      label = ((j) ? reverse_grammar_rules[j - 1] : GrammarXbwtTrie::kSpecialSymbol);
      if (label == GrammarXbwtTrie::kSpecialSymbol)
      {
        colex_to_lex_rank.push_back(grammar_trie_nodes[j]->rank);
      }
      if (reverse_grammar_rules[j + reduced_commom_prefix_lengths[i]] == GrammarXbwtTrie::kSpecialSymbol)
      {
        if (label_set.count(label) == 0 || label == GrammarXbwtTrie::kSpecialSymbol)
        {
          trie_bits.push_back(0);
          ++label_bucket_offsets.back();
        }
      }
      else
      {
        if (!trie_bits.empty() && trie_bits.back() == 0)
        {
          trie_bits.back() = 1;
          trie_bits.push_back(0);
        }
        else
        {
          trie_bits.push_back(1);
        }
        for (auto const label : label_set)
        {
          child_labels.push_back(label);
        }
        label_set.clear();
        if (reverse_grammar_rules[j] != reverse_grammar_rules[suffix_array[i - 1]])
        {
          label_bucket_offsets.push_back(1);
        }
        else
        {
          ++label_bucket_offsets.back();
        }
        counts.push_back(grammar_trie_nodes[j]->count);
      }
      label_set.insert(label);
    }
    {
      trie_bits.back() = 1;
      trie_bits_.resize(std::size(trie_bits));
      std::copy(std::begin(trie_bits), std::end(trie_bits), std::begin(trie_bits_));
      trie_rank_1_ = decltype(trie_rank_1_)(&trie_bits_);
      trie_select_1_ = decltype(trie_select_1_)(&trie_bits_);
      // {
      //   std::cout << "trie_bits: ";
      //   Print(trie_bits, std::cout);
      // }
    }
    {
      child_labels.push_back(reverse_grammar_rules[*std::prev(std::end(suffix_array)) - 1]);
      sdsl::int_vector<> temp_child_labels(std::size(child_labels), 0, grammar_trie.GetLabelWidth());
      std::copy(std::begin(child_labels), std::end(child_labels), std::begin(temp_child_labels));
      sdsl::construct_im(child_labels_, temp_child_labels);
      // {
      //   uint64_t divisor {1ULL << grammar_trie.GetLengthWidth()};
      //   std::cout << "child_labels: ";
      //   for (auto const &label : child_labels)
      //   {
      //     std::cout << "(" << (label / divisor) << "," << (label % divisor) << ")";
      //   }
      //   std::cout << "\n";
      // }
    }
    {
      std::partial_sum(std::begin(label_bucket_offsets), std::end(label_bucket_offsets), std::begin(label_bucket_offsets));
      label_bucket_offsets_ = decltype(label_bucket_offsets_)(label_bucket_offsets);
      // {
      //   std::cout << "label_bucket_offsets: ";
      //   Print(label_bucket_offsets, std::cout);
      // }
    }
    {
      std::partial_sum(std::begin(counts), std::end(counts), std::begin(counts));
      cumulative_counts_ = decltype(cumulative_counts_)(counts);
      // {
      //   std::cout << "counts: ";
      //   Print(counts, std::cout);
      // }
    }
    {
      colex_to_lex_rank_.width(sdsl::bits::hi(std::size(colex_to_lex_rank)) + 1);
      colex_to_lex_rank_.resize(std::size(colex_to_lex_rank));
      std::copy(std::begin(colex_to_lex_rank), std::end(colex_to_lex_rank), std::begin(colex_to_lex_rank_));
      // {
      //   std::cout << "colex_to_lex_rank: ";
      //   Print(colex_to_lex_rank, std::cout);
      // }
    }
  }
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
    child_labels_ = std::move(grammar_xbwt_trie.child_labels_);
    label_bucket_offsets_ = std::move(grammar_xbwt_trie.label_bucket_offsets_);
    cumulative_counts_ = std::move(grammar_xbwt_trie.cumulative_counts_);
    colex_to_lex_rank_ = std::move(grammar_xbwt_trie.colex_to_lex_rank_);
  }
  return *this;
}

uint64_t GrammarXbwtTrie::Serialize
(
  std::ostream& out,
  std::shared_ptr<SpaceNode> parent,
  std::string const name
)
{
  uint64_t size_in_bytes {};
  if (!parent)
  {
    sdsl::serialize(trie_bits_, out);
    sdsl::serialize(trie_rank_1_, out);
    sdsl::serialize(trie_select_1_, out);
    sdsl::serialize(child_labels_, out);
    label_bucket_offsets_.Serialize(out);
    cumulative_counts_.Serialize(out);
    sdsl::serialize(colex_to_lex_rank_, out);
  }
  else
  {
    auto node {std::make_shared<SpaceNode>(name)};
    node->AddLeaf("trie_bits_", sdsl::serialize(trie_bits_, out));
    node->AddLeaf("trie_rank_1_", sdsl::serialize(trie_rank_1_, out));
    node->AddLeaf("trie_select_1_", sdsl::serialize(trie_select_1_, out));
    node->AddLeaf("child_labels_", sdsl::serialize(child_labels_, out));
    node->AccumalateSizeInBytes(label_bucket_offsets_.Serialize(out, node, "label_bucket_offsets_"));
    node->AccumalateSizeInBytes(cumulative_counts_.Serialize(out, node, "cumulative_counts_"));
    node->AddLeaf("colex_to_lex_rank_", sdsl::serialize(colex_to_lex_rank_, out));
    parent->AddChild(node);
    size_in_bytes = node->GetSizeInBytes();
  }
  return size_in_bytes;
}

void GrammarXbwtTrie::Load (std::istream& in)
{
  trie_bits_.load(in);
  trie_rank_1_.load(in);
  trie_rank_1_.set_vector(&trie_bits_);
  trie_select_1_.load(in);
  trie_select_1_.set_vector(&trie_bits_);
  child_labels_.load(in);
  label_bucket_offsets_.Load(in);
  cumulative_counts_.Load(in);
  colex_to_lex_rank_.load(in);;
  return;
}

// template <typename Iterator>
// uint64_t GrammarXbwtTrie::Count (Iterator it, Iterator end)
// {
// }

std::ostream& operator<< (std::ostream& out, GrammarXbwtTrie const& grammar_xbwt_trie)
{
  // out << "xbwt_trie:\n";
  // out << "space:\n";
  // out << "trie_bits_: " << ProperSizeRepresentation(sdsl::size_in_bytes(grammar_xbwt_trie.trie_bits_)) << "B\n";
  // out << "trie_rank_1_: " << ProperSizeRepresentation(sdsl::size_in_bytes(grammar_xbwt_trie.trie_rank_1_)) << "B\n";
  // out << "trie_select_1_: " << ProperSizeRepresentation(sdsl::size_in_bytes(grammar_xbwt_trie.trie_select_1_)) << "B\n";
  // out << "child_labels_: " << ProperSizeRepresentation(sdsl::size_in_bytes(grammar_xbwt_trie.child_labels_)) << "B\n";
  // out << "label_bucket_offsets_:\n";
  // out << grammar_xbwt_trie.label_bucket_offsets_;
  // out << "cumulative_counts_:\n";
  // out << grammar_xbwt_trie.cumulative_counts_;
  // out << "colex_to_lex_rank_: " << ProperSizeRepresentation(sdsl::size_in_bytes(grammar_xbwt_trie.colex_to_lex_rank_)) << "B\n";
  return out;
}

void GrammarXbwtTrie::CalculateReverseGrammarRulesAndGrammarTrieNodes
(
  GrammarTrie const& grammar_trie,
  sdsl::int_vector<>& reverse_grammar_rules,
  std::deque<std::shared_ptr<GrammarTrie::Node>>& grammar_trie_nodes
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
      if (node->rank)
      {
        for (auto it {std::rbegin(ancestors)}; it != std::rend(ancestors); ++it)
        {
          temp_reverse_grammar_rules.push_back((*it)->label);
          grammar_trie_nodes.push_back(*it);
        }
        temp_reverse_grammar_rules.push_back(GrammarXbwtTrie::kSpecialSymbol);
        grammar_trie_nodes.push_back(node);
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
  reverse_grammar_rules.width(grammar_trie.GetLabelWidth());
  reverse_grammar_rules.resize(std::size(temp_reverse_grammar_rules) + 1);
  std::copy(std::begin(temp_reverse_grammar_rules), std::end(temp_reverse_grammar_rules), std::begin(reverse_grammar_rules));
  return;
}

void GrammarXbwtTrie::CalculateReducedCommomPrefixLengths
(
  sdsl::int_vector<> const& reverse_grammar_rules,
  sdsl::int_vector<> const& suffix_array,
  sdsl::int_vector<>& reduced_commom_prefix_lengths
)
{
  sdsl::int_vector<> inverse_suffix_array(std::size(suffix_array));
  for (uint64_t i {}; i != std::size(inverse_suffix_array); ++i)
  {
    inverse_suffix_array[suffix_array[i]] = i;
  }
  {
    reduced_commom_prefix_lengths.resize(std::size(suffix_array));
    *std::prev(std::end(reduced_commom_prefix_lengths)) = 0;
    uint64_t length {};
    for (uint64_t i {}; i != (std::size(reverse_grammar_rules) - 1); ++i)
    {
      if (length) { --length; }
      auto suffix_it {std::next(std::begin(reverse_grammar_rules), i + length)};
      uint64_t j {suffix_array[inverse_suffix_array[i] - 1]};
      auto prev_suffix_it {std::next(std::begin(reverse_grammar_rules), j + length)};
      while
      (
        (suffix_it != std::end(reverse_grammar_rules))
        && (prev_suffix_it != std::end(reverse_grammar_rules))
        && (*suffix_it == *prev_suffix_it)
        && (*suffix_it != GrammarXbwtTrie::kSpecialSymbol)
      )
      {
        ++length;
        ++suffix_it;
        ++prev_suffix_it;
      }
      reduced_commom_prefix_lengths[inverse_suffix_array[i]] = length;
    }
  }
  return;
}

}
