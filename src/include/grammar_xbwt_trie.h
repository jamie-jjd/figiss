#pragma once

#include <deque>
#include <map>
#include <memory>

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

  inline uint8_t GetLengthWidth () const
  {
    return length_width_;
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
  uint64_t const divisor {1ULL << grammar_trie.length_width_};
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
  out << "|nodes|:" << size << "\n";
  out << "length_width:" << static_cast<uint64_t>(grammar_trie.length_width_) << "\n";
  return out;
}

class GrammarXbwtTrie
{
public:

  static constexpr uint64_t kSpecialSymbol {1};

  GrammarXbwtTrie () = default;
  GrammarXbwtTrie (GrammarTrie const& grammar_trie, IntegerAlphabet&& run_length_alphabet);
  GrammarXbwtTrie (GrammarXbwtTrie const&) = delete;
  GrammarXbwtTrie (GrammarXbwtTrie&&);
  GrammarXbwtTrie& operator= (GrammarXbwtTrie const &) = delete;
  GrammarXbwtTrie& operator= (GrammarXbwtTrie&&);
  ~GrammarXbwtTrie () = default;

  void Swap (GrammarXbwtTrie&);

  uint64_t Serialize
  (
    std::ostream& out,
    std::shared_ptr<SpaceNode> parent = nullptr,
    std::string const name = ""
  );
  void Load (std::istream& in);

  template <typename Iterator>
  uint64_t Rank (Iterator it, Iterator end, bool const is_lex = false) const; // true: colex_rank
  template <typename Iterator>
  std::pair<uint64_t, uint64_t> GetLexRankRange (Iterator it, Iterator end) const;
  template <typename Iterator>
  std::pair<uint64_t, uint64_t> GetColexRankRange (Iterator it, Iterator end) const;
  template <typename Iterator>
  uint64_t Count (Iterator it, Iterator end) const;

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

  inline uint64_t GetLexAlphabetSize () const
  {
    return (std::size(colex_to_lex_rank_) + 1);
  }

  inline uint64_t ColexToLexRank (uint64_t const colex_rank) const
  {
    if (colex_rank && colex_rank <= std::size(colex_to_lex_rank_))
    {
      return colex_to_lex_rank_[colex_rank - 1];
    }
    return 0;
  }

  friend std::ostream& operator<< (std::ostream& out, GrammarXbwtTrie const& xbwt);

private:

  void CalculateReverseGrammarRulesAndGrammarTrieNodesAndCumulativeByteRankCounts
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
  uint64_t CalculateChildBeginOffset
  (
    std::pair<uint64_t, uint64_t> const& children_offset_range,
    uint64_t label_rank
  ) const;
  uint64_t CalculateChildRbeginOffset
  (
    std::pair<uint64_t, uint64_t> const& children_offset_range,
    uint64_t label_rank
  ) const;

  template <typename Range> // [,]
  inline bool IsEmptyRange (Range const& range) const
  {
    return (std::get<0>(range) > std::get<1>(range));
  }

  inline uint64_t InvalidOffset () const
  {
    return std::size(children_label_ranks_);
  }

  sdsl::bit_vector children_bits_;
  sdsl::bit_vector::rank_1_type children_rank_1_;
  sdsl::bit_vector::select_1_type children_select_1_;
  uint8_t length_width_;
  IntegerAlphabet run_length_alphabet_;
  sdsl::wm_int<> children_label_ranks_;
  SparsePrefixSum label_rank_bucket_offsets_;
  SparsePrefixSum cumulative_counts_;
  SparsePrefixSum cumulative_byte_rank_counts_;
  sdsl::int_vector<> colex_to_lex_rank_;

};

GrammarXbwtTrie::GrammarXbwtTrie (GrammarTrie const& grammar_trie, IntegerAlphabet&& run_length_alphabet)
{
  {
    length_width_ = grammar_trie.GetLengthWidth();
    // std::cout << "length_width_:" << static_cast<uint64_t>(length_width_) << "\n";
    run_length_alphabet_ = std::move(run_length_alphabet);
    // std::cout << "run_length_alphabet_:\n" << run_length_alphabet_;
  }
  sdsl::int_vector<> reverse_grammar_rules;
  std::deque<std::shared_ptr<GrammarTrie::Node>> grammar_trie_nodes;
  CalculateReverseGrammarRulesAndGrammarTrieNodesAndCumulativeByteRankCounts
  (
    grammar_trie,
    reverse_grammar_rules,
    grammar_trie_nodes
  );
  // {
  //   uint64_t const divisor {1ULL << length_width_};
  //   for (auto const& label : reverse_grammar_rules)
  //   {
  //     std::cout << "(" << (label / divisor) << "," << (label % divisor) << ")";
  //   }
  //   std::cout << "\n";
  //   std::cout << "(byte_rank,length)(count)(lex_rank)\n";
  //   for (auto node : grammar_trie_nodes)
  //   {
  //     std::cout
  //     << "(" << (node->label / divisor) << "," << (node->label % divisor) << ")"
  //     << "(" << node->count << ")"
  //     << "(" << node->lex_rank << ")\n";
  //   }
  //   std::cout << "cumulative_byte_rank_counts_:" << cumulative_byte_rank_counts_;
  // }
  sdsl::int_vector<> suffix_array;
  sdsl::qsufsort::construct_sa(suffix_array, reverse_grammar_rules);
  sdsl::int_vector<> reduced_commom_prefix_lengths;
  CalculateReducedCommomPrefixLengths(reverse_grammar_rules, suffix_array, reduced_commom_prefix_lengths);
  // {
  //   uint64_t const divisor {1ULL << length_width_};
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
    std::deque<uint8_t> children_bits;
    std::set<uint64_t> label_set;
    std::deque<uint64_t> children_labels;
    std::deque<uint64_t> label_rank_bucket_offsets(1, 0);
    std::deque<uint64_t> counts;
    std::deque<uint64_t> colex_to_lex_rank;
    for (uint64_t i {1}, j {suffix_array[i]}; i != std::size(suffix_array); ++i, j = suffix_array[i])
    {
      uint64_t label {};
      label = ((j) ? reverse_grammar_rules[j - 1] : GrammarXbwtTrie::kSpecialSymbol);
      if (label == GrammarXbwtTrie::kSpecialSymbol)
      {
        colex_to_lex_rank.push_back(grammar_trie_nodes[j]->lex_rank);
      }
      if (reverse_grammar_rules[j + reduced_commom_prefix_lengths[i]] == GrammarXbwtTrie::kSpecialSymbol)
      {
        if (label_set.count(label) == 0 || label == GrammarXbwtTrie::kSpecialSymbol)
        {
          children_bits.push_back(0);
          ++label_rank_bucket_offsets.back();
        }
      }
      else
      {
        if (!children_bits.empty() && children_bits.back() == 0)
        {
          children_bits.back() = 1;
          children_bits.push_back(0);
        }
        else
        {
          children_bits.push_back(1);
        }
        for (auto const label : label_set)
        {
          children_labels.push_back(label);
        }
        label_set.clear();
        if (reverse_grammar_rules[j] != reverse_grammar_rules[suffix_array[i - 1]])
        {
          label_rank_bucket_offsets.push_back(1);
        }
        else
        {
          ++label_rank_bucket_offsets.back();
        }
        counts.push_back(grammar_trie_nodes[j]->count);
      }
      label_set.insert(label);
    }
    {
      children_bits.back() = 1;
      children_bits_.resize(std::size(children_bits));
      std::copy(std::begin(children_bits), std::end(children_bits), std::begin(children_bits_));
      children_rank_1_ = decltype(children_rank_1_)(&children_bits_);
      children_select_1_ = decltype(children_select_1_)(&children_bits_);
      // {
      //   std::cout << "children_bits_:";
      //   Print(children_bits_, std::cout);
      // }
    }
    {
      for (auto const& label : label_set)
      {
        children_labels.push_back(label);
      }
      // {
      //   uint64_t const divisor {1ULL << length_width_};
      //   std::cout << "children_labels:";
      //   for (auto const &label : children_labels)
      //   {
      //     std::cout << "(" << (label / divisor) << "," << (label % divisor) << ")";
      //   }
      //   std::cout << "\n";
      // }
      sdsl::int_vector<> children_label_ranks(std::size(children_labels), 0, run_length_alphabet_.GetEffectiveAlphabetWidth());
      std::transform
      (
        std::begin(children_labels),
        std::end(children_labels),
        std::begin(children_label_ranks),
        [&] (auto const& label)
        {
          return run_length_alphabet_[label];
        }
      );
      sdsl::construct_im(children_label_ranks_, children_label_ranks);
      // {
      //   std::cout << "children_label_ranks_:";
      //   Print(children_label_ranks_, std::cout);
      // }
    }
    {
      std::partial_sum(std::begin(label_rank_bucket_offsets), std::end(label_rank_bucket_offsets), std::begin(label_rank_bucket_offsets));
      label_rank_bucket_offsets_ = decltype(label_rank_bucket_offsets_)(label_rank_bucket_offsets);
      // std::cout << "label_rank_bucket_offsets_:" << label_rank_bucket_offsets_;
    }
    {
      std::partial_sum(std::begin(counts), std::end(counts), std::begin(counts));
      cumulative_counts_ = decltype(cumulative_counts_)(counts);
      // std::cout << "cumulative_counts:" << cumulative_counts_;
    }
    {
      colex_to_lex_rank_.width(sdsl::bits::hi(std::size(colex_to_lex_rank)) + 1);
      colex_to_lex_rank_.resize(std::size(colex_to_lex_rank));
      std::copy(std::begin(colex_to_lex_rank), std::end(colex_to_lex_rank), std::begin(colex_to_lex_rank_));
      // {
      //   std::cout << "colex_to_lex_rank:";
      //   Print(colex_to_lex_rank, std::cout);
      // }
    }
  }
}

GrammarXbwtTrie::GrammarXbwtTrie (GrammarXbwtTrie&& grammar_xbwt_trie)
{
  if (this != &grammar_xbwt_trie)
  {
    this->Swap(grammar_xbwt_trie);
  }
}

GrammarXbwtTrie& GrammarXbwtTrie::operator= (GrammarXbwtTrie&& grammar_xbwt_trie)
{
  if (this != &grammar_xbwt_trie)
  {
    GrammarXbwtTrie temp {std::move(grammar_xbwt_trie)};
    this->Swap(temp);
  }
  return *this;
}

void GrammarXbwtTrie::Swap (GrammarXbwtTrie& grammar_xbwt_trie)
{
  if (this != &grammar_xbwt_trie)
  {
    children_bits_.swap(grammar_xbwt_trie.children_bits_);
    children_rank_1_.swap(grammar_xbwt_trie.children_rank_1_);
    children_rank_1_.set_vector(&children_bits_);
    grammar_xbwt_trie.children_rank_1_.set_vector(&grammar_xbwt_trie.children_bits_);
    children_select_1_.swap(grammar_xbwt_trie.children_select_1_);
    children_select_1_.set_vector(&children_bits_);
    grammar_xbwt_trie.children_select_1_.set_vector(&grammar_xbwt_trie.children_bits_);
    std::swap(length_width_, grammar_xbwt_trie.length_width_);
    run_length_alphabet_.Swap(grammar_xbwt_trie.run_length_alphabet_);
    children_label_ranks_.swap(grammar_xbwt_trie.children_label_ranks_);
    label_rank_bucket_offsets_.Swap(grammar_xbwt_trie.label_rank_bucket_offsets_);
    cumulative_counts_.Swap(grammar_xbwt_trie.cumulative_counts_);
    colex_to_lex_rank_.swap(grammar_xbwt_trie.colex_to_lex_rank_);
    cumulative_byte_rank_counts_.Swap(grammar_xbwt_trie.cumulative_byte_rank_counts_);
  }
  return;
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
    sdsl::serialize(children_bits_, out);
    sdsl::serialize(children_rank_1_, out);
    sdsl::serialize(children_select_1_, out);
    sdsl::write_member(length_width_, out);
    run_length_alphabet_.Serialize(out);
    sdsl::serialize(children_label_ranks_, out);
    label_rank_bucket_offsets_.Serialize(out);
    cumulative_counts_.Serialize(out);
    sdsl::serialize(colex_to_lex_rank_, out);
    cumulative_byte_rank_counts_.Serialize(out);
  }
  else
  {
    auto node {std::make_shared<SpaceNode>(name)};
    node->AddLeaf("children_bits_", sdsl::serialize(children_bits_, out));
    node->AddLeaf("children_rank_1_", sdsl::serialize(children_rank_1_, out));
    node->AddLeaf("children_select_1_", sdsl::serialize(children_select_1_, out));
    node->AddLeaf("length_width_", sdsl::write_member(length_width_, out));
    node->AccumalateSizeInBytes(run_length_alphabet_.Serialize(out, node, "run_length_alphabet_"));
    node->AddLeaf("children_label_ranks_", sdsl::serialize(children_label_ranks_, out));
    node->AccumalateSizeInBytes(label_rank_bucket_offsets_.Serialize(out, node, "label_rank_bucket_offsets_"));
    node->AccumalateSizeInBytes(cumulative_counts_.Serialize(out, node, "cumulative_counts_"));
    node->AddLeaf("colex_to_lex_rank_", sdsl::serialize(colex_to_lex_rank_, out));
    node->AccumalateSizeInBytes(cumulative_byte_rank_counts_.Serialize(out, node, "cumulative_byte_rank_counts_"));
    parent->AddChild(node);
    size_in_bytes = node->GetSizeInBytes();
  }
  return size_in_bytes;
}

void GrammarXbwtTrie::Load (std::istream& in)
{
  children_bits_.load(in);
  children_rank_1_.load(in);
  children_rank_1_.set_vector(&children_bits_);
  children_select_1_.load(in);
  children_select_1_.set_vector(&children_bits_);
  sdsl::read_member(length_width_, in);
  run_length_alphabet_.Load(in);
  children_label_ranks_.load(in);
  label_rank_bucket_offsets_.Load(in);
  cumulative_counts_.Load(in);
  colex_to_lex_rank_.load(in);;
  cumulative_byte_rank_counts_.Load(in);
  return;
}

template <typename Iterator>
std::pair<uint64_t, uint64_t> GrammarXbwtTrie::GetLexRankRange (Iterator it, Iterator end) const
{
  uint64_t const divisor {1ULL << length_width_};
  auto lex_rank_range {std::make_pair<uint64_t, uint64_t>(1, 0)};
  uint64_t offset {};
  uint64_t label_rank {};
  while (it != std::prev(end))
  {
    label_rank = run_length_alphabet_[*it++];
    if (label_rank)
    {
      offset = CalculateChildBeginOffset({offset, children_select_1_(children_rank_1_(offset) + 1)}, label_rank);
      if (offset == InvalidOffset()) { return {1, 0}; }
    }
    else { return {1, 0}; }
  }
  {
    auto temp_offset {offset};
    label_rank = run_length_alphabet_.Successor(*it - 1);
    if (label_rank && (*it / divisor) == (run_length_alphabet_.Select(label_rank) / divisor))
    {
      temp_offset = CalculateChildBeginOffset({temp_offset, children_select_1_(children_rank_1_(temp_offset) + 1)}, label_rank);
      if (temp_offset == InvalidOffset()) { return {1, 0}; }
      while ((label_rank = children_label_ranks_[temp_offset]))
      {
        temp_offset = CalculateChildBeginOffset({temp_offset, children_select_1_(children_rank_1_(temp_offset) + 1)}, label_rank);
        if (temp_offset == InvalidOffset()) { return {1, 0}; }
      }
      std::get<0>(lex_rank_range) = colex_to_lex_rank_[children_label_ranks_.rank(temp_offset, 0)];
      {
        temp_offset = offset;
        label_rank = run_length_alphabet_.Predecessor((*it / divisor + 1) * divisor);
        temp_offset = CalculateChildRbeginOffset({temp_offset, children_select_1_(children_rank_1_(temp_offset) + 1)}, label_rank);
        while ((label_rank = children_label_ranks_[temp_offset]))
        {
          temp_offset = CalculateChildRbeginOffset({temp_offset, children_select_1_(children_rank_1_(temp_offset) + 1)}, label_rank);
          if (temp_offset == InvalidOffset()) { return {1, 0}; }
        }
        std::get<1>(lex_rank_range) = colex_to_lex_rank_[children_label_ranks_.rank(temp_offset, 0)];
      }
    }
    else { return {1, 0}; }
  }
  return lex_rank_range;
}

template <typename Iterator>
uint64_t GrammarXbwtTrie::Rank (Iterator it, Iterator end, bool const is_lex) const
{
  uint64_t offset {};
  while (it != end)
  {
    auto label_rank {run_length_alphabet_[*it++]};
    if (label_rank)
    {
      offset = CalculateChildBeginOffset({offset, children_select_1_(children_rank_1_(offset) + 1)}, label_rank);
      if (offset == InvalidOffset()) { return 0; }
    }
    else { return 0; }
  }
  if (children_label_ranks_[offset] == 0)
  {
    if (is_lex)
    {
      return colex_to_lex_rank_[children_label_ranks_.rank(offset, 0)];
    }
    return children_label_ranks_.rank(offset + 1, 0);
  }
  return 0;
}

template <typename Iterator>
std::pair<uint64_t, uint64_t> GrammarXbwtTrie::GetColexRankRange (Iterator it, Iterator end) const
{
  uint64_t const divisor {1ULL << length_width_};
  auto offset_range {std::make_pair<uint64_t, uint64_t>(1, 0)};
  auto label_rank_range {std::make_pair<uint64_t, uint64_t>(1, 0)};
  std::get<0>(label_rank_range) = run_length_alphabet_.Successor(*it - 1);
  if
  (
    std::get<0>(label_rank_range)
    && (*it / divisor) == (run_length_alphabet_.Select(std::get<0>(label_rank_range)) / divisor)
  )
  {
    std::get<1>(label_rank_range) = run_length_alphabet_.Predecessor((*it / divisor + 1) * divisor);
    offset_range =
    {
      label_rank_bucket_offsets_[std::get<0>(label_rank_range)],
      label_rank_bucket_offsets_[std::get<1>(label_rank_range) + 1] - 1
    };
    if (IsEmptyRange(offset_range)) { return {1, 0}; }
  }
  else { return {1, 0}; }
  while (++it != end)
  {
    auto label_rank {run_length_alphabet_[*it]};
    if (label_rank)
    {
      auto begin_offset {label_rank_bucket_offsets_[label_rank]};
      offset_range =
      {
        begin_offset + children_label_ranks_.rank(std::get<0>(offset_range), label_rank),
        begin_offset + children_label_ranks_.rank(std::get<1>(offset_range) + 1, label_rank) - 1
      };
      if (IsEmptyRange(offset_range)) { return {1, 0}; }
    }
    else { return {1, 0}; }
  }
  return
  {
    children_label_ranks_.rank(std::get<0>(offset_range), 0) + 1,
    children_label_ranks_.rank(std::get<1>(offset_range) + 1, 0)
  };
}

template <typename Iterator>
uint64_t GrammarXbwtTrie::Count (Iterator it, Iterator end) const
{
  uint64_t const divisor {1ULL << length_width_};
  auto offset_range {std::make_pair<uint64_t, uint64_t>(1, 0)};
  auto label_rank_range {std::make_pair<uint64_t, uint64_t>(1, 0)};
  std::get<0>(label_rank_range) = run_length_alphabet_.Successor(*it - 1);
  if
  (
    std::get<0>(label_rank_range)
    && (*it / divisor) == (run_length_alphabet_.Select(std::get<0>(label_rank_range)) / divisor)
  )
  {
    std::get<1>(label_rank_range) = run_length_alphabet_.Predecessor((*it / divisor + 1) * divisor);
    offset_range =
    {
      label_rank_bucket_offsets_[std::get<0>(label_rank_range)],
      label_rank_bucket_offsets_[std::get<1>(label_rank_range) + 1] - 1
    };
    if (IsEmptyRange(offset_range)) { return 0; }
    if (it == std::prev(end))
    {
      return
      (
        cumulative_byte_rank_counts_[std::get<1>(label_rank_range)]
        - cumulative_byte_rank_counts_[std::get<0>(label_rank_range) - 1]
      )
      - ((*it % divisor) - 1) *
      (
        cumulative_counts_[children_rank_1_(std::get<1>(offset_range))]
        - cumulative_counts_[children_rank_1_(std::get<0>(offset_range)) - 1]
      );
    }
  }
  else { return 0; }
  while (++it != std::prev(end))
  {
    auto label_rank {run_length_alphabet_[*it]};
    if (label_rank)
    {
      auto begin_offset {label_rank_bucket_offsets_[label_rank]};
      offset_range =
      {
        begin_offset + children_label_ranks_.rank(std::get<0>(offset_range), label_rank),
        begin_offset + children_label_ranks_.rank(std::get<1>(offset_range) + 1, label_rank) - 1
      };
      if (IsEmptyRange(offset_range)) { return 0; }
    } else { return 0; }
  }
  uint64_t count {};
  std::get<0>(label_rank_range) = run_length_alphabet_.Successor(*it - 1);
  if
  (
    std::get<0>(label_rank_range)
    && (*it / divisor) == (run_length_alphabet_.Select(std::get<0>(label_rank_range)) / divisor)
  )
  {
    auto label_rank {std::get<0>(label_rank_range)};
    std::get<1>(label_rank_range) = run_length_alphabet_.Predecessor((*it / divisor + 1) * divisor);
    while (label_rank <= std::get<1>(label_rank_range))
    {
      auto begin_offset {label_rank_bucket_offsets_[label_rank]};
      auto temp_offset_range
      {
        std::make_pair
        (
          begin_offset + children_label_ranks_.rank(std::get<0>(offset_range), label_rank),
          begin_offset + children_label_ranks_.rank(std::get<1>(offset_range) + 1, label_rank) - 1
        )
      };
      if (!IsEmptyRange(temp_offset_range))
      {
        count +=
        (
          cumulative_counts_[children_rank_1_(std::get<1>(temp_offset_range))]
          - cumulative_counts_[children_rank_1_(std::get<0>(temp_offset_range)) - 1]
        );
      }
      ++label_rank;
    }
  }
  return count;
}

std::ostream& operator<< (std::ostream& out, GrammarXbwtTrie const& xbwt)
{
  uint64_t size {};
  out << "depth:(byte_rank,length)(count)(lex_rank)\n";
  std::deque<std::pair<uint64_t, uint64_t>> offsets;
  offsets.emplace_back(0, 0);
  while (!offsets.empty())
  {
    ++size;
    auto offset {std::get<0>(offsets.back())};
    auto depth {std::get<1>(offsets.back())};
    offsets.pop_back();
    if (depth)
    {
      auto label_rank {xbwt.label_rank_bucket_offsets_.Predecessor(offset + 1)};
      auto label {xbwt.run_length_alphabet_.Select(label_rank)};
      uint64_t const divisor {1ULL << xbwt.length_width_};
      auto children_rank {xbwt.children_rank_1_(offset)};
      uint64_t count {xbwt.cumulative_counts_[children_rank] - xbwt.cumulative_counts_[children_rank - 1]};
      uint64_t lex_rank {};
      if (xbwt.children_label_ranks_[offset] == 0)
      {
        lex_rank = xbwt.colex_to_lex_rank_[xbwt.children_label_ranks_.rank(offset, 0)];
      }
      std::cout
      << depth << ":"
      << "(" << (label / divisor) << "," << (label % divisor) << ")"
      << "(" << count << ")"
      << "(" << lex_rank << ")\n";
    }
    {
      auto children_offset {xbwt.children_select_1_(xbwt.children_rank_1_(offset) + 1)};
      auto children_rend_offset {offset - 1};
      while (children_offset != children_rend_offset)
      {
        auto pair {xbwt.children_label_ranks_.inverse_select(children_offset)};
        auto occurrence {std::get<0>(pair)};
        auto label_rank {std::get<1>(pair)};
        if (label_rank)
        {
          auto bucket_offset {xbwt.label_rank_bucket_offsets_[label_rank]};
          offsets.emplace_back(xbwt.children_select_1_(xbwt.children_rank_1_(bucket_offset) + occurrence) + 1, depth + 1);
        }
        --children_offset;
      }
    }
  }
  out << "|nodes|:" << size << "\n";
  out << "length_width:" << static_cast<uint64_t>(xbwt.GetLengthWidth()) << "\n";
  out << "alphabet_width:" << static_cast<uint64_t>(xbwt.GetAlphabetWidth()) << "\n";
  out << "effetive_alphabet_width:" << static_cast<uint64_t>(xbwt.GetEffectiveAlphabetWidth()) << "\n";
  return out;
}

void GrammarXbwtTrie::CalculateReverseGrammarRulesAndGrammarTrieNodesAndCumulativeByteRankCounts
(
  GrammarTrie const& grammar_trie,
  sdsl::int_vector<>& reverse_grammar_rules,
  std::deque<std::shared_ptr<GrammarTrie::Node>>& grammar_trie_nodes
)
{
  uint64_t const divisor {1ULL << grammar_trie.GetLengthWidth()};
  std::deque<std::pair<std::shared_ptr<GrammarTrie::Node>, bool>> nodes;
  std::deque<std::shared_ptr<GrammarTrie::Node>> ancestors;
  std::deque<uint64_t> temp_reverse_grammar_rules;
  std::map<uint64_t, uint64_t> byte_rank_counts;
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
      if (node->lex_rank)
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
      {
        uint64_t byte_rank_count {(node->label % divisor) * node->count};
        {
          auto it {byte_rank_counts.find(node->label)};
          if (it != byte_rank_counts.end())
          {
            std::get<1>(*it) += byte_rank_count;
          }
          else
          {
            byte_rank_counts[node->label] = byte_rank_count;
          }
        }
      }
    }
    else
    {
      nodes.pop_back();
      ancestors.pop_back();
    }
  }
  {
    reverse_grammar_rules.width(run_length_alphabet_.GetAlphabetWidth());
    reverse_grammar_rules.resize(std::size(temp_reverse_grammar_rules) + 1);
    std::copy(std::begin(temp_reverse_grammar_rules), std::end(temp_reverse_grammar_rules), std::begin(reverse_grammar_rules));
  }
  {
    uint64_t prefix_sum {};
    std::deque<uint64_t> cumulative_byte_rank_counts;
    for (auto const& pair : byte_rank_counts)
    {
      prefix_sum += std::get<1>(pair);
      cumulative_byte_rank_counts.push_back(prefix_sum);
    }
    cumulative_byte_rank_counts_ = decltype(cumulative_byte_rank_counts_)(cumulative_byte_rank_counts);
  }
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

uint64_t GrammarXbwtTrie::CalculateChildBeginOffset
(
  std::pair<uint64_t, uint64_t> const& children_offset_range,
  uint64_t const label_rank
) const
{
  auto offset {InvalidOffset()};
  auto occurrence_range
  {
    std::make_pair
    (
      children_label_ranks_.rank(std::get<0>(children_offset_range), label_rank),
      children_label_ranks_.rank(std::get<1>(children_offset_range) + 1, label_rank) - 1
    )
  };
  if (!IsEmptyRange(occurrence_range))
  {
    auto bucket_offset {label_rank_bucket_offsets_[label_rank]};
    offset = children_select_1_(children_rank_1_(bucket_offset) + std::get<0>(occurrence_range)) + 1;
  }
  return offset;
}

uint64_t GrammarXbwtTrie::CalculateChildRbeginOffset
(
  std::pair<uint64_t, uint64_t> const& children_offset_range,
  uint64_t const label_rank
) const
{
  auto offset {InvalidOffset()};
  auto occurrence_range
  {
    std::make_pair
    (
      children_label_ranks_.rank(std::get<0>(children_offset_range), label_rank),
      children_label_ranks_.rank(std::get<1>(children_offset_range) + 1, label_rank) - 1
    )
  };
  if (!IsEmptyRange(occurrence_range))
  {
    auto bucket_offset {label_rank_bucket_offsets_[label_rank]};
    offset = children_select_1_(children_rank_1_(bucket_offset) + std::get<0>(occurrence_range) + 1);
  }
  return offset;
}
}
