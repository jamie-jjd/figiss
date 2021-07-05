#pragma once

#include <deque>
#include <map>
#include <memory>
#include <set>

#include <sdsl/wavelet_trees.hpp>

#include "utility.h"

namespace project
{
class SymbolTable
{
public:

  SymbolTable () = default;
  SymbolTable (sdsl::int_vector<8> const &byte_text) noexcept;
  SymbolTable& operator= (SymbolTable &&) noexcept;

  uint64_t Serialize
  (
    std::ostream &out,
    std::shared_ptr<SpaceNode> parent = nullptr,
    std::string const name = ""
  );
  void Load (std::istream &in);

  inline uint16_t GetEffectiveAlphabetSize () const noexcept
  {
    return effective_alphabet_size_;
  }

  inline uint8_t GetEffectiveAlphabetWidth () const noexcept
  {
    return effective_alphabet_width_;
  }

  inline uint64_t ToSymbol (uint64_t const byte) const noexcept
  {
    if ((byte < std::size(alphabet_bits_)) && alphabet_bits_[byte])
    {
      return byte_to_symbol_(byte);
    }
    return 0;
  }

  inline uint64_t ToByte (uint64_t const symbol) const noexcept
  {
    if (symbol < effective_alphabet_size_)
    {
      return symbol_to_byte_(symbol + 1);
    }
    return 0;
  }

  friend std::ostream& operator<< (std::ostream &out, SymbolTable const &symbol_table);

private:

  uint16_t effective_alphabet_size_;
  uint8_t effective_alphabet_width_;
  sdsl::bit_vector alphabet_bits_;
  sdsl::bit_vector::rank_1_type byte_to_symbol_;
  sdsl::bit_vector::select_1_type symbol_to_byte_;

};

SymbolTable::SymbolTable (sdsl::int_vector<8> const &byte_text) noexcept
{
  alphabet_bits_.resize(*std::max_element(std::begin(byte_text), std::end(byte_text)) + 1);
  sdsl::util::set_to_value(alphabet_bits_, 0);
  for (auto const byte : byte_text)
  {
    alphabet_bits_[byte] = 1;
  }
  byte_to_symbol_ = decltype(byte_to_symbol_)(&alphabet_bits_);
  effective_alphabet_size_ = byte_to_symbol_(std::size(alphabet_bits_));
  effective_alphabet_width_ = sdsl::bits::hi(effective_alphabet_size_ - 1) + 1;
  symbol_to_byte_ = decltype(symbol_to_byte_)(&alphabet_bits_);
}

SymbolTable& SymbolTable::operator= (SymbolTable &&symbol_table) noexcept
{
  if (this != &symbol_table)
  {
    effective_alphabet_size_ = std::move(symbol_table.effective_alphabet_size_);
    effective_alphabet_width_ = std::move(symbol_table.effective_alphabet_width_);
    alphabet_bits_ = std::move(symbol_table.alphabet_bits_);
    byte_to_symbol_ = std::move(symbol_table.byte_to_symbol_);
    byte_to_symbol_.set_vector(&alphabet_bits_);
    symbol_to_byte_ = std::move(symbol_table.symbol_to_byte_);
    symbol_to_byte_.set_vector(&alphabet_bits_);
  }
  return *this;
}

uint64_t SymbolTable::Serialize
(
  std::ostream &out,
  std::shared_ptr<SpaceNode> parent,
  std::string const name
)
{
  uint64_t size_in_bytes {};
  if (!parent)
  {
    sdsl::write_member(effective_alphabet_size_, out);
    sdsl::write_member(effective_alphabet_width_, out);
    sdsl::serialize(alphabet_bits_, out);
    sdsl::serialize(byte_to_symbol_, out);
    sdsl::serialize(symbol_to_byte_, out);
  }
  else
  {
    auto node {std::make_shared<SpaceNode>(name)};
    node->AddLeaf("effective_alphabet_size_", sdsl::write_member(effective_alphabet_size_, out));
    node->AddLeaf("effective_alphabet_width_", sdsl::write_member(effective_alphabet_width_, out));
    node->AddLeaf("alphabet_bits_", sdsl::serialize(alphabet_bits_, out));
    node->AddLeaf("byte_to_symbol_", sdsl::serialize(byte_to_symbol_, out));
    node->AddLeaf("symbol_to_byte_", sdsl::serialize(symbol_to_byte_, out));
    parent->AddChild(node);
    size_in_bytes = node->GetSizeInBytes();
  }
  return size_in_bytes;
}

void SymbolTable::Load (std::istream &in)
{
  sdsl::read_member(effective_alphabet_size_, in);
  sdsl::read_member(effective_alphabet_width_, in);
  alphabet_bits_.load(in);
  byte_to_symbol_.load(in);
  byte_to_symbol_.set_vector(&alphabet_bits_);
  symbol_to_byte_.load(in);
  symbol_to_byte_.set_vector(&alphabet_bits_);
  return;
}

std::ostream& operator<< (std::ostream &out, SymbolTable const &symbol_table)
{
  {
    out << "value:\n";
    out << "effective_alphabet_size_:\n";
    out << static_cast<uint64_t>(symbol_table.effective_alphabet_size_) << "\n";
    out << "effective_alphabet_width_:\n";
    out << static_cast<uint64_t>(symbol_table.effective_alphabet_width_) << "\n";
    out << "byte alphabet:\n";
    for (uint16_t symbol {}; symbol != symbol_table.effective_alphabet_size_; ++symbol)
    {
      out << symbol_table.ToByte(symbol);
      out << ((symbol != (symbol_table.effective_alphabet_size_ - 1)) ? " " : "\n");
    }
  }
  {
    out << "space:\n";
    out << "effective_alphabet_size_: " << sizeof(symbol_table.effective_alphabet_size_) << "B\n";
    out << "effective_alphabet_width_: " << sizeof(symbol_table.effective_alphabet_width_) << "B\n";
    out << "alphabet_bits_: " << ProperSizeRepresentation(sdsl::size_in_bytes(symbol_table.alphabet_bits_)) << "B\n";
    out << "byte_to_symbol_: " << ProperSizeRepresentation(sdsl::size_in_bytes(symbol_table.byte_to_symbol_)) << "B\n";
    out << "symbol_to_byte_: " << ProperSizeRepresentation(sdsl::size_in_bytes(symbol_table.symbol_to_byte_)) << "B\n";
  }
  return out;
}

class Trie
{
public:

  class Node
  {
  public:
    std::map<uint8_t, std::shared_ptr<Node>> children;
    uint64_t count;
  };

  Trie (): root_ {std::make_shared<Node>()} {}

  template <typename Iterator>
  void Insert (Iterator it, Iterator last);

  inline auto GetRoot () const noexcept
  {
    return root_;
  }

  friend std::ostream& operator<< (std::ostream &out, std::pair<Trie, bool> const &pair);

private:

  std::shared_ptr<Node> root_;

};

template <typename Iterator>
void Trie::Insert (Iterator it, Iterator end)
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
  return;
}

std::ostream& operator<< (std::ostream &out, std::pair<Trie, bool> const &pair)
{
  auto const &trie {std::get<0>(pair)};
  auto const is_level_order {std::get<1>(pair)};
  std::deque<std::tuple<std::shared_ptr<Trie::Node>, uint64_t, uint64_t>> nodes;
  nodes.emplace_back(trie.root_, 0, 0);
  uint64_t size {};
  while (!nodes.empty())
  {
    ++size;
    auto node {std::get<0>(nodes.front())};
    auto label {std::get<1>(nodes.front())};
    auto depth {std::get<2>(nodes.front())};
    if (is_level_order)
    {
      nodes.pop_front();
    }
    else
    {
      node = std::get<0>(nodes.back());
      label = std::get<1>(nodes.back());
      depth = std::get<2>(nodes.back());
      nodes.pop_back();
    }
    if (depth != 0)
    {
      out << depth << ":" << label << "(" << node->count << ")\n";
    }
    if (is_level_order)
    {
      for (auto it {std::begin(node->children)}; it != std::end(node->children); ++it)
      {
        nodes.emplace_back(std::get<1>(*it), std::get<0>(*it), depth + 1);
      }
    }
    else
    {
      for (auto it {std::rbegin(node->children)}; it != std::rend(node->children); ++it)
      {
        nodes.emplace_back(std::get<1>(*it), std::get<0>(*it), depth + 1);
      }
    }
  }
  out << (size - 1) << "\n";
  return out;
}

class SubFactorTrie
{
public:

  SubFactorTrie () = default;
  SubFactorTrie (Trie const &trie) noexcept;
  SubFactorTrie (SubFactorTrie const &) = default;
  SubFactorTrie& operator= (SubFactorTrie &&) noexcept;

  uint64_t Serialize
  (
    std::ostream &out,
    std::shared_ptr<SpaceNode> parent = nullptr,
    std::string const name = ""
  );
  void Load (std::istream &in);

  template <typename Iterator>
  uint64_t Count (Iterator it, Iterator end);

  friend std::ostream& operator<< (std::ostream &out, std::pair<SubFactorTrie, bool> const &pair);

private:

  sdsl::bit_vector level_order_bits_;
  sdsl::bit_vector::select_1_type level_order_select_1_;
  sdsl::int_vector<> labels_;
  sdsl::int_vector<> counts_;

};

SubFactorTrie::SubFactorTrie (Trie const &trie) noexcept
{
  std::vector<bool> level_order_bits;
  std::vector<uint8_t> labels;
  std::vector<uint64_t> counts;
  std::deque<std::shared_ptr<Trie::Node>> nodes;
  nodes.emplace_back(trie.GetRoot());
  while (!nodes.empty())
  {
    auto node {nodes.front()};
    nodes.pop_front();
    for (auto const &pair : node->children)
    {
      auto label {std::get<0>(pair)};
      auto child {std::get<1>(pair)};
      nodes.emplace_back(child);
      level_order_bits.emplace_back(0);
      labels.emplace_back(label);
      counts.emplace_back(child->count);
    }
    level_order_bits.emplace_back(1);
  }
  {
    level_order_bits_.resize(std::size(level_order_bits));
    std::copy(std::begin(level_order_bits), std::end(level_order_bits), std::begin(level_order_bits_));
    level_order_select_1_ = decltype(level_order_select_1_)(&level_order_bits_);
  }
  {
    labels_.resize(std::size(labels));
    std::copy(std::begin(labels), std::end(labels), std::begin(labels_));
    sdsl::util::bit_compress(labels_);
  }
  {
    counts_.resize(std::size(counts));
    std::copy(std::begin(counts), std::end(counts), std::begin(counts_));
    sdsl::util::bit_compress(counts_);
  }
  return;
}

SubFactorTrie& SubFactorTrie::operator= (SubFactorTrie &&sub_factor_trie) noexcept
{
  if (this != &sub_factor_trie)
  {
    level_order_bits_ = std::move(sub_factor_trie.level_order_bits_);
    level_order_select_1_ = std::move(sub_factor_trie.level_order_select_1_);
    level_order_select_1_.set_vector(&level_order_bits_);
    labels_ = std::move(sub_factor_trie.labels_);
    counts_ = std::move(sub_factor_trie.counts_);
  }
  return *this;
}

uint64_t SubFactorTrie::Serialize
(
  std::ostream &out,
  std::shared_ptr<SpaceNode> parent,
  std::string const name
)
{
  uint64_t size_in_bytes {};
  if (!parent)
  {
    sdsl::serialize(level_order_bits_, out);
    sdsl::serialize(level_order_select_1_, out);
    sdsl::serialize(labels_, out);
    sdsl::serialize(counts_, out);
  }
  else
  {
    auto node {std::make_shared<SpaceNode>(name)};
    node->AddLeaf("level_order_bits_", sdsl::serialize(level_order_bits_, out));
    node->AddLeaf("level_order_select_1_", sdsl::serialize(level_order_select_1_, out));
    node->AddLeaf("labels_", sdsl::serialize(labels_, out));
    node->AddLeaf("counts_", sdsl::serialize(counts_, out));
    parent->AddChild(node);
    size_in_bytes = node->GetSizeInBytes();
  }
  return size_in_bytes;
}

void SubFactorTrie::Load (std::istream &in)
{
  level_order_bits_.load(in);
  level_order_select_1_.load(in);
  level_order_select_1_.set_vector(&level_order_bits_);
  labels_.load(in);
  counts_.load(in);
  return;
}

template <typename Iterator>
uint64_t SubFactorTrie::Count (Iterator it, Iterator end)
{
  uint64_t offset {};
  uint64_t end_offset {level_order_select_1_(1)};
  while (true)
  {
    while ((offset != end_offset) && (*it != labels_[offset]))
    {
      ++offset;
    }
    if (offset != end_offset)
    {
      if (it != end)
      {
        end_offset = level_order_select_1_(offset + 2) - (offset + 2) + 1;
        offset = level_order_select_1_(offset + 1) - (offset + 1) + 1;
        ++it;
      }
      else
      {
        return counts_[offset];
      }
    }
    else
    {
      return 0;
    }
  }
  return 0;
}

std::ostream& operator<< (std::ostream &out, std::pair<SubFactorTrie, bool> const &pair)
{
  auto const &trie {std::get<0>(pair)};
  auto const is_level_order {std::get<1>(pair)};
  std::deque<std::pair<uint64_t, uint64_t>> offsets;
  offsets.emplace_back(0, 0);
  while (!offsets.empty())
  {
    auto offset {std::get<0>(offsets.front())};
    auto depth {std::get<1>(offsets.front())};
    if (is_level_order)
    {
      offsets.pop_front();
    }
    else
    {
      offset = std::get<0>(offsets.back());
      depth = std::get<1>(offsets.back());
      offsets.pop_back();
    }
    uint64_t end_offset {};
    if (depth != 0)
    {
      out << depth << ":" << trie.labels_[offset] << "(" << trie.counts_[offset] << ")\n";
      end_offset = trie.level_order_select_1_(offset + 2) - (offset + 2) + 1;
      offset = trie.level_order_select_1_(offset + 1) - (offset + 1) + 1;
    }
    else
    {
      end_offset = trie.level_order_select_1_(1);
    }
    if (is_level_order)
    {
      while (offset != end_offset)
      {
        offsets.emplace_back(offset, depth + 1);
        ++offset;
      }
    }
    else
    {
      std::swap(--offset, --end_offset);
      while (offset != end_offset)
      {
        offsets.emplace_back(offset, depth + 1);
        --offset;
      }
    }
  }
  out << std::size(trie.labels_) << "\n";
  out << "space:\n";
  out << "level_order_bits_: " << ProperSizeRepresentation(sdsl::size_in_bytes(trie.level_order_bits_)) << "B\n";
  out << "level_order_select_1_: " << ProperSizeRepresentation(sdsl::size_in_bytes(trie.level_order_select_1_)) << "B\n";
  out << "labels_: " << ProperSizeRepresentation(sdsl::size_in_bytes(trie.labels_)) << "B\n";
  out << "counts_: " << ProperSizeRepresentation(sdsl::size_in_bytes(trie.counts_)) << "B\n";
  return out;
}

class GrammarSymbolTable
{
public:

  GrammarSymbolTable () = default;
  GrammarSymbolTable (std::vector<uint64_t> &sorted_factor_integers) noexcept;
  GrammarSymbolTable& operator= (GrammarSymbolTable &&) noexcept;

  uint64_t Serialize
  (
    std::ostream &out,
    std::shared_ptr<SpaceNode> parent = nullptr,
    std::string const name = ""
  );
  void Load (std::istream &in);

  uint64_t operator[] (uint64_t const factor_integer) const noexcept;

  std::pair<uint64_t, uint64_t> SymbolRange
  (
    uint64_t const smallest_factor_integer,
    uint64_t const largest_factor_integer
  ) const noexcept;

  friend std::ostream& operator<< (std::ostream &out, GrammarSymbolTable const &symbol_table);

private:

  sdsl::sd_vector<> factor_integer_bits_;
  sdsl::sd_vector<>::rank_1_type factor_integer_rank_1_;
  sdsl::sd_vector<>::select_1_type factor_integer_select_1_;

};

GrammarSymbolTable::GrammarSymbolTable (std::vector<uint64_t> &sorted_factor_integers) noexcept
{
  factor_integer_bits_ = decltype(factor_integer_bits_)(std::begin(sorted_factor_integers), std::end(sorted_factor_integers));
  factor_integer_rank_1_.set_vector(&factor_integer_bits_);
}

GrammarSymbolTable& GrammarSymbolTable::operator= (GrammarSymbolTable &&grammar_symbol_table) noexcept
{
  if (this != &grammar_symbol_table)
  {
    factor_integer_bits_ = std::move(grammar_symbol_table.factor_integer_bits_);
    factor_integer_rank_1_.set_vector(&factor_integer_bits_);
  }
  return *this;
}

uint64_t GrammarSymbolTable::Serialize
(
  std::ostream &out,
  std::shared_ptr<SpaceNode> parent,
  std::string const name
)
{
  uint64_t size_in_bytes {};
  if (!parent)
  {
    sdsl::serialize(factor_integer_bits_, out);
    sdsl::serialize(factor_integer_rank_1_, out);
  }
  else
  {
    auto node {std::make_shared<SpaceNode>(name)};
    node->AddLeaf("factor_integer_bits_", sdsl::serialize(factor_integer_bits_, out));
    node->AddLeaf("factor_integer_rank_1_", sdsl::serialize(factor_integer_rank_1_, out));
    parent->AddChild(node);
    size_in_bytes = node->GetSizeInBytes();
  }
  return size_in_bytes;
}

void GrammarSymbolTable::Load (std::istream &in)
{
  factor_integer_bits_.load(in);
  factor_integer_rank_1_.load(in);
  factor_integer_rank_1_.set_vector(&factor_integer_bits_);
  return;
}

uint64_t GrammarSymbolTable::operator[] (uint64_t const factor_integer) const noexcept
{
  if ((factor_integer < std::size(factor_integer_bits_)) && factor_integer_bits_[factor_integer])
  {
    return factor_integer_rank_1_(factor_integer);
  }
  return 0;
}

std::pair<uint64_t, uint64_t> GrammarSymbolTable::SymbolRange
(
  uint64_t const smallest_factor_integer,
  uint64_t const largest_factor_integer
) const noexcept
{
  std::pair<uint64_t, uint64_t> symbol_range {1, 0};
  if (smallest_factor_integer < std::size(factor_integer_bits_))
  {
    std::get<0>(symbol_range) = factor_integer_rank_1_(smallest_factor_integer);
    if (factor_integer_bits_[smallest_factor_integer] == 0)
    {
      ++std::get<0>(symbol_range);
    }
    if (largest_factor_integer < std::size(factor_integer_bits_))
    {
      std::get<1>(symbol_range) = factor_integer_rank_1_(largest_factor_integer);
    }
    else
    {
      std::get<1>(symbol_range) = factor_integer_rank_1_(std::size(factor_integer_bits_));
    }
  }
  return symbol_range;
}

std::ostream& operator<< (std::ostream &out, GrammarSymbolTable const &grammar_symbol_table)
{
  out << "space:\n";
  out << "factor_integer_bits_: " << ProperSizeRepresentation(sdsl::size_in_bytes(grammar_symbol_table.factor_integer_bits_)) << "B\n";
  return out;
}

template <uint8_t max_factor_size>
class Index;

template <uint8_t max_factor_size>
std::ostream& operator<< (std::ostream &out, Index<max_factor_size> const &index);

template <uint8_t max_factor_size = 4>
class Index
{
public:

  static constexpr uint8_t kMaxFactorSize {max_factor_size};
  static constexpr uint8_t kS {1};
  static constexpr uint8_t kL {0};

  Index () = default;
  Index (std::filesystem::path const &byte_text_path);

  uint64_t Serialize
  (
    std::filesystem::path const &index_path,
    std::shared_ptr<SpaceNode> root = nullptr
  );
  void Load (std::filesystem::path const &index_path);

  template <typename Iterator>
  uint64_t Count (Iterator begin, Iterator end);

  friend std::ostream& operator<< <max_factor_size> (std::ostream &out, Index<max_factor_size> const &index);

private:

  SymbolTable symbol_table_;
  SubFactorTrie sub_factor_trie_;
  GrammarSymbolTable lex_symbol_table_;
  GrammarSymbolTable colex_symbol_table_;
  sdsl::int_vector<> colex_to_lex_;
  sdsl::int_vector<> lex_symbol_bucket_offsets_;
  sdsl::wt_rlmn
  <
    sdsl::sd_vector<>,
    typename sdsl::sd_vector<>::rank_1_type,
    typename sdsl::sd_vector<>::select_1_type,
    sdsl::wt_ap<>
  >
  lex_bwt_;

  void CalculateSymbolTableAndText
  (
    std::filesystem::path const &byte_text_path,
    sdsl::int_vector<> &text
  );

  void CalculateSubFactorTrie (sdsl::int_vector<> const &text, Trie &trie);
  template <typename Iterator>
  uint64_t SymbolsToInteger (Iterator it, Iterator end, int8_t const step = 1);
  std::deque<uint64_t> IntegerToSymbols (uint64_t integer);
  void CalculateLexTextSizeAndLexFactorIntegersAndTemporaryColexToLex
  (
    sdsl::int_vector<> const &text,
    uint64_t &lex_text_size,
    std::set<uint64_t> &lex_factor_integers,
    std::map<uint64_t, uint64_t> &temporary_colex_to_lex
  );
  void CalculateLexSymbolTable (std::set<uint64_t> const &lex_factor_integers);
  void CalculateLexText (sdsl::int_vector<> const &text, sdsl::int_vector<> &lex_text);
  void CalculateColexSymbolTable (std::map<uint64_t, uint64_t> const &temporary_colex_to_lex);
  void CalculateColexToLex (std::map<uint64_t, uint64_t> const &temporary_colex_to_lex);
  void CalculateLexSymbolBucketOffsets
  (
    uint64_t const lex_text_alphabet_size,
    sdsl::int_vector<> const &lex_text
  );
  void CalculateLexBwt (sdsl::int_vector<> const &lex_text);

  template <typename Range> // [,]
  inline bool IsNotEmptyRange (Range const &range) const noexcept
  {
    return (std::get<0>(range) <= std::get<1>(range));
  }

  template <typename Range> // [,]
  inline uint64_t RangeSize (Range const &range) const noexcept
  {
    if (std::get<0>(range) <= std::get<1>(range))
    {
      return (std::get<1>(range) - std::get<0>(range) + 1);
    }
    return 0;
  }

  template <typename Iterator>
  void CalculateSlFactor (Iterator const rend, Iterator &rfirst, Iterator &rlast);

  template <typename Iterator>
  Iterator ClassifyPatternSuffix (Iterator rbegin, Iterator rend);

  template <typename Iterator, typename Range>
  Iterator BackwardSearchPatternSuffix
  (
    Iterator rbegin,
    Iterator rend,
    Range &pattern_range,
    Range &pattern_range_ls
  );

  template <typename Iterator, typename Range>
  void BackwardSearchPatternSuffixSlFactor (Iterator rbegin, Iterator rend, Range &pattern_range);

  template <typename Iterator>
  std::pair<uint64_t, uint64_t> CalculateSymbolRange
  (
    Iterator begin,
    Iterator end,
    uint8_t const step = 1
  );

  template <typename Iterator, typename Range>
  void BackwardSearchPatternInfixFactors (Iterator rbegin, Iterator rend, Range &pattern_range);

  template <typename Range>
  uint64_t RangeCount
  (
    Range const pattern_range,
    Range const colex_symbol_range
  );

  template <typename Iterator, typename Range>
  void BackwardSearchPatternInfixSlFactor (Iterator rbegin, Iterator rend, Range &pattern_range);

  template <typename Iterator, typename Range>
  void BackwardSearchPatternPrefixSlFactor (Iterator rbegin, Iterator rend, Range &pattern_range);

  template <typename Iterator, typename Range>
  void CountWithinSlFactor (Iterator rbegin, Iterator rend, Range &pattern_range);

  template <typename Iterator, typename Range>
  void BackwardSearchPatternInfix
  (
    Iterator rbegin,
    Iterator rend,
    Range &pattern_range,
    Range &pattern_range_ls
  );

  template <typename Iterator, typename Range>
  void BackwardSearchPatternPrefix
  (
    Iterator rbegin,
    Iterator rend,
    Range &pattern_range,
    Range &pattern_range_ls
  );

};

template <uint8_t max_factor_size>
Index<max_factor_size>::Index (std::filesystem::path const &byte_text_path)
{
  std::cout << "construct index of " << std::filesystem::canonical(byte_text_path) << "\n";
  sdsl::int_vector<> text;
  CalculateSymbolTableAndText(byte_text_path, text);
  {
    Trie trie;
    CalculateSubFactorTrie(text, trie);
    // std::cout << std::make_pair(trie, /*is_level_order=*/true);
    sub_factor_trie_ = decltype(sub_factor_trie_)(trie);
    // std::cout << std::make_pair(sub_factor_trie_, /*is_level_order=*/true);
  }
  uint64_t lex_text_size {};
  std::set<uint64_t> lex_factor_integers;
  std::map<uint64_t, uint64_t> temporary_colex_to_lex;
  {
    CalculateLexTextSizeAndLexFactorIntegersAndTemporaryColexToLex
    (
      text,
      lex_text_size,
      lex_factor_integers,
      temporary_colex_to_lex
    );
    // std::cout << "lex_text_size:\n" << lex_text_size << "\n";
    // std::cout << "|lex_factor_integers|: " << std::size(lex_factor_integers) << "\n";
    // std::cout << "lex_factor_integers:\n";
    // for (auto const &factor : lex_factor_integers)
    // {
    //   Print(IntegerToSymbols(factor), std::cout);
    // }
    // std::cout << "temporary_colex_to_lex:\n";
    // for (auto const &pair : temporary_colex_to_lex)
    // {
    //   Print(IntegerToSymbols(std::get<0>(pair)), std::cout, 1, " ", "");
    //   std::cout << ":";
    //   Print(IntegerToSymbols(std::get<1>(pair)), std::cout, 1, " ", "");
    //   std::cout << "\n";
    // }
  }
  uint64_t lex_text_alphabet_size {std::size(lex_factor_integers)};
  uint64_t lex_text_width {sdsl::bits::hi(lex_text_alphabet_size - 1) + 1};
  {
    CalculateLexSymbolTable(lex_factor_integers);
    // std::cout << lex_symbol_table_;
  }
  sdsl::int_vector<> lex_text(lex_text_size, 0, lex_text_width);
  {
    CalculateLexText(text, lex_text);
    // Print(lex_text, std::cout);
  }
  {
    CalculateColexSymbolTable(temporary_colex_to_lex);
    // std::cout << colex_symbol_table_;
  }
  {
    CalculateColexToLex(temporary_colex_to_lex);
    // Print(colex_to_lex_, std::cout);
  }
  {
    CalculateLexSymbolBucketOffsets(lex_text_alphabet_size, lex_text);
    // Print(lex_symbol_bucket_offsets_, std::cout);
  }
  {
    CalculateLexBwt(lex_text);
    // std::cout << lex_bwt_ << "\n";
  }
}

template <uint8_t max_factor_size>
uint64_t Index<max_factor_size>::Serialize
(
  std::filesystem::path const &index_path,
  std::shared_ptr<SpaceNode> root
)
{
  std::fstream index_file {index_path, std::ios_base::out | std::ios_base::trunc};
  std::cout << "serialize index to " << std::filesystem::canonical(index_path) << "\n";
  uint64_t size_in_bytes {};
  if (!root)
  {
    sdsl::write_member(kMaxFactorSize, index_file);
    symbol_table_.Serialize(index_file);
    sub_factor_trie_.Serialize(index_file);
    lex_symbol_table_.Serialize(index_file);
    colex_symbol_table_.Serialize(index_file);;
    sdsl::serialize(colex_to_lex_, index_file);
    sdsl::serialize(lex_symbol_bucket_offsets_, index_file);
    sdsl::serialize(lex_bwt_, index_file);
  }
  else
  {
    root->AddLeaf("kMaxFactorSize", sdsl::write_member(kMaxFactorSize, index_file));
    root->AccumalateSizeInBytes(symbol_table_.Serialize(index_file, root, "symbol_table_"));
    root->AccumalateSizeInBytes(sub_factor_trie_.Serialize(index_file, root, "sub_factor_trie_"));
    root->AccumalateSizeInBytes(lex_symbol_table_.Serialize(index_file, root, "lex_symbol_table_"));
    root->AccumalateSizeInBytes(colex_symbol_table_.Serialize(index_file, root, "colex_symbol_table_"));
    root->AddLeaf("colex_to_lex_", sdsl::serialize(colex_to_lex_, index_file));
    root->AddLeaf("lex_symbol_bucket_offsets_", sdsl::serialize(lex_symbol_bucket_offsets_, index_file));
    root->AddLeaf("lex_bwt_", sdsl::serialize(lex_bwt_, index_file));
    size_in_bytes = root->GetSizeInBytes();
  }
  return size_in_bytes;
}

template <uint8_t max_factor_size>
template <typename Iterator>
uint64_t Index<max_factor_size>::Count (Iterator begin, Iterator end)
{
  // Print(std::cout, begin, end);
  std::pair<uint64_t, uint64_t> pattern_range {1, 0};
  std::pair<uint64_t, uint64_t> pattern_range_ls {1, 0};
  auto rbegin {std::prev(end)};
  auto rend {std::prev(begin)};
  auto rfirst {rbegin};
  auto rlast {rbegin};
  rlast = BackwardSearchPatternSuffix(rbegin, rend, pattern_range, pattern_range_ls);
  if (rlast != rend)
  {
    while (IsNotEmptyRange(pattern_range) || IsNotEmptyRange(pattern_range_ls))
    {
      CalculateSlFactor(rend, rfirst, rlast);
      if (rlast != rend)
      {
        BackwardSearchPatternInfix(rfirst, rlast, pattern_range, pattern_range_ls);
      }
      else
      {
        BackwardSearchPatternPrefix(rfirst, rlast, pattern_range, pattern_range_ls);
        break;
      }
    }
  }
  return (RangeSize(pattern_range) + RangeSize(pattern_range_ls));
}

template <uint8_t max_factor_size>
void Index<max_factor_size>::Load (std::filesystem::path const &index_path)
{
  std::ifstream index_file {index_path};
  std::cout << "load index from " << std::filesystem::canonical(index_path) << "\n";
  uint8_t loaded_max_factor_size {};
  sdsl::read_member(loaded_max_factor_size, index_file);
  if (loaded_max_factor_size != kMaxFactorSize)
  {
    throw std::runtime_error("wrong max_factor_size");
  }
  symbol_table_.Load(index_file);
  sub_factor_trie_.Load(index_file);
  lex_symbol_table_.Load(index_file);
  colex_symbol_table_.Load(index_file);
  colex_to_lex_.load(index_file);
  lex_symbol_bucket_offsets_.load(index_file);
  lex_bwt_.load(index_file);
  return;
}

template <uint8_t max_factor_size>
std::ostream& operator<< (std::ostream &out, Index<max_factor_size> const &index)
{
  out << index.symbol_table_;
  out << std::make_pair(index.sub_factor_trie_, true);
  out << index.lex_symbol_table_;
  out << index.colex_symbol_table_;
  out << index.colex_to_lex_ << "\n";
  out << index.lex_symbol_bucket_offsets_ << "\n";
  out << index.lex_bwt_ << "\n";
  return out;
}

template <uint8_t max_factor_size>
void Index<max_factor_size>::CalculateSymbolTableAndText
(
  std::filesystem::path const &byte_text_path,
  sdsl::int_vector<> &text
)
{
  sdsl::int_vector<8> byte_text;
  sdsl::load_vector_from_file(byte_text, byte_text_path);
  {
    for (auto byte : byte_text)
    {
      if (byte == 0)
      {
        throw std::runtime_error("byte_text contains 0");
      }
    }
    sdsl::append_zero_symbol(byte_text);
  }
  symbol_table_ = decltype(symbol_table_)(byte_text);
  // std::cout << symbol_table_;
  text.width(symbol_table_.GetEffectiveAlphabetWidth());
  text.resize(std::size(byte_text));
  std::transform
  (
    std::begin(byte_text),
    std::end(byte_text),
    std::begin(text),
    [&] (auto const byte)
    {
      return symbol_table_.ToSymbol(byte);
    }
  );
  // Print(text, std::cout);
}

template <uint8_t max_factor_size>
void Index<max_factor_size>::CalculateSubFactorTrie (sdsl::int_vector<> const &text, Trie &trie)
{
  auto text_prev_begin {std::prev(std::begin(text))};
  auto text_it {std::prev(std::end(text), 3)};
  auto next_symbol {*std::prev(std::end(text), 2)};
  uint8_t sl_type {};
  uint8_t next_sl_type {Index::kL};
  auto sl_factor_end {std::prev(std::end(text))};
  while (text_it != text_prev_begin)
  {
    if (*text_it == next_symbol)
    {
      sl_type = next_sl_type;
    }
    else if (*text_it < next_symbol)
    {
      sl_type = Index::kS;
    }
    else
    {
      sl_type = Index::kL;
    }
    if ((sl_type == Index::kL) && (next_sl_type == Index::kS))
    {
      auto sl_factor_it {std::next(text_it)};
      while (std::distance(sl_factor_it, sl_factor_end) >= (Index::kMaxFactorSize - 1))
      {
        auto it {sl_factor_it};
        auto last {std::next(it, Index::kMaxFactorSize - 1)};
        while (it != last)
        {
          trie.Insert(it, last);
          ++it;
        }
        ++sl_factor_it;
      }
      while (sl_factor_it != sl_factor_end)
      {
        trie.Insert(sl_factor_it, sl_factor_end);
        ++sl_factor_it;
      }
      sl_factor_end = std::next(text_it);
    }
    next_symbol = *text_it--;
    next_sl_type = sl_type;
  }
  {
    auto sl_factor_it {std::next(text_prev_begin)};
    while (std::distance(sl_factor_it, sl_factor_end) >= (Index::kMaxFactorSize - 1))
    {
      auto it {sl_factor_it};
      auto last {std::next(it, Index::kMaxFactorSize - 1)};
      while (it != last)
      {
        trie.Insert(it, last);
        ++it;
      }
      ++sl_factor_it;
    }
    while (sl_factor_it != sl_factor_end)
    {
      trie.Insert(sl_factor_it, sl_factor_end);
      ++sl_factor_it;
    }
  }
  return;
}

template <uint8_t max_factor_size>
template <typename Iterator>
uint64_t Index<max_factor_size>::SymbolsToInteger (Iterator it, Iterator end, int8_t const step)
{
  uint64_t result {};
  for (auto k {Index::kMaxFactorSize}; (it != end) && (k != 0); it += step, --k)
  {
    result += *it * (1ULL << (symbol_table_.GetEffectiveAlphabetWidth() * (k - 1)));
  }
  return result;
}

template <uint8_t max_factor_size>
std::deque<uint64_t> Index<max_factor_size>::IntegerToSymbols (uint64_t integer)
{
  std::deque<uint64_t> symbols;
  auto const base {1ULL << symbol_table_.GetEffectiveAlphabetWidth()};
  if (integer == 0)
  {
    symbols.emplace_back(0);
  }
  else
  {
    while (integer != 0)
    {
      symbols.emplace_front(integer % base);
      integer /= base;
    }
  }
  return symbols;
}

template <uint8_t max_factor_size>
void Index<max_factor_size>::CalculateLexTextSizeAndLexFactorIntegersAndTemporaryColexToLex
(
  sdsl::int_vector<> const &text,
  uint64_t &lex_text_size,
  std::set<uint64_t> &lex_factor_integers,
  std::map<uint64_t, uint64_t> &temporary_colex_to_lex
)
{
  {
    ++lex_text_size;
    lex_factor_integers.insert(0);
    temporary_colex_to_lex[0] = 0;
  }
  auto text_prev_begin {std::prev(std::begin(text))};
  auto text_it {std::prev(std::end(text), 3)};
  auto next_symbol {*std::prev(std::end(text), 2)};
  uint8_t sl_type {};
  uint8_t next_sl_type {Index::kL};
  auto sl_factor_end {std::prev(std::end(text))};
  uint64_t lex_factor_integer {};
  uint64_t colex_factor_integer {};
  while (text_it != text_prev_begin)
  {
    if (*text_it == next_symbol)
    {
      sl_type = next_sl_type;
    }
    else if (*text_it < next_symbol)
    {
      sl_type = Index::kS;
    }
    else
    {
      sl_type = Index::kL;
    }
    if ((sl_type == Index::kL) && (next_sl_type == Index::kS))
    {
      auto factor_it {std::next(text_it)};
      while (std::distance(factor_it, sl_factor_end) >= Index::kMaxFactorSize)
      {
        auto factor_end {std::next(factor_it, Index::kMaxFactorSize)};
        ++lex_text_size;
        lex_factor_integer = SymbolsToInteger(factor_it, factor_end);
        colex_factor_integer = SymbolsToInteger(std::prev(factor_end), std::prev(factor_it), -1);
        lex_factor_integers.insert(lex_factor_integer);
        temporary_colex_to_lex[colex_factor_integer] = lex_factor_integer;
        factor_it = factor_end;
      }
      if (std::distance(factor_it, sl_factor_end) != 0)
      {
        ++lex_text_size;
        lex_factor_integer = SymbolsToInteger(factor_it, sl_factor_end);
        colex_factor_integer = SymbolsToInteger(std::prev(sl_factor_end), std::prev(factor_it), -1);
        lex_factor_integers.insert(lex_factor_integer);
        temporary_colex_to_lex[colex_factor_integer] = lex_factor_integer;
      }
      sl_factor_end = std::next(text_it);
    }
    next_symbol = *text_it--;
    next_sl_type = sl_type;
  }
  ++lex_text_size;
  lex_factor_integer = SymbolsToInteger(std::next(text_prev_begin), sl_factor_end);
  colex_factor_integer = SymbolsToInteger(std::prev(sl_factor_end), text_prev_begin, -1);
  lex_factor_integers.insert(lex_factor_integer);
  temporary_colex_to_lex[colex_factor_integer] = lex_factor_integer;
  return;
}

template <uint8_t max_factor_size>
void Index<max_factor_size>::CalculateLexSymbolTable (std::set<uint64_t> const &lex_factor_integers)
{
  std::vector<uint64_t> sorted_lex_factor_integers(std::begin(lex_factor_integers), std::end(lex_factor_integers));
  lex_symbol_table_ = decltype(lex_symbol_table_)(sorted_lex_factor_integers);
  return;
}

template <uint8_t max_factor_size>
void Index<max_factor_size>::CalculateLexText
(
  sdsl::int_vector<> const &text,
  sdsl::int_vector<> &lex_text
)
{
  auto text_prev_begin {std::prev(std::begin(text))};
  auto text_it {std::prev(std::end(text), 3)};
  auto next_symbol {*std::prev(std::end(text), 2)};
  uint8_t sl_type {};
  uint8_t next_sl_type {Index::kL};
  auto sl_factor_end {std::prev(std::end(text))};
  auto lex_text_it {std::prev(std::end(lex_text), 2)};
  std::vector<uint64_t> sl_factor;
  while (text_it != text_prev_begin)
  {
    if (*text_it == next_symbol)
    {
      sl_type = next_sl_type;
    }
    else if (*text_it < next_symbol)
    {
      sl_type = Index::kS;
    }
    else
    {
      sl_type = Index::kL;
    }
    if ((sl_type == Index::kL) && (next_sl_type == Index::kS))
    {
      auto factor_it {std::next(text_it)};
      while (std::distance(factor_it, sl_factor_end) >= Index::kMaxFactorSize)
      {
        auto factor_end {std::next(factor_it, Index::kMaxFactorSize)};
        sl_factor.emplace_back(lex_symbol_table_[SymbolsToInteger(factor_it, factor_end)]);
        factor_it = factor_end;
      }
      if (std::distance(factor_it, sl_factor_end) != 0)
      {
        sl_factor.emplace_back(lex_symbol_table_[SymbolsToInteger(factor_it, sl_factor_end)]);
      }
      sl_factor_end = std::next(text_it);
      for (auto it {std::rbegin(sl_factor)}; it != std::rend(sl_factor); ++it)
      {
        *lex_text_it-- = *it;
      }
      sl_factor.clear();
    }
    next_symbol = *text_it--;
    next_sl_type = sl_type;
  }
  *lex_text_it = lex_symbol_table_[SymbolsToInteger(std::next(text_prev_begin), sl_factor_end)];
  *std::prev(std::end(lex_text)) = 0;
  return;
}

template <uint8_t max_factor_size>
void Index<max_factor_size>::CalculateColexSymbolTable (std::map<uint64_t, uint64_t> const &temporary_colex_to_lex)
{
  std::vector<uint64_t> sorted_colex_factor_integers(std::size(temporary_colex_to_lex));
  auto it {std::begin(sorted_colex_factor_integers)};
  for (auto const &pair : temporary_colex_to_lex)
  {
    *it++ = std::get<0>(pair);
  }
  colex_symbol_table_ = decltype(colex_symbol_table_)(sorted_colex_factor_integers);
  return;
}

template <uint8_t max_factor_size>
void Index<max_factor_size>::CalculateColexToLex (std::map<uint64_t, uint64_t> const &temporary_colex_to_lex)
{
  colex_to_lex_.width(sdsl::bits::hi(std::size(temporary_colex_to_lex) - 1) + 1);
  colex_to_lex_.resize(std::size(temporary_colex_to_lex));
  auto it {std::begin(colex_to_lex_)};
  for (auto const &pair : temporary_colex_to_lex)
  {
    *it++ = lex_symbol_table_[std::get<1>(pair)];
  }
  return;
}

template <uint8_t max_factor_size>
void Index<max_factor_size>::CalculateLexSymbolBucketOffsets
(
  uint64_t const lex_text_alphabet_size,
  sdsl::int_vector<> const &lex_text
)
{
  lex_symbol_bucket_offsets_.width(sdsl::bits::hi(std::size(lex_text)) + 1);
  lex_symbol_bucket_offsets_.resize(lex_text_alphabet_size + 1);
  sdsl::util::set_to_value(lex_symbol_bucket_offsets_, 0);
  for (auto const lex_symbol : lex_text)
  {
    ++(lex_symbol_bucket_offsets_[lex_symbol]);
  }
  std::partial_sum
  (
    std::begin(lex_symbol_bucket_offsets_),
    std::end(lex_symbol_bucket_offsets_),
    std::begin(lex_symbol_bucket_offsets_)
  );
  auto it {std::prev(std::end(lex_symbol_bucket_offsets_))};
  auto begin {std::begin(lex_symbol_bucket_offsets_)};
  while (it != begin)
  {
    *it-- = *std::prev(it);
  }
  *begin = 0;
  return;
}

template <uint8_t max_factor_size>
void Index<max_factor_size>::CalculateLexBwt (sdsl::int_vector<> const &lex_text)
{
  sdsl::int_vector<> buffer;
  sdsl::qsufsort::construct_sa(buffer, lex_text);
  for (auto it {std::begin(buffer)}; it != std::end(buffer); ++it)
  {
    if (*it != 0)
    {
      *it = lex_text[*it - 1];
    }
  }
  sdsl::construct_im(lex_bwt_, buffer);
  return;
}

template <uint8_t max_factor_size>
template <typename Iterator, typename Range>
Iterator Index<max_factor_size>::BackwardSearchPatternSuffix
(
  Iterator rbegin,
  Iterator rend,
  Range &pattern_range,
  Range &pattern_range_ls
)
{
  auto rfirst {rbegin};
  auto rlast {ClassifyPatternSuffix(rbegin, rend)};
  if (rlast != rend)
  {
    if (rlast != rbegin)
    {
      BackwardSearchPatternSuffixSlFactor(rfirst, rlast, pattern_range_ls);
      CalculateSlFactor(rend, rfirst, rlast);
      if (IsNotEmptyRange(pattern_range_ls))
      {
        if (rlast != rend)
        {
          BackwardSearchPatternInfixSlFactor(rfirst, rlast, pattern_range_ls);
        }
        else
        {
          BackwardSearchPatternPrefixSlFactor(rfirst, rlast, pattern_range_ls);
        }
      }
    }
    else
    {
      CalculateSlFactor(rend, rfirst, rlast);
    }
    if (rlast != rend)
    {
      BackwardSearchPatternSuffixSlFactor(rbegin, rlast, pattern_range);
    }
    else
    {
      CountWithinSlFactor(rbegin, rend, pattern_range);
    }
  }
  else
  {
    CountWithinSlFactor(rbegin, rend, pattern_range);
  }
  return rlast;
}

template <uint8_t max_factor_size>
template <typename Iterator>
Iterator Index<max_factor_size>::ClassifyPatternSuffix (Iterator rbegin, Iterator rend)
{
  auto it {rbegin};
  while ((std::prev(it) != rend) && (*std::prev(it) == *it))
  {
    --it;
  }
  if ((std::prev(it) != rend) && (*std::prev(it) < *it))
  {
    return rbegin;
  }
  return std::prev(it); // LS
}

template <uint8_t max_factor_size>
template <typename Iterator>
void Index<max_factor_size>::CalculateSlFactor
(
  Iterator const rend,
  Iterator &rfirst,
  Iterator &rlast
)
{
  auto prev_sl_type {Index::kL};
  rfirst = rlast--;
  while ((rlast != rend) && !((prev_sl_type == Index::kS) && (*rlast > *std::next(rlast))))
  {
    if((prev_sl_type == Index::kL) && (*rlast < *std::next(rlast)))
    {
      prev_sl_type = Index::kS;
    }
    --rlast;
  }
  return;
}

template <uint8_t max_factor_size>
template <typename Iterator, typename Range>
void Index<max_factor_size>::BackwardSearchPatternSuffixSlFactor
(
  Iterator rbegin,
  Iterator rend,
  Range &pattern_range
)
{
  auto rfirst {rbegin};
  auto rlast {std::prev(rbegin, std::distance(rend, rbegin) % Index::kMaxFactorSize)};
  if (rlast == rfirst)
  {
    rlast = std::prev(rbegin, Index::kMaxFactorSize);
  }
  auto lex_symbol_range {CalculateSymbolRange(std::next(rlast), std::next(rfirst))};
  if (IsNotEmptyRange(lex_symbol_range))
  {
    pattern_range =
    {
      lex_symbol_bucket_offsets_[std::get<0>(lex_symbol_range)],
      lex_symbol_bucket_offsets_[std::get<1>(lex_symbol_range) + 1] - 1
    };
    if (rlast != rend)
    {
      BackwardSearchPatternInfixFactors(rlast, rend, pattern_range);
    }
  }
  return;
}

template <uint8_t max_factor_size>
template <typename Iterator>
std::pair<uint64_t, uint64_t> Index<max_factor_size>::CalculateSymbolRange
(
  Iterator begin,
  Iterator end,
  uint8_t const step
)
{
  uint64_t factor_integer {};
  auto it {begin};
  auto k {Index::kMaxFactorSize};
  while (it != std::prev(end))
  {
    factor_integer += *it * (1ULL << (symbol_table_.GetEffectiveAlphabetWidth() * ((k--) - 1)));
    it += step;
  }
  auto base {1ULL << (symbol_table_.GetEffectiveAlphabetWidth() * (k - 1))};
  factor_integer += *it * base;
  if (step == -1)
  {
    return colex_symbol_table_.SymbolRange(factor_integer, factor_integer + base - 1);
  }
  return lex_symbol_table_.SymbolRange(factor_integer, factor_integer + base - 1);
}

template <uint8_t max_factor_size>
template <typename Iterator, typename Range>
void Index<max_factor_size>::BackwardSearchPatternInfixFactors
(
  Iterator rbegin,
  Iterator rend,
  Range &pattern_range
)
{
  auto rfirst {rbegin};
  auto rlast {rbegin};
  while ((rfirst != rend) && IsNotEmptyRange(pattern_range))
  {
    rlast = std::prev(rfirst, Index::kMaxFactorSize);
    auto lex_symbol {lex_symbol_table_[SymbolsToInteger(std::next(rlast), std::next(rfirst))]};
    if (lex_symbol != 0)
    {
      auto begin_offset {lex_symbol_bucket_offsets_[lex_symbol]};
      pattern_range =
      {
        begin_offset + lex_bwt_.rank(std::get<0>(pattern_range), lex_symbol),
        begin_offset + lex_bwt_.rank(std::get<1>(pattern_range) + 1, lex_symbol) - 1
      };
    }
    else
    {
      pattern_range = {1, 0};
    }
    rfirst = rlast;
  }
  return;
}

template <uint8_t max_factor_size>
template <typename Iterator, typename Range>
void Index<max_factor_size>::BackwardSearchPatternInfixSlFactor
(
  Iterator rbegin,
  Iterator rend,
  Range &pattern_range
)
{
  auto rfirst {rbegin};
  auto rlast {std::prev(rbegin, std::distance(rend, rbegin) % Index::kMaxFactorSize)};
  if (rfirst != rlast)
  {
    auto lex_symbol {lex_symbol_table_[SymbolsToInteger(std::next(rlast), std::next(rfirst))]};
    if (lex_symbol != 0)
    {
      auto begin_offset {lex_symbol_bucket_offsets_[lex_symbol]};
      pattern_range =
      {
        begin_offset + lex_bwt_.rank(std::get<0>(pattern_range), lex_symbol),
        begin_offset + lex_bwt_.rank(std::get<1>(pattern_range) + 1, lex_symbol) - 1
      };
    }
    else
    {
      pattern_range = {1, 0};
    }
  }
  if (IsNotEmptyRange(pattern_range))
  {
    BackwardSearchPatternInfixFactors(rlast, rend, pattern_range);
  }
  return;
}

template <uint8_t max_factor_size>
template <typename Iterator, typename Range>
void Index<max_factor_size>::BackwardSearchPatternPrefixSlFactor
(
  Iterator rbegin,
  Iterator rend,
  Range &pattern_range
)
{
  uint64_t count {};
  auto sl_factor_size {std::distance(rend, rbegin)};
  for (int64_t k {1}; k <= Index::kMaxFactorSize; ++k)
  {
    auto temporary_pattern_range {pattern_range};
    if (k != sl_factor_size)
    {
      auto rfirst {rbegin};
      auto rlast {std::prev(rbegin, k)};
      auto lex_symbol {lex_symbol_table_[SymbolsToInteger(std::next(rlast), std::next(rfirst))]};
      if (lex_symbol != 0)
      {
        auto begin_offset {lex_symbol_bucket_offsets_[lex_symbol]};
        pattern_range =
        {
          begin_offset + lex_bwt_.rank(std::get<0>(temporary_pattern_range), lex_symbol),
          begin_offset + lex_bwt_.rank(std::get<1>(temporary_pattern_range) + 1, lex_symbol) - 1
        };
      }
      else
      {
        pattern_range = {1, 0};
      }
      if (IsNotEmptyRange(temporary_pattern_range))
      {
        if (std::distance(rend, rlast) > Index::kMaxFactorSize)
        {
          rfirst = rlast;
          rlast = std::next(rend, (sl_factor_size - k) % Index::kMaxFactorSize);
          BackwardSearchPatternInfixFactors(rfirst, rlast, temporary_pattern_range);
        }
        if (IsNotEmptyRange(temporary_pattern_range))
        {
          auto colex_symbol_range {CalculateSymbolRange(rlast, rend, -1)};
          if (IsNotEmptyRange(colex_symbol_range))
          {
            count += RangeCount(temporary_pattern_range, colex_symbol_range);
          }
        }
      }
    }
    else
    {
      auto colex_symbol_range {CalculateSymbolRange(rbegin, rend, -1)};
      if (IsNotEmptyRange(colex_symbol_range))
      {
        count += RangeCount(temporary_pattern_range, colex_symbol_range);
      }
      break;
    }
  }
  pattern_range = {1, count};
  return;
}

template <uint8_t max_factor_size>
template <typename Range>
uint64_t Index<max_factor_size>::RangeCount
(
  Range const pattern_range,
  Range const colex_symbol_range
)
{
  uint64_t count {};
  for
  (
    auto colex_symbol {std::get<0>(colex_symbol_range)};
    colex_symbol != std::get<1>(colex_symbol_range);
    ++colex_symbol
  )
  {
    count +=
    (
      lex_bwt_.rank(std::get<1>(pattern_range) + 1, colex_to_lex_[colex_symbol])
      - lex_bwt_.rank(std::get<0>(pattern_range), colex_to_lex_[colex_symbol])
    );
  }
  return count;
}

template <uint8_t max_factor_size>
template <typename Iterator, typename Range>
void Index<max_factor_size>::CountWithinSlFactor
(
  Iterator rbegin,
  Iterator rend,
  Range &pattern_range
)
{
  if (std::distance(rend, rbegin) >= Index::kMaxFactorSize)
  {
    uint64_t count {};
    for (uint64_t k {1}; k <= Index::kMaxFactorSize; ++k)
    {
      auto rfirst {rbegin};
      auto rlast {std::prev(rbegin, k)};
      auto temporary_pattern_range {pattern_range};
      auto lex_symbol_range {CalculateSymbolRange(std::next(rlast), std::next(rfirst))};
      if (IsNotEmptyRange(lex_symbol_range))
      {
        temporary_pattern_range =
        {
          lex_symbol_bucket_offsets_[std::get<0>(lex_symbol_range)],
          lex_symbol_bucket_offsets_[std::get<1>(lex_symbol_range) + 1] - 1
        };
      }
      if ((rlast != rend) && IsNotEmptyRange(temporary_pattern_range))
      {
        if (std::distance(rend, rlast) >= Index::kMaxFactorSize)
        {
          rfirst = rlast;
          rlast = std::next(rend, (std::distance(rend, rbegin) - k) % Index::kMaxFactorSize);
          BackwardSearchPatternInfixFactors(rfirst, rlast, pattern_range);
        }
        if (IsNotEmptyRange(temporary_pattern_range))
        {
          auto colex_symbol_range {CalculateSymbolRange(rlast, rend, -1)};
          if (IsNotEmptyRange(colex_symbol_range))
          {
            count += RangeCount(temporary_pattern_range, colex_symbol_range);
          }
        }
      }
      count += RangeSize(temporary_pattern_range);
    }
    pattern_range = {1, count};
  }
  else
  {
    pattern_range = {1, sub_factor_trie_.Count(std::next(rend), std::next(rbegin))};
  }
  return;
}

template <uint8_t max_factor_size>
template <typename Iterator, typename Range>
void Index<max_factor_size>::BackwardSearchPatternInfix
(
  Iterator rbegin,
  Iterator rend,
  Range &pattern_range,
  Range &pattern_range_ls
)
{
  auto rfirst {rbegin};
  auto rlast {std::prev(rbegin, std::distance(rend, rfirst) % Index::kMaxFactorSize)};
  if (rfirst == rlast)
  {
    rlast = std::prev(rbegin, Index::kMaxFactorSize);
  }
  while ((rfirst != rend) && (IsNotEmptyRange(pattern_range) || IsNotEmptyRange(pattern_range_ls)))
  {
    auto lex_symbol {lex_symbol_table_[SymbolsToInteger(std::next(rlast), std::next(rfirst))]};
    if (lex_symbol != 0)
    {
      auto begin_offset {lex_symbol_bucket_offsets_[lex_symbol]};
      if (IsNotEmptyRange(pattern_range))
      {
        pattern_range =
        {
          begin_offset + lex_bwt_.rank(std::get<0>(pattern_range), lex_symbol),
          begin_offset + lex_bwt_.rank(std::get<1>(pattern_range) + 1, lex_symbol) - 1
        };
      }
      if (IsNotEmptyRange(pattern_range_ls))
      {
        pattern_range_ls =
        {
          begin_offset + lex_bwt_.rank(std::get<0>(pattern_range_ls), lex_symbol),
          begin_offset + lex_bwt_.rank(std::get<1>(pattern_range_ls) + 1, lex_symbol) - 1
        };
      }
    }
    else
    {
      pattern_range = pattern_range_ls = {1, 0};
    }
    rfirst = rlast;
    rlast = std::prev(rlast, Index::kMaxFactorSize);
  }
  return;
}

template <uint8_t max_factor_size>
template <typename Iterator, typename Range>
void Index<max_factor_size>::BackwardSearchPatternPrefix
(
  Iterator rbegin,
  Iterator rend,
  Range &pattern_range,
  Range &pattern_range_ls
)
{
  uint64_t count {};
  uint64_t count_ls {};
  auto sl_factor_size {std::distance(rend, rbegin)};
  for (int64_t k {1}; k <= Index::kMaxFactorSize; ++k)
  {
    auto temporary_pattern_range {pattern_range};
    auto temporary_pattern_range_ls {temporary_pattern_range_ls};
    if (k != sl_factor_size)
    {
      auto rfirst {rbegin};
      auto rlast {std::prev(rbegin, k)};
      auto lex_symbol {lex_symbol_table_[SymbolsToInteger(std::next(rlast), std::next(rfirst))]};
      if (lex_symbol != 0)
      {
        auto begin_offset {lex_symbol_bucket_offsets_[lex_symbol]};
        if (IsNotEmptyRange(temporary_pattern_range))
        {
          temporary_pattern_range =
          {
            begin_offset + lex_bwt_.rank(std::get<0>(temporary_pattern_range), lex_symbol),
            begin_offset + lex_bwt_.rank(std::get<1>(temporary_pattern_range) + 1, lex_symbol) - 1
          };
        }
        if (IsNotEmptyRange(temporary_pattern_range_ls))
        {
          temporary_pattern_range_ls =
          {
            begin_offset + lex_bwt_.rank(std::get<0>(temporary_pattern_range_ls), lex_symbol),
            begin_offset + lex_bwt_.rank(std::get<1>(temporary_pattern_range_ls) + 1, lex_symbol) - 1
          };
        }
      }
      else
      {
        temporary_pattern_range = temporary_pattern_range_ls = {1, 0};
      }
      if (IsNotEmptyRange(temporary_pattern_range) || IsNotEmptyRange(temporary_pattern_range_ls))
      {
        if (std::distance(rend, rlast) > Index::kMaxFactorSize)
        {
          rfirst = rlast;
          rlast = std::next(rend, (sl_factor_size - k) % Index::kMaxFactorSize);
          BackwardSearchPatternInfix(rfirst, rlast, temporary_pattern_range, temporary_pattern_range_ls);
        }
        if (IsNotEmptyRange(temporary_pattern_range) || IsNotEmptyRange(temporary_pattern_range_ls))
        {
          auto colex_symbol_range {CalculateSymbolRange(rlast, rend, -1)};
          if (IsNotEmptyRange(colex_symbol_range))
          {
            count += RangeCount(temporary_pattern_range, colex_symbol_range);
            count_ls += RangeCount(temporary_pattern_range_ls, colex_symbol_range);
          }
        }
      }
    }
    else
    {
      auto colex_symbol_range {CalculateSymbolRange(rbegin, rend, -1)};
      if (IsNotEmptyRange(colex_symbol_range))
      {
        count += RangeCount(temporary_pattern_range, colex_symbol_range);
        count_ls += RangeCount(temporary_pattern_range_ls, colex_symbol_range);
      }
      break;
    }
  }
  pattern_range = {1, count};
  pattern_range_ls = {1, count_ls};
  return;
}
}
