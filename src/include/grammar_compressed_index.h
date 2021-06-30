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

class TinyPatternTrie
{
public:

  TinyPatternTrie () = default;
  TinyPatternTrie (Trie const &trie) noexcept;
  TinyPatternTrie (TinyPatternTrie const &) = default;
  TinyPatternTrie& operator= (TinyPatternTrie &&) noexcept;

  uint64_t Serialize
  (
    std::ostream &out,
    std::shared_ptr<SpaceNode> parent = nullptr,
    std::string const name = ""
  );
  void Load (std::istream &in);

  inline auto GetOffsetRange (uint64_t const level_order) const noexcept
  {
    return std::make_pair
    (
      level_order_select_1_(level_order + 1) - (level_order + 1) + 1,
      level_order_select_1_(level_order + 2) - (level_order + 2) + 1
    );
  }

  inline auto GetLabel (uint64_t const offset) const noexcept
  {
    return labels_[offset];
  }

  inline auto GetCount (uint64_t const offset) const noexcept
  {
    return counts_[offset];
  }

  friend std::ostream& operator<< (std::ostream &out, std::pair<TinyPatternTrie, bool> const &pair);

private:

  sdsl::bit_vector level_order_bits_;
  sdsl::bit_vector::select_1_type level_order_select_1_;
  sdsl::int_vector<> labels_;
  sdsl::int_vector<> counts_;

};

TinyPatternTrie::TinyPatternTrie (Trie const &trie) noexcept
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

TinyPatternTrie& TinyPatternTrie::operator= (TinyPatternTrie &&tiny_pattern_trie) noexcept
{
  if (this != &tiny_pattern_trie)
  {
    level_order_bits_ = std::move(tiny_pattern_trie.level_order_bits_);
    level_order_select_1_ = std::move(tiny_pattern_trie.level_order_select_1_);
    level_order_select_1_.set_vector(&level_order_bits_);
    labels_ = std::move(tiny_pattern_trie.labels_);
    counts_ = std::move(tiny_pattern_trie.counts_);
  }
  return *this;
}

uint64_t TinyPatternTrie::Serialize
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

void TinyPatternTrie::Load (std::istream &in)
{
  level_order_bits_.load(in);
  level_order_select_1_.load(in);
  level_order_select_1_.set_vector(&level_order_bits_);
  labels_.load(in);
  counts_.load(in);
  return;
}

std::ostream& operator<< (std::ostream &out, std::pair<TinyPatternTrie, bool> const &pair)
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
    auto offset_range {trie.GetOffsetRange(offset)};
    if (depth != 0)
    {
      out << depth << ":" << trie.GetLabel(offset) << "(" << trie.GetCount(offset) << ")\n";
    }
    else
    {
      std::get<1>(offset_range) = std::get<0>(offset_range);
      std::get<0>(offset_range) = 0;
    }
    if (is_level_order)
    {
      for (auto offset {std::get<0>(offset_range)}; offset != std::get<1>(offset_range); ++offset)
      {
        offsets.emplace_back(offset, depth + 1);
      }
    }
    else
    {
      for (auto offset {std::get<1>(offset_range) - 1}; offset != (std::get<0>(offset_range) - 1); --offset)
      {
        offsets.emplace_back(offset, depth + 1);
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

  friend std::ostream& operator<< (std::ostream &out, GrammarSymbolTable const &symbol_table);

private:

  sdsl::sd_vector<> factor_integer_bits_;
  sdsl::sd_vector<>::rank_1_type factor_integer_rank_1_;

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
  if ((factor_integer <= std::size(factor_integer_bits_)) && factor_integer_bits_[factor_integer])
  {
    return factor_integer_rank_1_(factor_integer);
  }
  return 0;
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

  template <typename PatternIterator>
  uint64_t Count (PatternIterator begin, PatternIterator end);

  friend std::ostream& operator<< <max_factor_size> (std::ostream &out, Index<max_factor_size> const &index);

private:

  SymbolTable symbol_table_;
  TinyPatternTrie tiny_pattern_trie_;
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

  void CalculateTinyPatternTrie (sdsl::int_vector<> const &text, Trie &trie);
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

};

template <uint8_t max_factor_size>
Index<max_factor_size>::Index (std::filesystem::path const &byte_text_path)
{
  std::cout << "construct index of " << std::filesystem::canonical(byte_text_path) << "\n";
  sdsl::int_vector<> text;
  CalculateSymbolTableAndText(byte_text_path, text);
  {
    Trie trie;
    CalculateTinyPatternTrie(text, trie);
    // std::cout << std::make_pair(trie, /*is_level_order=*/true);
    tiny_pattern_trie_ = decltype(tiny_pattern_trie_)(trie);
    // std::cout << std::make_pair(tiny_pattern_trie_, /*is_level_order=*/true);
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
    tiny_pattern_trie_.Serialize(index_file);
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
    root->AccumalateSizeInBytes(tiny_pattern_trie_.Serialize(index_file, root, "tiny_pattern_trie_"));
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
  tiny_pattern_trie_.Load(index_file);
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
  out << std::make_pair(index.tiny_pattern_trie_, true);
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
void Index<max_factor_size>::CalculateTinyPatternTrie (sdsl::int_vector<> const &text, Trie &trie)
{
  auto text_it {std::begin(text)};
  auto text_last {std::prev(std::end(text))};
  while (std::distance(text_it, text_last) >= (Index::kMaxFactorSize - 1))
  {
    auto it {text_it};
    auto end {std::next(it, Index::kMaxFactorSize - 1)};
    while (it != end)
    {
      trie.Insert(it, end);
      ++it;
    }
    ++text_it;
  }
  while (text_it != text_last)
  {
    trie.Insert(text_it, text_last);
    ++text_it;
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

// template <typename Range>
// constexpr bool IsNotEmptyRange (Range const &range)
// {
//   return (std::get<0>(range) <= std::get<1>(range));
// }
//
// template <typename Range>
// constexpr uint64_t CalculateRangeSize (Range const &range)
// {
//   if (std::get<0>(range) <= std::get<1>(range))
//   {
//     return (std::get<1>(range) - std::get<0>(range) + 1);
//   }
//   return 0;
// }
//
// template <typename WaveletTree>
// uint64_t RangeCount
// (
//   WaveletTree const &wavelet_tree,
//   uint64_t const begin_offset,
//   uint64_t const end_offset,
//   uint64_t const begin_value,
//   uint64_t const end_value
// )
// {
//   uint64_t count {};
//   for (uint64_t value {begin_value}; value != end_value; ++value)
//   {
//     count += (wavelet_tree.rank(end_offset, value) - wavelet_tree.rank(begin_offset, value));
//   }
//   return count;
// }
//
// template <typename TextIterator>
// void CalculateSlFactor
// (
//   TextIterator const rend,
//   TextIterator &rfirst,
//   TextIterator &rlast
// )
// {
//   uint64_t prev_sl_type {L};
//   rfirst = rlast--;
//   while ((rlast != rend) && !((prev_sl_type == S) && (*rlast > *std::next(rlast))))
//   {
//     if((prev_sl_type == L) && (*rlast < *std::next(rlast)))
//     {
//       prev_sl_type = S;
//     }
//     --rlast;
//   }
//   return;
// }
//
// template <typename StaticGrammarTrie, typename SlFactorIterator>
// std::pair<uint64_t, uint64_t> LookUpSlFactor
// (
//   StaticGrammarTrie const &trie,
//   SlFactorIterator it,
//   SlFactorIterator last,
//   bool const is_exact = true, // false: prefix
//   bool const is_rank = true // false: count
// )
// {
//   auto labels_begin {std::get<0>(trie.labels_range)};
//   std::pair<uint64_t, uint64_t> pair {1, 0};
//   uint64_t begin_offset {};
//   uint64_t end_offset {trie.branch_bits_select(1)};
//   uint64_t offset {};
//   while (begin_offset != end_offset)
//   {
//     while (begin_offset != end_offset)
//     {
//       offset = begin_offset + (end_offset - begin_offset) / 2;
//       auto branch_character {*std::next(labels_begin, trie.edge_begin_offsets[offset])};
//       if (*it == branch_character)
//       {
//         break;
//       }
//       else if (*it < branch_character)
//       {
//         end_offset = offset;
//       }
//       else
//       {
//         begin_offset = offset + 1;
//       }
//     }
//     if (begin_offset != end_offset)
//     {
//       auto edge_it {std::next(labels_begin, trie.edge_begin_offsets[offset])};
//       auto edge_end {std::next(labels_begin, trie.edge_prev_end_offsets[offset] + trie.step)};
//       while ((it != last) && (edge_it != edge_end) && (*it == *edge_it))
//       {
//         it += trie.step;
//         edge_it += trie.step;
//       }
//       if (it != last)
//       {
//         if (edge_it != edge_end)
//         {
//           break;
//         }
//         else
//         {
//           begin_offset = trie.branch_bits_select(offset + 1) - (offset + 1) + 1;
//           end_offset = trie.branch_bits_select(offset + 2) - (offset + 2) + 1;
//         }
//       }
//       else
//       {
//         if (is_rank)
//         {
//           if (is_exact)
//           {
//             if (edge_it == edge_end)
//             {
//               begin_offset = trie.branch_bits_select(offset + 1) + 1;
//               auto bit {trie.branch_bits[begin_offset]};
//               if ((bit == 1) || ((bit == 0) && trie.leftmost_ranks[offset] != trie.leftmost_ranks[begin_offset - (offset + 1)]))
//               {
//                 std::get<1>(pair) = trie.leftmost_ranks[offset];
//               }
//             }
//           }
//           else
//           {
//             std::get<0>(pair) = trie.leftmost_ranks[offset];
//             std::get<1>(pair) = trie.rightmost_ranks[offset];
//           }
//         }
//         else
//         {
//           std::get<1>(pair) = trie.counts[offset];
//         }
//         break;
//       }
//     }
//   }
//   return pair;
// }
//
// template
// <
//   typename Index,
//   typename PatternRange,
//   typename PatternIterator
// >
// void BackwardSearchPatternPrefix
// (
//   Index const &index,
//   PatternRange &pattern_range_l,
//   PatternRange &pattern_range_s,
//   PatternIterator rfirst,
//   PatternIterator rlast
// )
// {
//   auto colex_rank_range {LookUpSlFactor(index.colex_grammar_rank_trie, rfirst, rlast, false)};
//   if (IsNotEmptyRange(colex_rank_range))
//   {
//     if (IsNotEmptyRange(pattern_range_l))
//     {
//       pattern_range_l =
//       {
//         1,
//         RangeCount
//         (
//           index.colex_bwt,
//           std::get<0>(pattern_range_l),
//           std::get<1>(pattern_range_l) + 1,
//           std::get<0>(colex_rank_range),
//           std::get<1>(colex_rank_range) + 1
//         )
//       };
//     }
//     if (IsNotEmptyRange(pattern_range_s))
//     {
//       pattern_range_s =
//       {
//         1,
//         RangeCount
//         (
//           index.colex_bwt,
//           std::get<0>(pattern_range_s),
//           std::get<1>(pattern_range_s) + 1,
//           std::get<0>(colex_rank_range),
//           std::get<1>(colex_rank_range) + 1
//         )
//       };
//     }
//   }
//   // {
//   //   if (IsNotEmptyRange(pattern_range_l))
//   //   {
//   //     std::cout << "pattern-l prefix:\n";
//   //   }
//   //   if (IsNotEmptyRange(pattern_range_s))
//   //   {
//   //     std::cout << "pattern-s prefix:\n";
//   //   }
//   //   Print(std::cout, rfirst, rlast, -1);
//   //   std::cout
//   //   << "->[" << std::get<0>(colex_rank_range) << "," << std::get<1>(colex_rank_range) << "]"
//   //   << "->L:[" << std::get<0>(pattern_range_l) << "," << std::get<1>(pattern_range_l) << "]"
//   //   << "->S:[" << std::get<0>(pattern_range_s) << "," << std::get<1>(pattern_range_s) << "]\n";
//   // }
//   return;
// }
//
// template
// <
//   typename Index,
//   typename PatternRange,
//   typename PatternIterator
// >
// auto BackwardSearchPatternInfix
// (
//   Index const &index,
//   PatternRange &pattern_range_l,
//   PatternRange &pattern_range_s,
//   PatternIterator rfirst,
//   PatternIterator rlast
// )
// {
//   auto colex_rank {std::get<1>(LookUpSlFactor(index.colex_grammar_rank_trie, rfirst, rlast))};
//   if (colex_rank != 0)
//   {
//     auto begin_offset {index.lex_rank_bucket_begin_offsets[index.colex_to_lex[colex_rank]]};
//     if (IsNotEmptyRange(pattern_range_l))
//     {
//       pattern_range_l =
//       {
//         begin_offset + index.colex_bwt.rank(std::get<0>(pattern_range_l), colex_rank),
//         begin_offset + index.colex_bwt.rank(std::get<1>(pattern_range_l) + 1, colex_rank) - 1
//       };
//     }
//     if (IsNotEmptyRange(pattern_range_s))
//     {
//       pattern_range_s =
//       {
//         begin_offset + index.colex_bwt.rank(std::get<0>(pattern_range_s), colex_rank),
//         begin_offset + index.colex_bwt.rank(std::get<1>(pattern_range_s) + 1, colex_rank) - 1
//       };
//     }
//   }
//   else
//   {
//     pattern_range_l = pattern_range_s = {1, 0};
//   }
//   // {
//   //   if (IsNotEmptyRange(pattern_range_l))
//   //   {
//   //     std::cout << "pattern-l sl-factor\n";
//   //   }
//   //   if (IsNotEmptyRange(pattern_range_s))
//   //   {
//   //     std::cout << "pattern-s sl-factor\n";
//   //   }
//   //   Print(std::cout, rfirst, rlast, -1);
//   //   std::cout
//   //   << "->[" << colex_rank << ":" << index.colex_to_lex[colex_rank] << "]"
//   //   << "->L:[" << std::get<0>(pattern_range_l) << "," << std::get<1>(pattern_range_l) << "]"
//   //   << "->S:[" << std::get<0>(pattern_range_s) << "," << std::get<1>(pattern_range_s) << "]\n";
//   // }
//   return rlast;
// }
//
// template <typename PatternIterator>
// auto CalculatePatternSuffixS (PatternIterator rbegin, PatternIterator rend)
// {
//   auto it {rbegin};
//   while ((std::prev(it) != rend) && (*std::prev(it) == *it))
//   {
//     --it;
//   }
//   if ((std::prev(it) != rend) && (*std::prev(it) < *it))
//   {
//     return rbegin;
//   }
//   return std::prev(it);
// }
//
// template
// <
//   typename Index,
//   typename PatternRange,
//   typename PatternIterator
// >
// auto BackwardSearchPatternSuffix
// (
//   Index const &index,
//   PatternRange &pattern_range_l,
//   PatternRange &pattern_range_s,
//   PatternIterator rbegin,
//   PatternIterator rend
// )
// {
//   auto rfirst {rbegin};
//   auto rlast {CalculatePatternSuffixS(rbegin, rend)};
//   if (rlast == rend)
//   {
//     std::get<1>(pattern_range_l) = std::get<1>
//     (
//       LookUpSlFactor
//       (
//         index.lex_grammar_count_trie,
//         std::next(rend),
//         std::next(rbegin),
//         false,
//         false
//       )
//     );
//     // {
//     //   std::cout << "pattern of single-run:\n";
//     //   Print(std::cout, std::next(rend), std::next(rbegin));
//     //   std::cout << "->L:[" << std::get<0>(pattern_range_l) << "," << std::get<1>(pattern_range_l) << "]\n";
//     // }
//   }
//   else
//   {
//     if (rlast != rbegin)
//     {
//       auto lex_rank_range {LookUpSlFactor(index.lex_grammar_rank_trie, std::next(rlast), std::next(rbegin), false)};
//       if (IsNotEmptyRange(lex_rank_range))
//       {
//         pattern_range_s =
//         {
//           index.lex_rank_bucket_begin_offsets[std::get<0>(lex_rank_range)],
//           index.lex_rank_bucket_begin_offsets[std::get<1>(lex_rank_range) + 1] - 1
//         };
//       }
//       // {
//       //   std::cout << "pattern-s suffix:\n";
//       //   Print(std::cout, std::next(rlast), std::next(rbegin));
//       //   std::cout << "->[" << std::get<0>(lex_rank_range) << "," << std::get<1>(lex_rank_range) << "]";
//       //   std::cout << "->S:[" << std::get<0>(pattern_range_s) << "," << std::get<1>(pattern_range_s) << "]\n";
//       // }
//       CalculateSlFactor(rend, rfirst, rlast);
//       if (IsNotEmptyRange(pattern_range_s))
//       {
//         if (rlast == rend)
//         {
//           auto colex_rank_range {LookUpSlFactor(index.colex_grammar_rank_trie, rfirst, rlast, false)};
//           pattern_range_s =
//           {
//             1,
//             RangeCount
//             (
//               index.colex_bwt,
//               std::get<0>(pattern_range_s),
//               std::get<1>(pattern_range_s) + 1,
//               std::get<0>(colex_rank_range),
//               std::get<1>(colex_rank_range) + 1
//             )
//           };
//           // {
//           //   std::cout << "pattern-s prefix:\n";
//           //   Print(std::cout, rfirst, rlast, -1);
//           //   std::cout
//           //   << "->[" << std::get<0>(colex_rank_range)
//           //   << ":" << index.colex_to_lex[std::get<0>(colex_rank_range)]
//           //   << "," << std::get<1>(colex_rank_range) << "]"
//           //   << "->S:[" << std::get<0>(pattern_range_s) << "," << std::get<1>(pattern_range_s) << "]\n";
//           // }
//         }
//         else
//         {
//           auto colex_rank {std::get<1>(LookUpSlFactor(index.colex_grammar_rank_trie, rfirst, rlast))};
//           if (colex_rank != 0)
//           {
//             auto begin_offset {index.lex_rank_bucket_begin_offsets[index.colex_to_lex[colex_rank]]};
//             pattern_range_s =
//             {
//               begin_offset + index.colex_bwt.rank(std::get<0>(pattern_range_s), colex_rank),
//               begin_offset + index.colex_bwt.rank(std::get<1>(pattern_range_s) + 1, colex_rank) - 1
//             };
//           }
//           else
//           {
//             pattern_range_s = {1, 0};
//           }
//           // {
//           //   std::cout << "pattern-s exact-sl-factor:\n";
//           //   Print(std::cout, rfirst, rlast, -1);
//           //   std::cout
//           //   << "->(" << colex_rank << ":" << index.colex_to_lex[colex_rank] << ")"
//           //   << "->S:[" << std::get<0>(pattern_range_s) << "," << std::get<1>(pattern_range_s) << "]\n";
//           // }
//         }
//       }
//     }
//     else
//     {
//       CalculateSlFactor(rend, rfirst, rlast);
//     }
//     if (rlast == rend)
//     {
//       std::get<1>(pattern_range_l) = std::get<1>
//       (
//         LookUpSlFactor
//         (
//           index.lex_grammar_count_trie,
//           std::next(rend),
//           std::next(rbegin),
//           false,
//           false
//         )
//       );
//       // {
//       //   std::cout << "pattern-l sub-sl-factor:\n";
//       //   Print(std::cout, std::next(rend), std::next(rbegin));
//       //   std::cout << "->L:[" << std::get<0>(pattern_range_l) << "," << std::get<1>(pattern_range_l) << "]\n";
//       // }
//     }
//     else
//     {
//       auto lex_rank_range {LookUpSlFactor(index.lex_grammar_rank_trie, std::next(rlast), std::next(rbegin), false)};
//       if (IsNotEmptyRange(lex_rank_range))
//       {
//         pattern_range_l =
//         {
//           index.lex_rank_bucket_begin_offsets[std::get<0>(lex_rank_range)],
//           (index.lex_rank_bucket_begin_offsets[std::get<1>(lex_rank_range) + 1] - 1)
//         };
//       }
//       // {
//       //   std::cout << "pattern-l suffix:\n";
//       //   Print(std::cout, std::next(rlast), std::next(rbegin));
//       //   std::cout << "->[" << std::get<0>(lex_rank_range) << "," << std::get<1>(lex_rank_range) << "]";
//       //   std::cout << "->L:[" << std::get<0>(pattern_range_l) << "," << std::get<1>(pattern_range_l) << "]\n";
//       // }
//     }
//   }
//   return rlast;
// }
//
// template <typename Index, typename PatternIterator>
// uint64_t Count
// (
//   Index const &index,
//   PatternIterator begin,
//   PatternIterator end
// )
// {
//   std::pair<uint64_t, uint64_t> pattern_range_l {1, 0};
//   std::pair<uint64_t, uint64_t> pattern_range_s {1, 0};
//   auto rbegin {std::prev(end)};
//   auto rend {std::prev(begin)};
//   auto rfirst {rbegin};
//   auto rlast {rbegin};
//   // Print(std::cout, begin, end);
//   rlast = BackwardSearchPatternSuffix(index, pattern_range_l, pattern_range_s, rbegin, rend);
//   if (rlast != rend)
//   {
//     while (IsNotEmptyRange(pattern_range_l) || IsNotEmptyRange(pattern_range_s))
//     {
//       CalculateSlFactor(rend, rfirst, rlast);
//       if (rlast != rend)
//       {
//         rlast = BackwardSearchPatternInfix(index, pattern_range_l, pattern_range_s, rfirst, rlast);
//       }
//       else
//       {
//         BackwardSearchPatternPrefix(index, pattern_range_l, pattern_range_s, rfirst, rlast);
//         break;
//       }
//     }
//   }
//   return (CalculateRangeSize(pattern_range_l) + CalculateRangeSize(pattern_range_s));
// }
}
