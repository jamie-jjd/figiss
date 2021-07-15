#pragma once

#include <deque>
#include <map>
#include <memory>
#include <set>

#include <sdsl/wavelet_trees.hpp>

#include "sub_factor_trie.h"
#include "symbol_table.h"
#include "utility.h"

namespace project
{
class GrammarSymbolTable
{
public:

  GrammarSymbolTable () = default;
  GrammarSymbolTable (std::vector<uint64_t> &sorted_factor_integers);
  GrammarSymbolTable& operator= (GrammarSymbolTable &&);

  uint64_t Serialize
  (
    std::ostream &out,
    std::shared_ptr<SpaceNode> parent = nullptr,
    std::string const name = ""
  );
  void Load (std::istream &in);

  uint64_t operator[] (uint64_t const factor_integer) const;

  std::pair<uint64_t, uint64_t> SymbolRange
  (
    uint64_t const smallest_factor_integer,
    uint64_t const largest_factor_integer
  ) const;

  friend std::ostream& operator<< (std::ostream &out, GrammarSymbolTable const &symbol_table);

private:

  sdsl::sd_vector<> factor_integer_bits_;
  sdsl::sd_vector<>::rank_1_type factor_integer_rank_1_;
  sdsl::sd_vector<>::select_1_type factor_integer_select_1_;

};

GrammarSymbolTable::GrammarSymbolTable (std::vector<uint64_t> &sorted_factor_integers)
{
  factor_integer_bits_ = decltype(factor_integer_bits_)(std::begin(sorted_factor_integers), std::end(sorted_factor_integers));
  factor_integer_rank_1_.set_vector(&factor_integer_bits_);
}

GrammarSymbolTable& GrammarSymbolTable::operator= (GrammarSymbolTable &&grammar_symbol_table)
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

uint64_t GrammarSymbolTable::operator[] (uint64_t const factor_integer) const
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
) const
{
  std::pair<uint64_t, uint64_t> symbol_range {1, 0};
  if (smallest_factor_integer < std::size(factor_integer_bits_))
  {
    std::get<0>(symbol_range) = factor_integer_rank_1_(smallest_factor_integer);
    if (largest_factor_integer < std::size(factor_integer_bits_))
    {
      std::get<1>(symbol_range) =
      (
        factor_integer_rank_1_(largest_factor_integer)
        - (factor_integer_bits_[largest_factor_integer] == 0)
      );
    }
    else
    {
      std::get<1>(symbol_range) = factor_integer_rank_1_(std::size(factor_integer_bits_)) - 1;
    }
  }
  return symbol_range;
}

std::ostream& operator<< (std::ostream &out, GrammarSymbolTable const &grammar_symbol_table)
{
  out << std::size(grammar_symbol_table.factor_integer_bits_) << "\n";
  out << "space:\n";
  out << "factor_integer_bits_: " << ProperSizeRepresentation(sdsl::size_in_bytes(grammar_symbol_table.factor_integer_bits_)) << "B\n";
  return out;
}

class SymbolBucketOffsets
{
public:

  SymbolBucketOffsets () = default;
  SymbolBucketOffsets (sdsl::int_vector<> const &offsets);
  SymbolBucketOffsets& operator= (SymbolBucketOffsets &&);

  uint64_t Serialize
  (
    std::ostream &out,
    std::shared_ptr<SpaceNode> parent = nullptr,
    std::string const name = ""
  );
  void Load (std::istream &in);

  uint64_t operator[] (uint64_t const symbol) const;

  friend std::ostream& operator<< (std::ostream &out, SymbolBucketOffsets const &symbol_bucket_offsets);

private:

  sdsl::sd_vector<> offset_bits_;
  sdsl::sd_vector<>::select_1_type offset_select_1_;

};

SymbolBucketOffsets::SymbolBucketOffsets (sdsl::int_vector<> const &offsets)
{
  offset_bits_ = decltype(offset_bits_)(std::begin(offsets), std::end(offsets));
  offset_select_1_.set_vector(&offset_bits_);
}

SymbolBucketOffsets& SymbolBucketOffsets::operator= (SymbolBucketOffsets &&symbol_bucket_offsets)
{
  if (this != &symbol_bucket_offsets)
  {
    offset_bits_ = std::move(symbol_bucket_offsets.offset_bits_);
    offset_select_1_.set_vector(&offset_bits_);
  }
  return *this;
}

uint64_t SymbolBucketOffsets::operator[] (uint64_t const symbol) const
{
  if (symbol < std::size(offset_bits_))
  {
    return offset_select_1_(symbol);
  }
  return std::size(offset_bits_);
}

std::ostream& operator<< (std::ostream &out, SymbolBucketOffsets const &symbol_bucket_offsets)
{
  for (uint64_t i {1}, offset {}; offset != (std::size(symbol_bucket_offsets.offset_bits_) - 1); ++i)
  {
    offset = symbol_bucket_offsets.offset_select_1_(i);
    out << offset << ((offset != (std::size(symbol_bucket_offsets.offset_bits_) - 1)) ? " " : "\n");
  }
  out << "space:\n";
  out << ProperSizeRepresentation(sdsl::size_in_bytes(symbol_bucket_offsets.offset_bits_)) << "B\n";
  return out;
}

uint64_t SymbolBucketOffsets::Serialize
(
  std::ostream &out,
  std::shared_ptr<SpaceNode> parent,
  std::string const name
)
{
  uint64_t size_in_bytes {};
  if (!parent)
  {
    sdsl::serialize(offset_bits_, out);
    sdsl::serialize(offset_select_1_, out);
  }
  else
  {
    auto node {std::make_shared<SpaceNode>(name)};
    node->AddLeaf("offset_bits_", sdsl::serialize(offset_bits_, out));
    node->AddLeaf("offset_select_1_", sdsl::serialize(offset_select_1_, out));
    parent->AddChild(node);
    size_in_bytes = node->GetSizeInBytes();
  }
  return size_in_bytes;
}

void SymbolBucketOffsets::Load (std::istream &in)
{
  offset_bits_.load(in);
  offset_select_1_.load(in);
  offset_select_1_.set_vector(&offset_bits_);
  return;
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
  SymbolBucketOffsets lex_symbol_bucket_offsets_;
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
  void InsertSubFactorsInSlFactor (Iterator begin, Iterator end, Trie &trie);
  template <typename Iterator>
  uint64_t SymbolsToInteger (Iterator it, Iterator end, int8_t const step = 1);
  std::deque<uint64_t> IntegerToSymbols (uint64_t integer);
  void CalculateLexTextSizeAndLexFactorIntegersAndTempColexToLex
  (
    sdsl::int_vector<> const &text,
    uint64_t &lex_text_size,
    std::set<uint64_t> &lex_factor_integers,
    std::map<uint64_t, uint64_t> &temp_colex_to_lex
  );
  template <typename Iterator>
  void CalculateLexTextSizeAndLexFactorIntegersAndTempColexToLexInSlFactor
  (
    Iterator rbegin,
    Iterator rend,
    uint64_t &lex_text_size,
    std::set<uint64_t> &lex_factor_integers,
    std::map<uint64_t, uint64_t> &temp_colex_to_lex
  );
  void CalculateLexSymbolTable (std::set<uint64_t> const &lex_factor_integers);
  void CalculateLexText (sdsl::int_vector<> const &text, sdsl::int_vector<> &lex_text);
  template <typename TextIterator, typename LexTextIterator>
  LexTextIterator CalculateLexTextInSlFactor
  (
    TextIterator rbegin,
    TextIterator rend,
    LexTextIterator lex_text_it
  );
  void CalculateColexSymbolTable (std::map<uint64_t, uint64_t> const &temp_colex_to_lex);
  void CalculateColexToLex (std::map<uint64_t, uint64_t> const &temp_colex_to_lex);
  void CalculateLexSymbolBucketOffsets
  (
    uint64_t const lex_text_alphabet_size,
    sdsl::int_vector<> const &lex_text
  );
  void CalculateLexBwt (sdsl::int_vector<> const &lex_text);

  template <typename Range> // [,]
  inline bool IsNotEmptyRange (Range const &range) const
  {
    return (std::get<0>(range) <= std::get<1>(range));
  }

  template <typename Range> // [,]
  inline uint64_t RangeSize (Range const &range) const
  {
    if (std::get<0>(range) <= std::get<1>(range))
    {
      return (std::get<1>(range) - std::get<0>(range) + 1);
    }
    return 0;
  }

  template <typename Iterator>
  bool CalculatePattern (Iterator begin, Iterator end, sdsl::int_vector<> &pattern);
  template <typename Iterator, typename Range>
  Iterator BackwardSearchSuffix
  (
    Iterator rbegin,
    Iterator rend,
    Range &pattern_range_s,
    Range &pattern_range_l
  );
  template <typename Iterator, typename Range>
  void BackwardSearchInfix
  (
    Iterator rbegin,
    Iterator rend,
    Range &pattern_range_s,
    Range &pattern_range_l
  );
  template <typename Iterator, typename Range>
  void BackwardSearchPrefix
  (
    Iterator rbegin,
    Iterator rend,
    Range &pattern_range_s,
    Range &pattern_range_l
  );
  template <typename Iterator>
  Iterator ClassifySuffix (Iterator rbegin, Iterator rend);
  template <typename Iterator>
  void CalculateSlFactor (Iterator const rend, Iterator &rfirst, Iterator &rlast);
  template <typename Iterator, typename Range>
  void BackwardSearchSuffixSlFactor (Iterator rbegin, Iterator rend, Range &pattern_range);
  template <typename Iterator, typename Range>
  void BackwardSearchInfixFactors (Iterator rbegin, Iterator rend, Range &pattern_range);
  template <typename Iterator, typename Range>
  void BackwardSearchPrefixSlFactor (Iterator rbegin, Iterator rend, Range &pattern_range);
  template <typename Iterator, typename Range>
  void CountWithinSlFactor
  (
    Iterator rbegin,
    Iterator rend,
    uint64_t const duplicate_factor_size,
    Range &pattern_range
  );
  template <typename Iterator>
  std::pair<uint64_t, uint64_t> CalculateSymbolRange
  (
    Iterator begin,
    Iterator end,
    int8_t const step = 1
  );
  template <typename Range>
  uint64_t RangeCount
  (
    Range const pattern_range,
    Range const colex_symbol_range
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
  std::map<uint64_t, uint64_t> temp_colex_to_lex;
  {
    CalculateLexTextSizeAndLexFactorIntegersAndTempColexToLex
    (
      text,
      lex_text_size,
      lex_factor_integers,
      temp_colex_to_lex
    );
    // std::cout << "lex_text_size:\n" << lex_text_size << "\n";
    // std::cout << "|lex_factor_integers|: " << std::size(lex_factor_integers) << "\n";
    // std::cout << "lex_factor_integers:\n";
    // uint64_t i {};
    // for (auto const &factor : lex_factor_integers)
    // {
    //   std::cout << i++ << ":";
    //   Print(IntegerToSymbols(factor), std::cout);
    // }
    // std::cout << "temp_colex_to_lex:\n";
    // for (auto const &pair : temp_colex_to_lex)
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
    CalculateColexSymbolTable(temp_colex_to_lex);
    // std::cout << colex_symbol_table_;
  }
  {
    CalculateColexToLex(temp_colex_to_lex);
    // Print(colex_to_lex_, std::cout);
  }
  {
    CalculateLexSymbolBucketOffsets(lex_text_alphabet_size, lex_text);
    // std::cout << lex_symbol_bucket_offsets_;
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
    lex_symbol_bucket_offsets_.Serialize(index_file);
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
    root->AccumalateSizeInBytes(lex_symbol_bucket_offsets_.Serialize(index_file, root, "lex_symbol_bucket_offsets_"));
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
  sub_factor_trie_.Load(index_file);
  lex_symbol_table_.Load(index_file);
  colex_symbol_table_.Load(index_file);
  colex_to_lex_.load(index_file);
  lex_symbol_bucket_offsets_.Load(index_file);
  lex_bwt_.load(index_file);
  return;
}

template <uint8_t max_factor_size>
template <typename Iterator>
uint64_t Index<max_factor_size>::Count (Iterator begin, Iterator end)
{
  uint64_t count {};
  sdsl::int_vector<> pattern;
  if (CalculatePattern(begin, end, pattern))
  {
    // Print(pattern, std::cout);
    std::pair<uint64_t, uint64_t> pattern_range_s {1, 0};
    std::pair<uint64_t, uint64_t> pattern_range_l {1, 0};
    auto rbegin {std::prev(std::end(pattern))};
    auto rend {std::prev(std::begin(pattern))};
    auto rfirst {rbegin};
    auto rlast {rbegin};
    rlast = BackwardSearchSuffix(rbegin, rend, pattern_range_s, pattern_range_l);
    if (rlast != rend)
    {
      while (IsNotEmptyRange(pattern_range_s) || IsNotEmptyRange(pattern_range_l))
      {
        CalculateSlFactor(rend, rfirst, rlast);
        if (rlast != rend)
        {
          BackwardSearchInfix(rfirst, rlast, pattern_range_s, pattern_range_l);
        }
        else
        {
          BackwardSearchPrefix(rfirst, rlast, pattern_range_s, pattern_range_l);
          break;
        }
      }
    }
    {
      count += RangeSize(pattern_range_s);
      count += RangeSize(pattern_range_l);
      auto size {std::distance(begin, end)};
      if (size < (Index::kMaxFactorSize - 1))
      {
        count += sub_factor_trie_.Count(std::next(rend), std::next(rbegin));
      }
      if (size < Index::kMaxFactorSize)
      {
        auto colex_symbol_range {CalculateSymbolRange(rbegin, rend, -1)};
        std::get<0>(colex_symbol_range) += (colex_symbol_table_[SymbolsToInteger(rbegin, rend, -1)] != 0);
        if (IsNotEmptyRange(colex_symbol_range))
        {
          count += RangeCount({0, std::size(lex_bwt_) - 1}, colex_symbol_range);
        }
      }
    }
  }
  return count;
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
      return symbol_table_[byte];
    }
  );
  // Print(text, std::cout);
}

template <uint8_t max_factor_size>
void Index<max_factor_size>::CalculateSubFactorTrie (sdsl::int_vector<> const &text, Trie &trie)
{
  auto text_rend {std::prev(std::begin(text))};
  auto text_it {std::prev(std::end(text), 3)};
  auto next_symbol {*std::prev(std::end(text), 2)};
  uint8_t sl_type {};
  uint8_t next_sl_type {Index::kL};
  auto sl_factor_end {std::prev(std::end(text))};
  while (text_it != text_rend)
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
      auto sl_factor_begin {std::next(text_it)};
      InsertSubFactorsInSlFactor(sl_factor_begin, sl_factor_end, trie);
      sl_factor_end = sl_factor_begin;
    }
    next_symbol = *text_it--;
    next_sl_type = sl_type;
  }
  InsertSubFactorsInSlFactor(std::begin(text), sl_factor_end, trie);
  return;
}

template <uint8_t max_factor_size>
template <typename Iterator>
void Index<max_factor_size>::InsertSubFactorsInSlFactor (Iterator begin, Iterator end, Trie &trie)
{
  auto rend {std::prev(begin)};
  auto rfirst {std::prev(end)};
  auto rlast {std::prev(rfirst, std::distance(begin, end) % Index::kMaxFactorSize)};
  if (rfirst == rlast)
  {
    rlast = std::prev(rfirst, Index::kMaxFactorSize);
  }
  while (rfirst != rend)
  {
    if (std::distance(rlast, rfirst) > 2)
    {
      for (auto it {std::next(rlast, 2)}; it != rfirst; ++it)
      {
        trie.Insert(it, rfirst);
      }
    }
    rfirst = rlast;
    rlast = std::prev(rlast, Index::kMaxFactorSize);
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
void Index<max_factor_size>::CalculateLexTextSizeAndLexFactorIntegersAndTempColexToLex
(
  sdsl::int_vector<> const &text,
  uint64_t &lex_text_size,
  std::set<uint64_t> &lex_factor_integers,
  std::map<uint64_t, uint64_t> &temp_colex_to_lex
)
{
  {
    ++lex_text_size;
    lex_factor_integers.insert(0);
    temp_colex_to_lex[0] = 0;
  }
  auto text_rend {std::prev(std::begin(text))};
  auto text_it {std::prev(std::end(text), 3)};
  auto next_symbol {*std::prev(std::end(text), 2)};
  uint8_t sl_type {};
  uint8_t next_sl_type {Index::kL};
  auto sl_factor_end {std::prev(std::end(text))};
  while (text_it != text_rend)
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
      auto sl_factor_begin {std::next(text_it)};
      CalculateLexTextSizeAndLexFactorIntegersAndTempColexToLexInSlFactor
      (
        std::prev(sl_factor_end),
        std::prev(sl_factor_begin),
        lex_text_size,
        lex_factor_integers,
        temp_colex_to_lex
      );
      sl_factor_end = sl_factor_begin;
    }
    next_symbol = *text_it--;
    next_sl_type = sl_type;
  }
  CalculateLexTextSizeAndLexFactorIntegersAndTempColexToLexInSlFactor
  (
    std::prev(sl_factor_end),
    text_rend,
    lex_text_size,
    lex_factor_integers,
    temp_colex_to_lex
  );
  return;
}

template <uint8_t max_factor_size>
template <typename Iterator>
void Index<max_factor_size>::CalculateLexTextSizeAndLexFactorIntegersAndTempColexToLexInSlFactor
(
  Iterator rbegin,
  Iterator rend,
  uint64_t &lex_text_size,
  std::set<uint64_t> &lex_factor_integers,
  std::map<uint64_t, uint64_t> &temp_colex_to_lex
)
{
  auto rfirst {rbegin};
  auto rlast {std::prev(rbegin, std::distance(rend, rbegin) % Index::kMaxFactorSize)};
  if (rfirst == rlast)
  {
    rlast = std::prev(rbegin, Index::kMaxFactorSize);
  }
  while (rfirst != rend)
  {
    ++lex_text_size;
    auto lex_factor_integer {SymbolsToInteger(std::next(rlast), std::next(rfirst))};
    auto colex_factor_integer {SymbolsToInteger(rfirst, rlast, -1)};
    lex_factor_integers.insert(lex_factor_integer);
    temp_colex_to_lex[colex_factor_integer] = lex_factor_integer;
    rfirst = rlast;
    rlast = std::prev(rlast, Index::kMaxFactorSize);
  }
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
  auto text_rend {std::prev(std::begin(text))};
  auto text_it {std::prev(std::end(text), 3)};
  auto next_symbol {*std::prev(std::end(text), 2)};
  uint8_t sl_type {};
  uint8_t next_sl_type {Index::kL};
  auto sl_factor_end {std::prev(std::end(text))};
  auto lex_text_it {std::prev(std::end(lex_text), 2)};
  while (text_it != text_rend)
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
      auto sl_factor_begin {std::next(text_it)};
      lex_text_it = CalculateLexTextInSlFactor(std::prev(sl_factor_end), std::prev(sl_factor_begin), lex_text_it);
      sl_factor_end = sl_factor_begin;
    }
    next_symbol = *text_it--;
    next_sl_type = sl_type;
  }
  CalculateLexTextInSlFactor(std::prev(sl_factor_end), text_rend, lex_text_it);
  return;
}

template <uint8_t max_factor_size>
template <typename TextIterator, typename LexTextIterator>
LexTextIterator Index<max_factor_size>::CalculateLexTextInSlFactor
(
  TextIterator rbegin,
  TextIterator rend,
  LexTextIterator lex_text_it
)
{
  auto rfirst {rbegin};
  auto rlast {std::prev(rbegin, std::distance(rend, rbegin) % Index::kMaxFactorSize)};
  if (rfirst == rlast)
  {
    rlast = std::prev(rbegin, Index::kMaxFactorSize);
  }
  while (rfirst != rend)
  {
    *lex_text_it-- = lex_symbol_table_[SymbolsToInteger(std::next(rlast), std::next(rfirst))];
    rfirst = rlast;
    rlast = std::prev(rlast, Index::kMaxFactorSize);
  }
  return lex_text_it;
}

template <uint8_t max_factor_size>
void Index<max_factor_size>::CalculateColexSymbolTable (std::map<uint64_t, uint64_t> const &temp_colex_to_lex)
{
  std::vector<uint64_t> sorted_colex_factor_integers(std::size(temp_colex_to_lex));
  auto it {std::begin(sorted_colex_factor_integers)};
  for (auto const &pair : temp_colex_to_lex)
  {
    *it++ = std::get<0>(pair);
  }
  colex_symbol_table_ = decltype(colex_symbol_table_)(sorted_colex_factor_integers);
  return;
}

template <uint8_t max_factor_size>
void Index<max_factor_size>::CalculateColexToLex (std::map<uint64_t, uint64_t> const &temp_colex_to_lex)
{
  colex_to_lex_.width(sdsl::bits::hi(std::size(temp_colex_to_lex) - 1) + 1);
  colex_to_lex_.resize(std::size(temp_colex_to_lex));
  auto it {std::begin(colex_to_lex_)};
  for (auto const &pair : temp_colex_to_lex)
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
  sdsl::int_vector<> offsets
  (
    lex_text_alphabet_size + 1,
    0,
    sdsl::bits::hi(std::size(lex_text)) + 1
  );
  for (auto const lex_symbol : lex_text)
  {
    ++offsets[lex_symbol];
  }
  std::partial_sum(std::begin(offsets), std::end(offsets), std::begin(offsets));
  lex_symbol_bucket_offsets_ = decltype(lex_symbol_bucket_offsets_)(offsets);
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
template <typename Iterator>
bool Index<max_factor_size>::CalculatePattern (Iterator begin, Iterator end, sdsl::int_vector<> &pattern)
{
  pattern.width(symbol_table_.GetEffectiveAlphabetWidth());
  pattern.resize(std::distance(begin, end));
  auto pattern_it {std::begin(pattern)};
  for (auto it {begin}; it != end; ++it, ++pattern_it)
  {
    *pattern_it = symbol_table_[*it];
    if (*pattern_it == 0)
    {
      return false;
    }
  }
  return true;
}

template <uint8_t max_factor_size>
template <typename Iterator, typename Range>
void Index<max_factor_size>::CountWithinSlFactor
(
  Iterator rbegin,
  Iterator rend,
  uint64_t const duplicate_factor_size,
  Range &pattern_range
)
{
  uint64_t count {};
  auto size {static_cast<uint64_t>(std::distance(rend, rbegin))};
  for (uint64_t k {1}; (k <= Index::kMaxFactorSize) && (k <= size); ++k)
  {
    if (k != duplicate_factor_size)
    {
      auto rfirst {rbegin};
      auto rlast {std::prev(rbegin, k)};
      Range temp_pattern_range {1, 0};
      {
        auto lex_symbol_range {CalculateSymbolRange(std::next(rlast), std::next(rfirst))};
        if (IsNotEmptyRange(lex_symbol_range))
        {
          temp_pattern_range =
          {
            lex_symbol_bucket_offsets_[std::get<0>(lex_symbol_range)],
            lex_symbol_bucket_offsets_[std::get<1>(lex_symbol_range) + 1] - 1
          };
        }
      }
      auto is_not_counted {true};
      if ((rlast != rend) && IsNotEmptyRange(temp_pattern_range))
      {
        auto prefix_infix_size {size - k};
        if (prefix_infix_size >= Index::kMaxFactorSize)
        {
          rfirst = rlast;
          rlast = std::next(rend, prefix_infix_size % Index::kMaxFactorSize);
          BackwardSearchInfixFactors(rfirst, rlast, temp_pattern_range);
        }
        if ((rlast != rend) && IsNotEmptyRange(temp_pattern_range))
        {
          auto colex_symbol_range {CalculateSymbolRange(rlast, rend, -1)};
          if (IsNotEmptyRange(colex_symbol_range))
          {
            count += RangeCount(temp_pattern_range, colex_symbol_range);
          }
          is_not_counted = false;
        }
      }
      if (is_not_counted)
      {
        count += RangeSize(temp_pattern_range);
      }
    }
  }
  pattern_range = {1, count};
  return;
}

template <uint8_t max_factor_size>
template <typename Iterator>
std::pair<uint64_t, uint64_t> Index<max_factor_size>::CalculateSymbolRange
(
  Iterator begin,
  Iterator end,
  int8_t const step
)
{
  uint64_t factor_integer {};
  auto it {begin};
  auto k {Index::kMaxFactorSize};
  while (it != std::prev(end, step))
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
void Index<max_factor_size>::BackwardSearchInfixFactors
(
  Iterator rbegin,
  Iterator rend,
  Range &pattern_range
)
{
  auto rfirst {rbegin};
  auto rlast {std::prev(rbegin, std::distance(rend, rbegin) % Index::kMaxFactorSize)};
  if (rfirst == rlast)
  {
    rlast = std::prev(rbegin, Index::kMaxFactorSize);
  }
  while (rfirst != rend)
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
    if (IsNotEmptyRange(pattern_range))
    {
      rfirst = rlast;
      rlast = std::prev(rlast, Index::kMaxFactorSize);
    }
    else
    {
      return;
    }
  }
  return;
}

template <uint8_t max_factor_size>
template <typename Iterator, typename Range>
void Index<max_factor_size>::BackwardSearchPrefixSlFactor
(
  Iterator rbegin,
  Iterator rend,
  Range &pattern_range
)
{
  uint64_t count {};
  auto size {static_cast<uint64_t>(std::distance(rend, rbegin))};
  for (uint64_t k {1}; (k <= Index::kMaxFactorSize) && (k < size); ++k)
  {
    auto temp_pattern_range {pattern_range};
    auto rfirst {rbegin};
    auto rlast {std::prev(rbegin, k)};
    auto lex_symbol {lex_symbol_table_[SymbolsToInteger(std::next(rlast), std::next(rfirst))]};
    if (lex_symbol != 0)
    {
      auto begin_offset {lex_symbol_bucket_offsets_[lex_symbol]};
      temp_pattern_range =
      {
        begin_offset + lex_bwt_.rank(std::get<0>(temp_pattern_range), lex_symbol),
        begin_offset + lex_bwt_.rank(std::get<1>(temp_pattern_range) + 1, lex_symbol) - 1
      };
    }
    else
    {
      temp_pattern_range = {1, 0};
    }
    auto is_not_counted {true};
    if (IsNotEmptyRange(temp_pattern_range))
    {
      auto prefix_infix_size {size - k};
      if (prefix_infix_size >= Index::kMaxFactorSize)
      {
        rfirst = rlast;
        rlast = std::next(rend, prefix_infix_size % Index::kMaxFactorSize);
        BackwardSearchInfixFactors(rfirst, rlast, temp_pattern_range);
      }
      if ((rlast != rend) && IsNotEmptyRange(temp_pattern_range))
      {
        auto colex_symbol_range {CalculateSymbolRange(rlast, rend, -1)};
        if (IsNotEmptyRange(colex_symbol_range))
        {
          count += RangeCount(temp_pattern_range, colex_symbol_range);
        }
        is_not_counted = false;
      }
    }
    if (is_not_counted)
    {
      count += RangeSize(temp_pattern_range);
    }
  }
  if (size <= Index::kMaxFactorSize)
  {
    auto temp_pattern_range {pattern_range};
    auto colex_symbol_range {CalculateSymbolRange(rbegin, rend, -1)};
    if (IsNotEmptyRange(colex_symbol_range))
    {
      count += RangeCount(temp_pattern_range, colex_symbol_range);
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
    colex_symbol <= std::get<1>(colex_symbol_range);
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
Iterator Index<max_factor_size>::BackwardSearchSuffix
(
  Iterator rbegin,
  Iterator rend,
  Range &pattern_range_s,
  Range &pattern_range_l
)
{
  auto rfirst {rbegin};
  auto rlast {ClassifySuffix(rbegin, rend)};
  if (rlast != rend)
  {
    CalculateSlFactor(rend, rfirst, rlast);
    auto sl_factor_size {static_cast<uint64_t>(std::distance(rlast, rbegin))};
    auto suffix_s_size {static_cast<uint64_t>(std::distance(rfirst, rbegin))};
    if (rfirst != rbegin)
    {
      BackwardSearchSuffixSlFactor(rbegin, rfirst, pattern_range_s);
      if (IsNotEmptyRange(pattern_range_s))
      {
        if (rlast != rend)
        {
          BackwardSearchInfixFactors(rfirst, rlast, pattern_range_s);
        }
        else
        {
          BackwardSearchPrefixSlFactor(rfirst, rlast, pattern_range_s);
        }
      }
      if (rlast != rend)
      {
        if ((sl_factor_size % Index::kMaxFactorSize) != (suffix_s_size % Index::kMaxFactorSize))
        {
          BackwardSearchSuffixSlFactor(rbegin, rlast, pattern_range_l);
        }
      }
      else
      {
        auto duplicate_factor_size {suffix_s_size % Index::kMaxFactorSize};
        if (duplicate_factor_size == 0)
        {
          duplicate_factor_size = Index::kMaxFactorSize;
        }
        CountWithinSlFactor(rbegin, rend, duplicate_factor_size, pattern_range_l);
      }
    }
    else
    {
      if (rlast != rend)
      {
        BackwardSearchSuffixSlFactor(rbegin, rlast, pattern_range_l);
      }
      else
      {
        CountWithinSlFactor(rbegin, rend, 0, pattern_range_l);
      }
    }
  }
  else
  {
    CountWithinSlFactor(rbegin, rend, 0, pattern_range_s);
  }
  return rlast;
}

template <uint8_t max_factor_size>
template <typename Iterator>
Iterator Index<max_factor_size>::ClassifySuffix (Iterator rbegin, Iterator rend)
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
  return std::prev(it);
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
  auto next_sl_type {Index::kL};
  rfirst = rlast--;
  while ((rlast != rend) && !((next_sl_type == Index::kS) && (*rlast > *std::next(rlast))))
  {
    if((next_sl_type == Index::kL) && (*rlast < *std::next(rlast)))
    {
      next_sl_type = Index::kS;
    }
    --rlast;
  }
  return;
}

template <uint8_t max_factor_size>
template <typename Iterator, typename Range>
void Index<max_factor_size>::BackwardSearchSuffixSlFactor
(
  Iterator rbegin,
  Iterator rend,
  Range &pattern_range
)
{
  auto rfirst {rbegin};
  auto rlast {std::prev(rbegin, std::distance(rend, rbegin) % Index::kMaxFactorSize)};
  if (rlast == rbegin)
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
      BackwardSearchInfixFactors(rlast, rend, pattern_range);
    }
  }
  return;
}

template <uint8_t max_factor_size>
template <typename Iterator, typename Range>
void Index<max_factor_size>::BackwardSearchInfix
(
  Iterator rbegin,
  Iterator rend,
  Range &pattern_range_s,
  Range &pattern_range_l
)
{
  auto rfirst {rbegin};
  auto rlast {std::prev(rbegin, std::distance(rend, rbegin) % Index::kMaxFactorSize)};
  if (rfirst == rlast)
  {
    rlast = std::prev(rbegin, Index::kMaxFactorSize);
  }
  while ((rfirst != rend) && (IsNotEmptyRange(pattern_range_s) || IsNotEmptyRange(pattern_range_l)))
  {
    auto lex_symbol {lex_symbol_table_[SymbolsToInteger(std::next(rlast), std::next(rfirst))]};
    if (lex_symbol != 0)
    {
      auto begin_offset {lex_symbol_bucket_offsets_[lex_symbol]};
      if (IsNotEmptyRange(pattern_range_s))
      {
        pattern_range_s =
        {
          begin_offset + lex_bwt_.rank(std::get<0>(pattern_range_s), lex_symbol),
          begin_offset + lex_bwt_.rank(std::get<1>(pattern_range_s) + 1, lex_symbol) - 1
        };
      }
      if (IsNotEmptyRange(pattern_range_l))
      {
        pattern_range_l =
        {
          begin_offset + lex_bwt_.rank(std::get<0>(pattern_range_l), lex_symbol),
          begin_offset + lex_bwt_.rank(std::get<1>(pattern_range_l) + 1, lex_symbol) - 1
        };
      }
    }
    else
    {
      pattern_range_s = pattern_range_l = {1, 0};
    }
    rfirst = rlast;
    rlast = std::prev(rlast, Index::kMaxFactorSize);
  }
  return;
}

template <uint8_t max_factor_size>
template <typename Iterator, typename Range>
void Index<max_factor_size>::BackwardSearchPrefix
(
  Iterator rbegin,
  Iterator rend,
  Range &pattern_range_s,
  Range &pattern_range_l
)
{
  uint64_t count_s {};
  uint64_t count_l {};
  auto size {static_cast<uint64_t>(std::distance(rend, rbegin))};
  for (uint64_t k {1}; (k <= Index::kMaxFactorSize) && (k < size); ++k)
  {
    auto temp_pattern_range_s {pattern_range_s};
    auto temp_pattern_range_l {pattern_range_l};
    auto rfirst {rbegin};
    auto rlast {std::prev(rbegin, k)};
    auto lex_symbol {lex_symbol_table_[SymbolsToInteger(std::next(rlast), std::next(rfirst))]};
    if (lex_symbol != 0)
    {
      auto begin_offset {lex_symbol_bucket_offsets_[lex_symbol]};
      if (IsNotEmptyRange(temp_pattern_range_s))
      {
        temp_pattern_range_s =
        {
          begin_offset + lex_bwt_.rank(std::get<0>(temp_pattern_range_s), lex_symbol),
          begin_offset + lex_bwt_.rank(std::get<1>(temp_pattern_range_s) + 1, lex_symbol) - 1
        };
      }
      if (IsNotEmptyRange(temp_pattern_range_l))
      {
        temp_pattern_range_l =
        {
          begin_offset + lex_bwt_.rank(std::get<0>(temp_pattern_range_l), lex_symbol),
          begin_offset + lex_bwt_.rank(std::get<1>(temp_pattern_range_l) + 1, lex_symbol) - 1
        };
      }
    }
    else
    {
      temp_pattern_range_s = temp_pattern_range_l = {1, 0};
    }
    auto is_not_counted {true};
    if (IsNotEmptyRange(temp_pattern_range_s) || IsNotEmptyRange(temp_pattern_range_l))
    {
      auto prefix_infix_size {size - k};
      if (prefix_infix_size >= Index::kMaxFactorSize)
      {
        rfirst = rlast;
        rlast = std::next(rend, prefix_infix_size % Index::kMaxFactorSize);
        BackwardSearchInfix(rfirst, rlast, temp_pattern_range_s, temp_pattern_range_l);
      }
      if ((rlast != rend) && (IsNotEmptyRange(temp_pattern_range_s) || IsNotEmptyRange(temp_pattern_range_l)))
      {
        auto colex_symbol_range {CalculateSymbolRange(rlast, rend, -1)};
        if (IsNotEmptyRange(colex_symbol_range))
        {
          if (IsNotEmptyRange(temp_pattern_range_s))
          {
            auto pattern_range {temp_pattern_range_s};
            count_s += RangeCount(pattern_range, colex_symbol_range);
          }
          if (IsNotEmptyRange(temp_pattern_range_l))
          {
            auto pattern_range {temp_pattern_range_l};
            count_l += RangeCount(pattern_range, colex_symbol_range);
          }
        }
        is_not_counted = false;
      }
      if (is_not_counted)
      {
        count_s += RangeSize(temp_pattern_range_s);
        count_l += RangeSize(temp_pattern_range_l);
      }
    }
  }
  if (size <= Index::kMaxFactorSize)
  {
    auto colex_symbol_range {CalculateSymbolRange(rbegin, rend, -1)};
    if (IsNotEmptyRange(colex_symbol_range))
    {
      if (IsNotEmptyRange(pattern_range_s))
      {
        auto pattern_range {pattern_range_s};
        count_s += RangeCount(pattern_range, colex_symbol_range);
      }
      if (IsNotEmptyRange(pattern_range_l))
      {
        auto pattern_range {pattern_range_l};
        count_l += RangeCount(pattern_range, colex_symbol_range);
      }
    }
  }
  pattern_range_s = {1, count_s};
  pattern_range_l = {1, count_l};
  return;
}
}
