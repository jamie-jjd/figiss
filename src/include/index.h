#pragma once

#include <deque>
#include <map>
#include <memory>
#include <set>

#include <sdsl/wavelet_trees.hpp>

#include "grammar_symbol_table.h"
#include "sub_factor_trie.h"
#include "symbol_bucket_offsets.h"
#include "symbol_table.h"
#include "utility.h"

namespace figiss
{
class Index
{
public:

  static constexpr uint8_t kS {1};
  static constexpr uint8_t kL {0};

  Index () = default;
  Index (Index const&) = delete;
  Index (Index&&);
  Index
  (
    std::filesystem::path const &byte_text_path,
    uint8_t const max_factor_size_ = 4,
    bool const is_parameter_written = false
  );
  Index& operator= (Index const&) = delete;
  Index& operator= (Index&&);

  void Swap (Index&);

  uint64_t Serialize
  (
    std::filesystem::path const &index_path,
    std::shared_ptr<SpaceNode> root = nullptr
  );
  void Load (std::filesystem::path const &index_path);

  template <typename Iterator>
  uint64_t Count (Iterator begin, Iterator end);

  friend std::ostream& operator<<  (std::ostream &out, Index const &index);

private:

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

  uint8_t max_factor_size_;
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

};

Index::Index (Index&& index)
{
  if (this != &index)
  {
    this->Swap(index);
  }
}

void Index::Swap (Index& index)
{
  if (this != &index)
  {
    std::swap(max_factor_size_, index.max_factor_size_);
    symbol_table_.Swap(index.symbol_table_);
    sub_factor_trie_.Swap(index.sub_factor_trie_);
    lex_symbol_table_.Swap(index.lex_symbol_table_);
    colex_symbol_table_.Swap(index.colex_symbol_table_);
    colex_to_lex_.swap(index.colex_to_lex_);
    lex_symbol_bucket_offsets_.Swap(index.lex_symbol_bucket_offsets_);
    lex_bwt_.swap(index.lex_bwt_);
  }
  return;
}

Index::Index
(
  std::filesystem::path const &byte_text_path,
  uint8_t const max_factor_size,
  bool is_parameter_written
)
{
  if (max_factor_size == 0 || max_factor_size > 8)
  {
    std::cerr
    << "\033[31mfalied at constructing index of " << std::filesystem::canonical(byte_text_path)
    << "\nmax factor size must be integer within [1..8]\033[0m\n";
    return;
  }
  max_factor_size_ = max_factor_size;
  if (!is_parameter_written)
  {
    std::cout << "construct index of " << std::filesystem::canonical(byte_text_path).string() << "\n";
  }
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
  if (is_parameter_written)
  {
    auto parameter_path {std::filesystem::path{"../data/parameter/" + byte_text_path.filename().string() + ".figiss.parameter"}};
    if (!std::filesystem::exists(parameter_path.parent_path()))
    {
      std::filesystem::create_directories(parameter_path.parent_path());
    }
    std::ofstream out {parameter_path};
    std::cout << "write figiss parameters to " << std::filesystem::canonical(parameter_path).string() << "\n";
    out << "n," << std::size(text) << "\n";
    out << "n'," << std::size(lex_text) << "\n";
    {
      sdsl::int_vector<> buffer;
      sdsl::qsufsort::construct_sa(buffer, text);
      uint64_t r {};
      uint64_t prev_character {std::size(text)};
      for (auto it {std::begin(buffer)}; it != std::end(buffer); ++it)
      {
        if (*it != 0)
        {
          *it = text[*it - 1];
        }
        if (*it != prev_character)
        {
          prev_character = *it;
          ++r;
        }
      }
      out << "r," << r << "\n";
    }
    {
      uint64_t r {};
      uint64_t prev_lex_rank {std::size(lex_bwt_)};
      for (auto const& lex_rank : lex_bwt_)
      {
        if (prev_lex_rank != lex_rank)
        {
          prev_lex_rank = lex_rank;
          ++r;
        }
      }
      out << "r'," << r << "\n";
    }
    out << "sigma," << (*std::max_element(std::begin(text), std::end(text)) + 1) << "\n";
    out << "sigma'," << std::size(colex_to_lex_) << "\n";
  }
}

Index& Index::operator= (Index&& index)
{
  if (this != &index)
  {
    Index temp {std::move(index)};
    this->Swap(temp);
  }
  return *this;
}

uint64_t Index::Serialize
(
  std::filesystem::path const &index_path,
  std::shared_ptr<SpaceNode> root
)
{
  if
  (
    !index_path.parent_path().filename().string().empty()
    && !std::filesystem::exists(index_path.parent_path())
  )
  {
    std::filesystem::create_directories(index_path.parent_path());
  }
  std::fstream out {index_path, std::ios_base::out | std::ios_base::trunc};
  uint64_t size_in_bytes {};
  if (!root)
  {
    std::cout << "serialize index to " << std::filesystem::canonical(index_path).string() << "\n";
    sdsl::write_member(max_factor_size_, out);
    symbol_table_.Serialize(out);
    sub_factor_trie_.Serialize(out);
    lex_symbol_table_.Serialize(out);
    colex_symbol_table_.Serialize(out);;
    sdsl::serialize(colex_to_lex_, out);
    lex_symbol_bucket_offsets_.Serialize(out);
    sdsl::serialize(lex_bwt_, out);
  }
  else
  {
    root->AddLeaf("max_factor_size_", sdsl::write_member(max_factor_size_, out));
    root->AccumalateSizeInBytes(symbol_table_.Serialize(out, root, "symbol_table_"));
    root->AccumalateSizeInBytes(sub_factor_trie_.Serialize(out, root, "sub_factor_trie_"));
    root->AccumalateSizeInBytes(lex_symbol_table_.Serialize(out, root, "lex_symbol_table_"));
    root->AccumalateSizeInBytes(colex_symbol_table_.Serialize(out, root, "colex_symbol_table_"));
    root->AddLeaf("colex_to_lex_", sdsl::serialize(colex_to_lex_, out));
    root->AccumalateSizeInBytes(lex_symbol_bucket_offsets_.Serialize(out, root, "lex_symbol_bucket_offsets_"));
    root->AddLeaf("lex_bwt_", sdsl::serialize(lex_bwt_, out));
    size_in_bytes = root->GetSizeInBytes();
  }
  return size_in_bytes;
}

void Index::Load (std::filesystem::path const &index_path)
{
  std::ifstream in {index_path};
  std::cout << "load index from " << std::filesystem::canonical(index_path).string() << "\n";
  sdsl::read_member(max_factor_size_, in);
  symbol_table_.Load(in);
  sub_factor_trie_.Load(in);
  lex_symbol_table_.Load(in);
  colex_symbol_table_.Load(in);
  colex_to_lex_.load(in);
  lex_symbol_bucket_offsets_.Load(in);
  lex_bwt_.load(in);
  return;
}

template <typename Iterator>
uint64_t Index::Count (Iterator begin, Iterator end)
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
      if (size < (max_factor_size_ - 1))
      {
        count += sub_factor_trie_.Count(std::next(rend), std::next(rbegin));
      }
      if (size < max_factor_size_)
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

std::ostream& operator<< (std::ostream &out, Index const &index)
{
  out << index.symbol_table_;
  out << std::make_pair(index.sub_factor_trie_, true);
  out << index.lex_symbol_table_;
  out << index.colex_symbol_table_;
  out << index.colex_to_lex_ << "\n";
  out << index.lex_symbol_bucket_offsets_;
  out << index.lex_bwt_ << "\n";
  return out;
}

void Index::CalculateSymbolTableAndText
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
  for (uint64_t i {}; i != std::size(byte_text); ++i)
  {
    text[i] = symbol_table_[byte_text[i]];
  }
  // Print(text, std::cout);
}

void Index::CalculateSubFactorTrie (sdsl::int_vector<> const &text, Trie &trie)
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

template <typename Iterator>
void Index::InsertSubFactorsInSlFactor (Iterator begin, Iterator end, Trie &trie)
{
  auto rend {std::prev(begin)};
  auto rfirst {std::prev(end)};
  auto rlast {std::prev(rfirst, std::distance(begin, end) % max_factor_size_)};
  if (rfirst == rlast)
  {
    rlast = std::prev(rfirst, max_factor_size_);
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
    rlast = std::prev(rlast, max_factor_size_);
  }
  return;
}

template <typename Iterator>
uint64_t Index::SymbolsToInteger (Iterator it, Iterator end, int8_t const step)
{
  uint64_t result {};
  for (auto k {max_factor_size_}; (it != end) && (k != 0); it += step, --k)
  {
    result += *it * (1ULL << (symbol_table_.GetEffectiveAlphabetWidth() * (k - 1)));
  }
  return result;
}

std::deque<uint64_t> Index::IntegerToSymbols (uint64_t integer)
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

void Index::CalculateLexTextSizeAndLexFactorIntegersAndTempColexToLex
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

template <typename Iterator>
void Index::CalculateLexTextSizeAndLexFactorIntegersAndTempColexToLexInSlFactor
(
  Iterator rbegin,
  Iterator rend,
  uint64_t &lex_text_size,
  std::set<uint64_t> &lex_factor_integers,
  std::map<uint64_t, uint64_t> &temp_colex_to_lex
)
{
  auto rfirst {rbegin};
  auto rlast {std::prev(rbegin, std::distance(rend, rbegin) % max_factor_size_)};
  if (rfirst == rlast)
  {
    rlast = std::prev(rbegin, max_factor_size_);
  }
  while (rfirst != rend)
  {
    ++lex_text_size;
    auto lex_factor_integer {SymbolsToInteger(std::next(rlast), std::next(rfirst))};
    auto colex_factor_integer {SymbolsToInteger(rfirst, rlast, -1)};
    lex_factor_integers.insert(lex_factor_integer);
    temp_colex_to_lex[colex_factor_integer] = lex_factor_integer;
    rfirst = rlast;
    rlast = std::prev(rlast, max_factor_size_);
  }
  return;
}

void Index::CalculateLexSymbolTable (std::set<uint64_t> const &lex_factor_integers)
{
  std::vector<uint64_t> sorted_lex_factor_integers(std::begin(lex_factor_integers), std::end(lex_factor_integers));
  lex_symbol_table_ = decltype(lex_symbol_table_)(sorted_lex_factor_integers);
  return;
}

void Index::CalculateLexText
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

template <typename TextIterator, typename LexTextIterator>
LexTextIterator Index::CalculateLexTextInSlFactor
(
  TextIterator rbegin,
  TextIterator rend,
  LexTextIterator lex_text_it
)
{
  auto rfirst {rbegin};
  auto rlast {std::prev(rbegin, std::distance(rend, rbegin) % max_factor_size_)};
  if (rfirst == rlast)
  {
    rlast = std::prev(rbegin, max_factor_size_);
  }
  while (rfirst != rend)
  {
    *lex_text_it-- = lex_symbol_table_[SymbolsToInteger(std::next(rlast), std::next(rfirst))];
    rfirst = rlast;
    rlast = std::prev(rlast, max_factor_size_);
  }
  return lex_text_it;
}

void Index::CalculateColexSymbolTable (std::map<uint64_t, uint64_t> const &temp_colex_to_lex)
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

void Index::CalculateColexToLex (std::map<uint64_t, uint64_t> const &temp_colex_to_lex)
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

void Index::CalculateLexSymbolBucketOffsets
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
  for (uint64_t i {1}; i != std::size(offsets); ++i)
  {
    offsets[i] += offsets[i - 1];
  }
  lex_symbol_bucket_offsets_ = decltype(lex_symbol_bucket_offsets_)(offsets);
  return;
}

void Index::CalculateLexBwt (sdsl::int_vector<> const &lex_text)
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

template <typename Iterator>
bool Index::CalculatePattern (Iterator begin, Iterator end, sdsl::int_vector<> &pattern)
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

template <typename Iterator, typename Range>
void Index::CountWithinSlFactor
(
  Iterator rbegin,
  Iterator rend,
  uint64_t const duplicate_factor_size,
  Range &pattern_range
)
{
  uint64_t count {};
  auto size {static_cast<uint64_t>(std::distance(rend, rbegin))};
  for (uint64_t k {1}; (k <= max_factor_size_) && (k <= size); ++k)
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
        if (prefix_infix_size >= max_factor_size_)
        {
          rfirst = rlast;
          rlast = std::next(rend, prefix_infix_size % max_factor_size_);
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

template <typename Iterator>
std::pair<uint64_t, uint64_t> Index::CalculateSymbolRange
(
  Iterator begin,
  Iterator end,
  int8_t const step
)
{
  uint64_t factor_integer {};
  auto it {begin};
  auto k {max_factor_size_};
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

template <typename Iterator, typename Range>
void Index::BackwardSearchInfixFactors
(
  Iterator rbegin,
  Iterator rend,
  Range &pattern_range
)
{
  auto rfirst {rbegin};
  auto rlast {std::prev(rbegin, std::distance(rend, rbegin) % max_factor_size_)};
  if (rfirst == rlast)
  {
    rlast = std::prev(rbegin, max_factor_size_);
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
      rlast = std::prev(rlast, max_factor_size_);
    }
    else
    {
      return;
    }
  }
  return;
}

template <typename Iterator, typename Range>
void Index::BackwardSearchPrefixSlFactor
(
  Iterator rbegin,
  Iterator rend,
  Range &pattern_range
)
{
  uint64_t count {};
  auto size {static_cast<uint64_t>(std::distance(rend, rbegin))};
  for (uint64_t k {1}; (k <= max_factor_size_) && (k < size); ++k)
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
      if (prefix_infix_size >= max_factor_size_)
      {
        rfirst = rlast;
        rlast = std::next(rend, prefix_infix_size % max_factor_size_);
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
  if (size <= max_factor_size_)
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

template <typename Range>
uint64_t Index::RangeCount
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

template <typename Iterator, typename Range>
Iterator Index::BackwardSearchSuffix
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
        if ((sl_factor_size % max_factor_size_) != (suffix_s_size % max_factor_size_))
        {
          BackwardSearchSuffixSlFactor(rbegin, rlast, pattern_range_l);
        }
      }
      else
      {
        auto duplicate_factor_size {suffix_s_size % max_factor_size_};
        if (duplicate_factor_size == 0)
        {
          duplicate_factor_size = max_factor_size_;
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

template <typename Iterator>
Iterator Index::ClassifySuffix (Iterator rbegin, Iterator rend)
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

template <typename Iterator>
void Index::CalculateSlFactor
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

template <typename Iterator, typename Range>
void Index::BackwardSearchSuffixSlFactor
(
  Iterator rbegin,
  Iterator rend,
  Range &pattern_range
)
{
  auto rfirst {rbegin};
  auto rlast {std::prev(rbegin, std::distance(rend, rbegin) % max_factor_size_)};
  if (rlast == rbegin)
  {
    rlast = std::prev(rbegin, max_factor_size_);
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

template <typename Iterator, typename Range>
void Index::BackwardSearchInfix
(
  Iterator rbegin,
  Iterator rend,
  Range &pattern_range_s,
  Range &pattern_range_l
)
{
  auto rfirst {rbegin};
  auto rlast {std::prev(rbegin, std::distance(rend, rbegin) % max_factor_size_)};
  if (rfirst == rlast)
  {
    rlast = std::prev(rbegin, max_factor_size_);
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
    rlast = std::prev(rlast, max_factor_size_);
  }
  return;
}

template <typename Iterator, typename Range>
void Index::BackwardSearchPrefix
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
  for (uint64_t k {1}; (k <= max_factor_size_) && (k < size); ++k)
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
      if (prefix_infix_size >= max_factor_size_)
      {
        rfirst = rlast;
        rlast = std::next(rend, prefix_infix_size % max_factor_size_);
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
  if (size <= max_factor_size_)
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
