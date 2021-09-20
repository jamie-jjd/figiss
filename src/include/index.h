#pragma once

// #include <deque>
// #include <map>
// #include <memory>
// #include <set>

// #include <sdsl/wavelet_trees.hpp>

#include "byte_alphabet.h"
#include "grammar_xbwt_trie.h"
#include "integer_alphabet.h"
// #include "sparse_prefix_sum.h"
// #include "utility.h"

namespace figiss
{
class Index
{
public:

  static constexpr uint8_t kS {1};
  static constexpr uint8_t kL {0};

  Index () = default;
  Index (std::filesystem::path const& byte_text_path);

  // uint64_t Serialize
  // (
  //   std::filesystem::path const& index_path,
  //   std::shared_ptr<SpaceNode> root = nullptr
  // );
  // void Load (std::filesystem::path const& index_path);
  //
  // template <typename Iterator>
  // uint64_t Count (Iterator begin, Iterator end);
  //
  // friend std::ostream& operator<< (std::ostream& out, Index const& index);

private:

  void CalculateRunLengthText
  (
    sdsl::int_vector<8> const& byte_text,
    sdsl::int_vector<>& run_length_text
  );
  void InsertGrammarRulesIntoGrammarTrie
  (
    sdsl::int_vector<> const& run_length_text,
    GrammarTrie& grammar_trie
  );
  void CalculateRunLengthAlphabetAndLexRankOnGrammarTrie
  (
    IntegerAlphabet& run_length_alphabet,
    GrammarTrie& grammar_trie
  );

  // template <typename Iterator>
  // void CalculateLexTextSizeAndLexFactorIntegersAndTempColexToLex
  // (
  //   sdsl::int_vector<> const& text,
  //   uint64_t& lex_text_size,
  //   std::set<uint64_t>& lex_factor_integers,
  //   std::map<uint64_t, uint64_t>& temp_colex_to_lex
  // );
  // template <typename Iterator>
  // void CalculateLexTextSizeAndLexFactorIntegersAndTempColexToLexInSlFactor
  // (
  //   Iterator rbegin,
  //   Iterator rend,
  //   uint64_t& lex_text_size,
  //   std::set<uint64_t>& lex_factor_integers,
  //   std::map<uint64_t, uint64_t>& temp_colex_to_lex
  // );
  // void CalculateLexSymbolTable (std::set<uint64_t> const& lex_factor_integers);
  // void CalculateLexText (sdsl::int_vector<> const& text, sdsl::int_vector<>& lex_text);
  // template <typename TextIterator, typename LexTextIterator>
  // LexTextIterator CalculateLexTextInSlFactor
  // (
  //   TextIterator rbegin,
  //   TextIterator rend,
  //   LexTextIterator lex_text_it
  // );
  // void CalculateColexSymbolTable (std::map<uint64_t, uint64_t> const& temp_colex_to_lex);
  // void CalculateColexToLex (std::map<uint64_t, uint64_t> const& temp_colex_to_lex);
  // void CalculateLexSymbolBucketOffsets
  // (
  //   uint64_t const lex_text_alphabet_size,
  //   sdsl::int_vector<> const& lex_text
  // );
  // void CalculateLexBwt (sdsl::int_vector<> const& lex_text);
  //
  // template <typename Range> // [,]
  // inline bool IsNotEmptyRange (Range const& range) const
  // {
  //   return (std::get<0>(range) <= std::get<1>(range));
  // }
  //
  // template <typename Range> // [,]
  // inline uint64_t RangeSize (Range const& range) const
  // {
  //   if (std::get<0>(range) <= std::get<1>(range))
  //   {
  //     return (std::get<1>(range) - std::get<0>(range) + 1);
  //   }
  //   return 0;
  // }
  //
  // template <typename Iterator>
  // bool CalculatePattern (Iterator begin, Iterator end, sdsl::int_vector<>& pattern);
  // template <typename Iterator, typename Range>
  // Iterator BackwardSearchSuffix
  // (
  //   Iterator rbegin,
  //   Iterator rend,
  //   Range& pattern_range_s,
  //   Range& pattern_range_l
  // );
  // template <typename Iterator, typename Range>
  // void BackwardSearchInfix
  // (
  //   Iterator rbegin,
  //   Iterator rend,
  //   Range& pattern_range_s,
  //   Range& pattern_range_l
  // );
  // template <typename Iterator, typename Range>
  // void BackwardSearchPrefix
  // (
  //   Iterator rbegin,
  //   Iterator rend,
  //   Range& pattern_range_s,
  //   Range& pattern_range_l
  // );
  // template <typename Iterator>
  // Iterator ClassifySuffix (Iterator rbegin, Iterator rend);
  // template <typename Iterator>
  // void CalculateSlFactor (Iterator const rend, Iterator& rfirst, Iterator& rlast);
  // template <typename Iterator, typename Range>
  // void BackwardSearchSuffixSlFactor (Iterator rbegin, Iterator rend, Range& pattern_range);
  // template <typename Iterator, typename Range>
  // void BackwardSearchInfixFactors (Iterator rbegin, Iterator rend, Range& pattern_range);
  // template <typename Iterator, typename Range>
  // void BackwardSearchPrefixSlFactor (Iterator rbegin, Iterator rend, Range& pattern_range);
  // template <typename Iterator, typename Range>
  // void CountWithinSlFactor
  // (
  //   Iterator rbegin,
  //   Iterator rend,
  //   uint64_t const duplicate_factor_size,
  //   Range& pattern_range
  // );
  // template <typename Iterator>
  // std::pair<uint64_t, uint64_t> CalculateSymbolRange
  // (
  //   Iterator begin,
  //   Iterator end,
  //   int8_t const step = 1
  // );
  // template <typename Range>
  // uint64_t RangeCount
  // (
  //   Range const pattern_range,
  //   Range const colex_symbol_range
  // );

  ByteAlphabet byte_alphabet_;
  // GrammarXbwtTrie grammar_xbwt_trie_;
  // SymbolBucketOffsets lex_symbol_bucket_offsets_;
  // sdsl::wt_rlmn
  // <
  //   sdsl::sd_vector<>,
  //   typename sdsl::sd_vector<>::rank_1_type,
  //   typename sdsl::sd_vector<>::select_1_type,
  //   sdsl::wt_ap<>
  // >
  // lex_bwt_;

};

Index::Index (std::filesystem::path const& byte_text_path)
{
  std::cout << "construct index of " << std::filesystem::canonical(byte_text_path) << "\n";
  sdsl::int_vector<8> byte_text;
  sdsl::load_vector_from_file(byte_text, byte_text_path);
  {
    for (auto const byte : byte_text)
    {
      if (byte == 0)
      {
        throw std::runtime_error("byte_text contains 0");
      }
    }
    sdsl::append_zero_symbol(byte_text);
  }
  byte_alphabet_ = decltype(byte_alphabet_)(byte_text);
  // std::cout << byte_alphabet_;
  sdsl::int_vector<> run_length_text;
  CalculateRunLengthText(byte_text, run_length_text);
  // {
  //   auto length_width {run_length_text.width() - byte_alphabet_.GetEffectiveAlphabetWidth()};
  //   uint64_t divisor {1ULL << length_width};
  //   for (auto const& integer : run_length_text)
  //   {
  //     std::cout << "(" << (integer / divisor) << "," << (integer % divisor) << ")";
  //   }
  //   std::cout << "\n";
  // }
  // std::cout << run_length_alphabet;
  {
    auto length_width {static_cast<uint8_t>(run_length_text.width() - byte_alphabet_.GetEffectiveAlphabetWidth())};
    GrammarTrie grammar_trie {length_width};
    InsertGrammarRulesIntoGrammarTrie(run_length_text, grammar_trie);
    // std::cout << grammar_trie;
    IntegerAlphabet run_length_alphabet;
    CalculateRunLengthAlphabetAndLexRankOnGrammarTrie(run_length_alphabet, grammar_trie);
    // std::cout << grammar_trie;
    // std::cout << run_length_alphabet;
    // grammar_xbwt_trie_ = decltype(grammar_xbwt_trie_)(grammar_trie, run_length_alphabet);
    // std::cout << grammar_xbwt_trie_;
  }
  // uint64_t lex_text_size {};
  // std::set<uint64_t> lex_factor_integers;
  // std::map<uint64_t, uint64_t> temp_colex_to_lex;
  // {
  //   CalculateLexTextSizeAndLexFactorIntegersAndTempColexToLex
  //   (
  //     text,
  //     lex_text_size,
  //     lex_factor_integers,
  //     temp_colex_to_lex
  //   );
  //   // std::cout << "lex_text_size:\n" << lex_text_size << "\n";
  //   // std::cout << "|lex_factor_integers|: " << std::size(lex_factor_integers) << "\n";
  //   // std::cout << "lex_factor_integers:\n";
  //   // uint64_t i {};
  //   // for (auto const& factor : lex_factor_integers)
  //   // {
  //   //   std::cout << i++ << ":";
  //   //   Print(IntegerToSymbols(factor), std::cout);
  //   // }
  //   // std::cout << "temp_colex_to_lex:\n";
  //   // for (auto const& pair : temp_colex_to_lex)
  //   // {
  //   //   Print(IntegerToSymbols(std::get<0>(pair)), std::cout, 1, " ", "");
  //   //   std::cout << ":";
  //   //   Print(IntegerToSymbols(std::get<1>(pair)), std::cout, 1, " ", "");
  //   //   std::cout << "\n";
  //   // }
  // }
  // uint64_t lex_text_alphabet_size {std::size(lex_factor_integers)};
  // uint64_t lex_text_width {sdsl::bits::hi(lex_text_alphabet_size - 1) + 1};
  // {
  //   CalculateLexSymbolTable(lex_factor_integers);
  //   // std::cout << lex_symbol_table_;
  // }
  // sdsl::int_vector<> lex_text(lex_text_size, 0, lex_text_width);
  // {
  //   CalculateLexText(text, lex_text);
  //   // Print(lex_text, std::cout);
  // }
  // {
  //   CalculateColexSymbolTable(temp_colex_to_lex);
  //   // std::cout << colex_symbol_table_;
  // }
  // {
  //   CalculateColexToLex(temp_colex_to_lex);
  //   // Print(colex_to_lex_, std::cout);
  // }
  // {
  //   CalculateLexSymbolBucketOffsets(lex_text_alphabet_size, lex_text);
  //   // std::cout << lex_symbol_bucket_offsets_;
  // }
  // {
  //   CalculateLexBwt(lex_text);
  //   // std::cout << lex_bwt_ << "\n";
  // }
}

// uint64_t Index::Serialize
// (
//   std::filesystem::path const& index_path,
//   std::shared_ptr<SpaceNode> root
// )
// {
//   if
//   (
//     !index_path.parent_path().filename().string().empty()
//     && !std::filesystem::exists(index_path.parent_path())
//   )
//   {
//     std::filesystem::create_directories(index_path.parent_path());
//   }
//   std::fstream index_file {index_path, std::ios_base::out | std::ios_base::trunc};
//   std::cout << "serialize index to " << std::filesystem::canonical(index_path) << "\n";
//   uint64_t size_in_bytes {};
//   if (!root)
//   {
//     sdsl::write_member(kMaxFactorSize, index_file);
//     symbol_table_.Serialize(index_file);
//     sub_factor_trie_.Serialize(index_file);
//     lex_symbol_table_.Serialize(index_file);
//     colex_symbol_table_.Serialize(index_file);;
//     sdsl::serialize(colex_to_lex_, index_file);
//     lex_symbol_bucket_offsets_.Serialize(index_file);
//     sdsl::serialize(lex_bwt_, index_file);
//   }
//   else
//   {
//     root->AddLeaf("kMaxFactorSize", sdsl::write_member(kMaxFactorSize, index_file));
//     root->AccumalateSizeInBytes(symbol_table_.Serialize(index_file, root, "symbol_table_"));
//     root->AccumalateSizeInBytes(sub_factor_trie_.Serialize(index_file, root, "sub_factor_trie_"));
//     root->AccumalateSizeInBytes(lex_symbol_table_.Serialize(index_file, root, "lex_symbol_table_"));
//     root->AccumalateSizeInBytes(colex_symbol_table_.Serialize(index_file, root, "colex_symbol_table_"));
//     root->AddLeaf("colex_to_lex_", sdsl::serialize(colex_to_lex_, index_file));
//     root->AccumalateSizeInBytes(lex_symbol_bucket_offsets_.Serialize(index_file, root, "lex_symbol_bucket_offsets_"));
//     root->AddLeaf("lex_bwt_", sdsl::serialize(lex_bwt_, index_file));
//     size_in_bytes = root->GetSizeInBytes();
//   }
//   return size_in_bytes;
// }

// void Index::Load (std::filesystem::path const& index_path)
// {
//   std::ifstream index_file {index_path};
//   std::cout << "load index from " << std::filesystem::canonical(index_path) << "\n";
//   symbol_table_.Load(index_file);
//   sub_factor_trie_.Load(index_file);
//   lex_symbol_table_.Load(index_file);
//   colex_symbol_table_.Load(index_file);
//   colex_to_lex_.load(index_file);
//   lex_symbol_bucket_offsets_.Load(index_file);
//   lex_bwt_.load(index_file);
//   return;
// }

// template <typename Iterator>
// uint64_t Index::Count (Iterator begin, Iterator end)
// {
//   uint64_t count {};
//   sdsl::int_vector<> pattern;
//   if (CalculatePattern(begin, end, pattern))
//   {
//     // Print(pattern, std::cout);
//     std::pair<uint64_t, uint64_t> pattern_range_s {1, 0};
//     std::pair<uint64_t, uint64_t> pattern_range_l {1, 0};
//     auto rbegin {std::prev(std::end(pattern))};
//     auto rend {std::prev(std::begin(pattern))};
//     auto rfirst {rbegin};
//     auto rlast {rbegin};
//     rlast = BackwardSearchSuffix(rbegin, rend, pattern_range_s, pattern_range_l);
//     if (rlast != rend)
//     {
//       while (IsNotEmptyRange(pattern_range_s) || IsNotEmptyRange(pattern_range_l))
//       {
//         CalculateSlFactor(rend, rfirst, rlast);
//         if (rlast != rend)
//         {
//           BackwardSearchInfix(rfirst, rlast, pattern_range_s, pattern_range_l);
//         }
//         else
//         {
//           BackwardSearchPrefix(rfirst, rlast, pattern_range_s, pattern_range_l);
//           break;
//         }
//       }
//     }
//     {
//       count += RangeSize(pattern_range_s);
//       count += RangeSize(pattern_range_l);
//       auto size {std::distance(begin, end)};
//       if (size < (Index::kMaxFactorSize - 1))
//       {
//         count += sub_factor_trie_.Count(std::next(rend), std::next(rbegin));
//       }
//       if (size < Index::kMaxFactorSize)
//       {
//         auto colex_symbol_range {CalculateSymbolRange(rbegin, rend, -1)};
//         std::get<0>(colex_symbol_range) += (colex_symbol_table_[SymbolsToInteger(rbegin, rend, -1)] != 0);
//         if (IsNotEmptyRange(colex_symbol_range))
//         {
//           count += RangeCount({0, std::size(lex_bwt_) - 1}, colex_symbol_range);
//         }
//       }
//     }
//   }
//   return count;
// }

// std::ostream& operator<< (std::ostream& out, Index const& index)
// {
//   out << index.symbol_table_;
//   out << std::make_pair(index.sub_factor_trie_, true);
//   out << index.lex_symbol_table_;
//   out << index.colex_symbol_table_;
//   out << index.colex_to_lex_ << "\n";
//   out << index.lex_symbol_bucket_offsets_ << "\n";
//   out << index.lex_bwt_ << "\n";
//   return out;
// }

void Index::CalculateRunLengthText
(
  sdsl::int_vector<8> const& byte_text,
  sdsl::int_vector<>& run_length_text
)
{
  uint8_t prev_byte {byte_text[0]};
  uint64_t run_length_text_size {1};
  uint64_t length {};
  uint64_t max_length {1};
  for (auto const byte : byte_text)
  {
    if (prev_byte != byte)
    {
      if (max_length < length)
      {
        max_length = length;
      }
      prev_byte = byte;
      length = 0;
      ++run_length_text_size;
    }
    ++length;
  }
  {
    auto length_width {sdsl::bits::hi(max_length) + 1};
    run_length_text.width(byte_alphabet_.GetEffectiveAlphabetWidth() + length_width);
    run_length_text.resize(run_length_text_size);
    auto it {std::begin(run_length_text)};
    uint64_t divisor {1ULL << length_width};
    prev_byte = byte_text[0];
    length = 0;
    for (auto const byte : byte_text)
    {
      if (prev_byte != byte)
      {
        *it++ = (byte_alphabet_.GetRank(prev_byte) * divisor) + length;
        prev_byte = byte;
        length = 0;
      }
      ++length;
    }
    *it = 0;
  }
}

void Index::InsertGrammarRulesIntoGrammarTrie
(
  sdsl::int_vector<> const& run_length_text,
  GrammarTrie& grammar_trie
)
{
  auto text_rend {std::prev(std::begin(run_length_text))};
  auto text_it {std::prev(std::end(run_length_text), 3)};
  auto next_symbol {*std::prev(std::end(run_length_text), 2)};
  uint8_t sl_type {};
  uint8_t next_sl_type {Index::kL};
  auto sl_factor_end {std::prev(std::end(run_length_text))};
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
      grammar_trie.Insert(sl_factor_begin, sl_factor_end);
      sl_factor_end = sl_factor_begin;
    }
    next_symbol = *text_it--;
    next_sl_type = sl_type;
  }
  grammar_trie.Insert(std::begin(run_length_text), sl_factor_end);
  return;
}

void Index::CalculateRunLengthAlphabetAndLexRankOnGrammarTrie
(
  IntegerAlphabet& run_length_alphabet,
  GrammarTrie& grammar_trie
)
{
  std::set<uint64_t> alphabet;
  uint64_t lex_rank {};
  std::deque<std::shared_ptr<GrammarTrie::Node>> nodes;
  nodes.emplace_back(grammar_trie.GetRoot());
  while (!nodes.empty())
  {
    auto node {nodes.back()};
    alphabet.insert(node->label);
    if (node->lex_rank)
    {
      node->lex_rank = ++lex_rank;
    }
    nodes.pop_back();
    for (auto it {std::rbegin(node->children)}; it != std::rend(node->children); ++it)
    {
      nodes.emplace_back(std::get<1>(*it));
    }
  }
  run_length_alphabet = IntegerAlphabet(alphabet);
  return;
}

// template <typename Iterator>
// void Index::InsertSubFactorsInSlFactor (Iterator begin, Iterator end, Trie& trie)
// {
//   auto rend {std::prev(begin)};
//   auto rfirst {std::prev(end)};
//   auto rlast {std::prev(rfirst, std::distance(begin, end) % Index::kMaxFactorSize)};
//   if (rfirst == rlast)
//   {
//     rlast = std::prev(rfirst, Index::kMaxFactorSize);
//   }
//   while (rfirst != rend)
//   {
//     if (std::distance(rlast, rfirst) > 2)
//     {
//       for (auto it {std::next(rlast, 2)}; it != rfirst; ++it)
//       {
//         trie.Insert(it, rfirst);
//       }
//     }
//     rfirst = rlast;
//     rlast = std::prev(rlast, Index::kMaxFactorSize);
//   }
//   return;
// }

// void Index::CalculateLexTextSizeAndLexFactorIntegersAndTempColexToLex
// (
//   sdsl::int_vector<> const& text,
//   uint64_t& lex_text_size,
//   std::set<uint64_t>& lex_factor_integers,
//   std::map<uint64_t, uint64_t>& temp_colex_to_lex
// )
// {
//   {
//     ++lex_text_size;
//     lex_factor_integers.insert(0);
//     temp_colex_to_lex[0] = 0;
//   }
//   auto text_rend {std::prev(std::begin(text))};
//   auto text_it {std::prev(std::end(text), 3)};
//   auto next_symbol {*std::prev(std::end(text), 2)};
//   uint8_t sl_type {};
//   uint8_t next_sl_type {Index::kL};
//   auto sl_factor_end {std::prev(std::end(text))};
//   while (text_it != text_rend)
//   {
//     if (*text_it == next_symbol)
//     {
//       sl_type = next_sl_type;
//     }
//     else if (*text_it < next_symbol)
//     {
//       sl_type = Index::kS;
//     }
//     else
//     {
//       sl_type = Index::kL;
//     }
//     if ((sl_type == Index::kL) && (next_sl_type == Index::kS))
//     {
//       auto sl_factor_begin {std::next(text_it)};
//       CalculateLexTextSizeAndLexFactorIntegersAndTempColexToLexInSlFactor
//       (
//         std::prev(sl_factor_end),
//         std::prev(sl_factor_begin),
//         lex_text_size,
//         lex_factor_integers,
//         temp_colex_to_lex
//       );
//       sl_factor_end = sl_factor_begin;
//     }
//     next_symbol = *text_it--;
//     next_sl_type = sl_type;
//   }
//   CalculateLexTextSizeAndLexFactorIntegersAndTempColexToLexInSlFactor
//   (
//     std::prev(sl_factor_end),
//     text_rend,
//     lex_text_size,
//     lex_factor_integers,
//     temp_colex_to_lex
//   );
//   return;
// }

// template <typename Iterator>
// void Index::CalculateLexTextSizeAndLexFactorIntegersAndTempColexToLexInSlFactor
// (
//   Iterator rbegin,
//   Iterator rend,
//   uint64_t& lex_text_size,
//   std::set<uint64_t>& lex_factor_integers,
//   std::map<uint64_t, uint64_t>& temp_colex_to_lex
// )
// {
//   auto rfirst {rbegin};
//   auto rlast {std::prev(rbegin, std::distance(rend, rbegin) % Index::kMaxFactorSize)};
//   if (rfirst == rlast)
//   {
//     rlast = std::prev(rbegin, Index::kMaxFactorSize);
//   }
//   while (rfirst != rend)
//   {
//     ++lex_text_size;
//     auto lex_factor_integer {SymbolsToInteger(std::next(rlast), std::next(rfirst))};
//     auto colex_factor_integer {SymbolsToInteger(rfirst, rlast, -1)};
//     lex_factor_integers.insert(lex_factor_integer);
//     temp_colex_to_lex[colex_factor_integer] = lex_factor_integer;
//     rfirst = rlast;
//     rlast = std::prev(rlast, Index::kMaxFactorSize);
//   }
//   return;
// }

// void Index::CalculateLexSymbolTable (std::set<uint64_t> const& lex_factor_integers)
// {
//   std::vector<uint64_t> sorted_lex_factor_integers(std::begin(lex_factor_integers), std::end(lex_factor_integers));
//   lex_symbol_table_ = decltype(lex_symbol_table_)(sorted_lex_factor_integers);
//   return;
// }

// void Index::CalculateLexText
// (
//   sdsl::int_vector<> const& text,
//   sdsl::int_vector<>& lex_text
// )
// {
//   auto text_rend {std::prev(std::begin(text))};
//   auto text_it {std::prev(std::end(text), 3)};
//   auto next_symbol {*std::prev(std::end(text), 2)};
//   uint8_t sl_type {};
//   uint8_t next_sl_type {Index::kL};
//   auto sl_factor_end {std::prev(std::end(text))};
//   auto lex_text_it {std::prev(std::end(lex_text), 2)};
//   while (text_it != text_rend)
//   {
//     if (*text_it == next_symbol)
//     {
//       sl_type = next_sl_type;
//     }
//     else if (*text_it < next_symbol)
//     {
//       sl_type = Index::kS;
//     }
//     else
//     {
//       sl_type = Index::kL;
//     }
//     if ((sl_type == Index::kL) && (next_sl_type == Index::kS))
//     {
//       auto sl_factor_begin {std::next(text_it)};
//       lex_text_it = CalculateLexTextInSlFactor(std::prev(sl_factor_end), std::prev(sl_factor_begin), lex_text_it);
//       sl_factor_end = sl_factor_begin;
//     }
//     next_symbol = *text_it--;
//     next_sl_type = sl_type;
//   }
//   CalculateLexTextInSlFactor(std::prev(sl_factor_end), text_rend, lex_text_it);
//   return;
// }

// template <typename TextIterator, typename LexTextIterator>
// LexTextIterator Index::CalculateLexTextInSlFactor
// (
//   TextIterator rbegin,
//   TextIterator rend,
//   LexTextIterator lex_text_it
// )
// {
//   auto rfirst {rbegin};
//   auto rlast {std::prev(rbegin, std::distance(rend, rbegin) % Index::kMaxFactorSize)};
//   if (rfirst == rlast)
//   {
//     rlast = std::prev(rbegin, Index::kMaxFactorSize);
//   }
//   while (rfirst != rend)
//   {
//     *lex_text_it-- = lex_symbol_table_[SymbolsToInteger(std::next(rlast), std::next(rfirst))];
//     rfirst = rlast;
//     rlast = std::prev(rlast, Index::kMaxFactorSize);
//   }
//   return lex_text_it;
// }

// void Index::CalculateColexSymbolTable (std::map<uint64_t, uint64_t> const& temp_colex_to_lex)
// {
//   std::vector<uint64_t> sorted_colex_factor_integers(std::size(temp_colex_to_lex));
//   auto it {std::begin(sorted_colex_factor_integers)};
//   for (auto const& pair : temp_colex_to_lex)
//   {
//     *it++ = std::get<0>(pair);
//   }
//   colex_symbol_table_ = decltype(colex_symbol_table_)(sorted_colex_factor_integers);
//   return;
// }

// void Index::CalculateColexToLex (std::map<uint64_t, uint64_t> const& temp_colex_to_lex)
// {
//   colex_to_lex_.width(sdsl::bits::hi(std::size(temp_colex_to_lex) - 1) + 1);
//   colex_to_lex_.resize(std::size(temp_colex_to_lex));
//   auto it {std::begin(colex_to_lex_)};
//   for (auto const& pair : temp_colex_to_lex)
//   {
//     *it++ = lex_symbol_table_[std::get<1>(pair)];
//   }
//   return;
// }

// void Index::CalculateLexSymbolBucketOffsets
// (
//   uint64_t const lex_text_alphabet_size,
//   sdsl::int_vector<> const& lex_text
// )
// {
//   sdsl::int_vector<> offsets
//   (
//     lex_text_alphabet_size + 1,
//     0,
//     sdsl::bits::hi(std::size(lex_text)) + 1
//   );
//   for (auto const lex_symbol : lex_text)
//   {
//     ++offsets[lex_symbol];
//   }
//   std::partial_sum(std::begin(offsets), std::end(offsets), std::begin(offsets));
//   lex_symbol_bucket_offsets_ = decltype(lex_symbol_bucket_offsets_)(offsets);
//   return;
// }

// void Index::CalculateLexBwt (sdsl::int_vector<> const& lex_text)
// {
//   sdsl::int_vector<> buffer;
//   sdsl::qsufsort::construct_sa(buffer, lex_text);
//   for (auto it {std::begin(buffer)}; it != std::end(buffer); ++it)
//   {
//     if (*it != 0)
//     {
//       *it = lex_text[*it - 1];
//     }
//   }
//   sdsl::construct_im(lex_bwt_, buffer);
//   return;
// }

// template <typename Iterator>
// bool Index::CalculatePattern (Iterator begin, Iterator end, sdsl::int_vector<>& pattern)
// {
//   pattern.width(symbol_table_.GetEffectiveAlphabetWidth());
//   pattern.resize(std::distance(begin, end));
//   auto pattern_it {std::begin(pattern)};
//   for (auto it {begin}; it != end; ++it, ++pattern_it)
//   {
//     *pattern_it = symbol_table_[*it];
//     if (*pattern_it == 0)
//     {
//       return false;
//     }
//   }
//   return true;
// }

// template <typename Iterator, typename Range>
// void Index::CountWithinSlFactor
// (
//   Iterator rbegin,
//   Iterator rend,
//   uint64_t const duplicate_factor_size,
//   Range& pattern_range
// )
// {
//   uint64_t count {};
//   auto size {static_cast<uint64_t>(std::distance(rend, rbegin))};
//   for (uint64_t k {1}; (k <= Index::kMaxFactorSize) && (k <= size); ++k)
//   {
//     if (k != duplicate_factor_size)
//     {
//       auto rfirst {rbegin};
//       auto rlast {std::prev(rbegin, k)};
//       Range temp_pattern_range {1, 0};
//       {
//         auto lex_symbol_range {CalculateSymbolRange(std::next(rlast), std::next(rfirst))};
//         if (IsNotEmptyRange(lex_symbol_range))
//         {
//           temp_pattern_range =
//           {
//             lex_symbol_bucket_offsets_[std::get<0>(lex_symbol_range)],
//             lex_symbol_bucket_offsets_[std::get<1>(lex_symbol_range) + 1] - 1
//           };
//         }
//       }
//       auto is_not_counted {true};
//       if ((rlast != rend) && IsNotEmptyRange(temp_pattern_range))
//       {
//         auto prefix_infix_size {size - k};
//         if (prefix_infix_size >= Index::kMaxFactorSize)
//         {
//           rfirst = rlast;
//           rlast = std::next(rend, prefix_infix_size % Index::kMaxFactorSize);
//           BackwardSearchInfixFactors(rfirst, rlast, temp_pattern_range);
//         }
//         if ((rlast != rend) && IsNotEmptyRange(temp_pattern_range))
//         {
//           auto colex_symbol_range {CalculateSymbolRange(rlast, rend, -1)};
//           if (IsNotEmptyRange(colex_symbol_range))
//           {
//             count += RangeCount(temp_pattern_range, colex_symbol_range);
//           }
//           is_not_counted = false;
//         }
//       }
//       if (is_not_counted)
//       {
//         count += RangeSize(temp_pattern_range);
//       }
//     }
//   }
//   pattern_range = {1, count};
//   return;
// }

// template <typename Iterator>
// std::pair<uint64_t, uint64_t> Index::CalculateSymbolRange
// (
//   Iterator begin,
//   Iterator end,
//   int8_t const step
// )
// {
//   uint64_t factor_integer {};
//   auto it {begin};
//   auto k {Index::kMaxFactorSize};
//   while (it != std::prev(end, step))
//   {
//     factor_integer += *it * (1ULL << (symbol_table_.GetEffectiveAlphabetWidth() * ((k--) - 1)));
//     it += step;
//   }
//   auto base {1ULL << (symbol_table_.GetEffectiveAlphabetWidth() * (k - 1))};
//   factor_integer += *it * base;
//   if (step == -1)
//   {
//     return colex_symbol_table_.SymbolRange(factor_integer, factor_integer + base - 1);
//   }
//   return lex_symbol_table_.SymbolRange(factor_integer, factor_integer + base - 1);
// }

// template <typename Iterator, typename Range>
// void Index::BackwardSearchInfixFactors
// (
//   Iterator rbegin,
//   Iterator rend,
//   Range& pattern_range
// )
// {
//   auto rfirst {rbegin};
//   auto rlast {std::prev(rbegin, std::distance(rend, rbegin) % Index::kMaxFactorSize)};
//   if (rfirst == rlast)
//   {
//     rlast = std::prev(rbegin, Index::kMaxFactorSize);
//   }
//   while (rfirst != rend)
//   {
//     auto lex_symbol {lex_symbol_table_[SymbolsToInteger(std::next(rlast), std::next(rfirst))]};
//     if (lex_symbol != 0)
//     {
//       auto begin_offset {lex_symbol_bucket_offsets_[lex_symbol]};
//       pattern_range =
//       {
//         begin_offset + lex_bwt_.rank(std::get<0>(pattern_range), lex_symbol),
//         begin_offset + lex_bwt_.rank(std::get<1>(pattern_range) + 1, lex_symbol) - 1
//       };
//     }
//     else
//     {
//       pattern_range = {1, 0};
//     }
//     if (IsNotEmptyRange(pattern_range))
//     {
//       rfirst = rlast;
//       rlast = std::prev(rlast, Index::kMaxFactorSize);
//     }
//     else
//     {
//       return;
//     }
//   }
//   return;
// }

// template <typename Iterator, typename Range>
// void Index::BackwardSearchPrefixSlFactor
// (
//   Iterator rbegin,
//   Iterator rend,
//   Range& pattern_range
// )
// {
//   uint64_t count {};
//   auto size {static_cast<uint64_t>(std::distance(rend, rbegin))};
//   for (uint64_t k {1}; (k <= Index::kMaxFactorSize) && (k < size); ++k)
//   {
//     auto temp_pattern_range {pattern_range};
//     auto rfirst {rbegin};
//     auto rlast {std::prev(rbegin, k)};
//     auto lex_symbol {lex_symbol_table_[SymbolsToInteger(std::next(rlast), std::next(rfirst))]};
//     if (lex_symbol != 0)
//     {
//       auto begin_offset {lex_symbol_bucket_offsets_[lex_symbol]};
//       temp_pattern_range =
//       {
//         begin_offset + lex_bwt_.rank(std::get<0>(temp_pattern_range), lex_symbol),
//         begin_offset + lex_bwt_.rank(std::get<1>(temp_pattern_range) + 1, lex_symbol) - 1
//       };
//     }
//     else
//     {
//       temp_pattern_range = {1, 0};
//     }
//     auto is_not_counted {true};
//     if (IsNotEmptyRange(temp_pattern_range))
//     {
//       auto prefix_infix_size {size - k};
//       if (prefix_infix_size >= Index::kMaxFactorSize)
//       {
//         rfirst = rlast;
//         rlast = std::next(rend, prefix_infix_size % Index::kMaxFactorSize);
//         BackwardSearchInfixFactors(rfirst, rlast, temp_pattern_range);
//       }
//       if ((rlast != rend) && IsNotEmptyRange(temp_pattern_range))
//       {
//         auto colex_symbol_range {CalculateSymbolRange(rlast, rend, -1)};
//         if (IsNotEmptyRange(colex_symbol_range))
//         {
//           count += RangeCount(temp_pattern_range, colex_symbol_range);
//         }
//         is_not_counted = false;
//       }
//     }
//     if (is_not_counted)
//     {
//       count += RangeSize(temp_pattern_range);
//     }
//   }
//   if (size <= Index::kMaxFactorSize)
//   {
//     auto temp_pattern_range {pattern_range};
//     auto colex_symbol_range {CalculateSymbolRange(rbegin, rend, -1)};
//     if (IsNotEmptyRange(colex_symbol_range))
//     {
//       count += RangeCount(temp_pattern_range, colex_symbol_range);
//     }
//   }
//   pattern_range = {1, count};
//   return;
// }

// template <typename Range>
// uint64_t Index::RangeCount
// (
//   Range const pattern_range,
//   Range const colex_symbol_range
// )
// {
//   uint64_t count {};
//   for
//   (
//     auto colex_symbol {std::get<0>(colex_symbol_range)};
//     colex_symbol <= std::get<1>(colex_symbol_range);
//     ++colex_symbol
//   )
//   {
//     count +=
//     (
//       lex_bwt_.rank(std::get<1>(pattern_range) + 1, colex_to_lex_[colex_symbol])
//       - lex_bwt_.rank(std::get<0>(pattern_range), colex_to_lex_[colex_symbol])
//     );
//   }
//   return count;
// }

// template <typename Iterator, typename Range>
// Iterator Index::BackwardSearchSuffix
// (
//   Iterator rbegin,
//   Iterator rend,
//   Range& pattern_range_s,
//   Range& pattern_range_l
// )
// {
//   auto rfirst {rbegin};
//   auto rlast {ClassifySuffix(rbegin, rend)};
//   if (rlast != rend)
//   {
//     CalculateSlFactor(rend, rfirst, rlast);
//     auto sl_factor_size {static_cast<uint64_t>(std::distance(rlast, rbegin))};
//     auto suffix_s_size {static_cast<uint64_t>(std::distance(rfirst, rbegin))};
//     if (rfirst != rbegin)
//     {
//       BackwardSearchSuffixSlFactor(rbegin, rfirst, pattern_range_s);
//       if (IsNotEmptyRange(pattern_range_s))
//       {
//         if (rlast != rend)
//         {
//           BackwardSearchInfixFactors(rfirst, rlast, pattern_range_s);
//         }
//         else
//         {
//           BackwardSearchPrefixSlFactor(rfirst, rlast, pattern_range_s);
//         }
//       }
//       if (rlast != rend)
//       {
//         if ((sl_factor_size % Index::kMaxFactorSize) != (suffix_s_size % Index::kMaxFactorSize))
//         {
//           BackwardSearchSuffixSlFactor(rbegin, rlast, pattern_range_l);
//         }
//       }
//       else
//       {
//         auto duplicate_factor_size {suffix_s_size % Index::kMaxFactorSize};
//         if (duplicate_factor_size == 0)
//         {
//           duplicate_factor_size = Index::kMaxFactorSize;
//         }
//         CountWithinSlFactor(rbegin, rend, duplicate_factor_size, pattern_range_l);
//       }
//     }
//     else
//     {
//       if (rlast != rend)
//       {
//         BackwardSearchSuffixSlFactor(rbegin, rlast, pattern_range_l);
//       }
//       else
//       {
//         CountWithinSlFactor(rbegin, rend, 0, pattern_range_l);
//       }
//     }
//   }
//   else
//   {
//     CountWithinSlFactor(rbegin, rend, 0, pattern_range_s);
//   }
//   return rlast;
// }

// template <typename Iterator>
// Iterator Index::ClassifySuffix (Iterator rbegin, Iterator rend)
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

// template <typename Iterator>
// void Index::CalculateSlFactor
// (
//   Iterator const rend,
//   Iterator& rfirst,
//   Iterator& rlast
// )
// {
//   auto next_sl_type {Index::kL};
//   rfirst = rlast--;
//   while ((rlast != rend) && !((next_sl_type == Index::kS) && (*rlast > *std::next(rlast))))
//   {
//     if((next_sl_type == Index::kL) && (*rlast < *std::next(rlast)))
//     {
//       next_sl_type = Index::kS;
//     }
//     --rlast;
//   }
//   return;
// }

// template <typename Iterator, typename Range>
// void Index::BackwardSearchSuffixSlFactor
// (
//   Iterator rbegin,
//   Iterator rend,
//   Range& pattern_range
// )
// {
//   auto rfirst {rbegin};
//   auto rlast {std::prev(rbegin, std::distance(rend, rbegin) % Index::kMaxFactorSize)};
//   if (rlast == rbegin)
//   {
//     rlast = std::prev(rbegin, Index::kMaxFactorSize);
//   }
//   auto lex_symbol_range {CalculateSymbolRange(std::next(rlast), std::next(rfirst))};
//   if (IsNotEmptyRange(lex_symbol_range))
//   {
//     pattern_range =
//     {
//       lex_symbol_bucket_offsets_[std::get<0>(lex_symbol_range)],
//       lex_symbol_bucket_offsets_[std::get<1>(lex_symbol_range) + 1] - 1
//     };
//     if (rlast != rend)
//     {
//       BackwardSearchInfixFactors(rlast, rend, pattern_range);
//     }
//   }
//   return;
// }

// template <typename Iterator, typename Range>
// void Index::BackwardSearchInfix
// (
//   Iterator rbegin,
//   Iterator rend,
//   Range& pattern_range_s,
//   Range& pattern_range_l
// )
// {
//   auto rfirst {rbegin};
//   auto rlast {std::prev(rbegin, std::distance(rend, rbegin) % Index::kMaxFactorSize)};
//   if (rfirst == rlast)
//   {
//     rlast = std::prev(rbegin, Index::kMaxFactorSize);
//   }
//   while ((rfirst != rend) && (IsNotEmptyRange(pattern_range_s) || IsNotEmptyRange(pattern_range_l)))
//   {
//     auto lex_symbol {lex_symbol_table_[SymbolsToInteger(std::next(rlast), std::next(rfirst))]};
//     if (lex_symbol != 0)
//     {
//       auto begin_offset {lex_symbol_bucket_offsets_[lex_symbol]};
//       if (IsNotEmptyRange(pattern_range_s))
//       {
//         pattern_range_s =
//         {
//           begin_offset + lex_bwt_.rank(std::get<0>(pattern_range_s), lex_symbol),
//           begin_offset + lex_bwt_.rank(std::get<1>(pattern_range_s) + 1, lex_symbol) - 1
//         };
//       }
//       if (IsNotEmptyRange(pattern_range_l))
//       {
//         pattern_range_l =
//         {
//           begin_offset + lex_bwt_.rank(std::get<0>(pattern_range_l), lex_symbol),
//           begin_offset + lex_bwt_.rank(std::get<1>(pattern_range_l) + 1, lex_symbol) - 1
//         };
//       }
//     }
//     else
//     {
//       pattern_range_s = pattern_range_l = {1, 0};
//     }
//     rfirst = rlast;
//     rlast = std::prev(rlast, Index::kMaxFactorSize);
//   }
//   return;
// }

// template <typename Iterator, typename Range>
// void Index::BackwardSearchPrefix
// (
//   Iterator rbegin,
//   Iterator rend,
//   Range& pattern_range_s,
//   Range& pattern_range_l
// )
// {
//   uint64_t count_s {};
//   uint64_t count_l {};
//   auto size {static_cast<uint64_t>(std::distance(rend, rbegin))};
//   for (uint64_t k {1}; (k <= Index::kMaxFactorSize) && (k < size); ++k)
//   {
//     auto temp_pattern_range_s {pattern_range_s};
//     auto temp_pattern_range_l {pattern_range_l};
//     auto rfirst {rbegin};
//     auto rlast {std::prev(rbegin, k)};
//     auto lex_symbol {lex_symbol_table_[SymbolsToInteger(std::next(rlast), std::next(rfirst))]};
//     if (lex_symbol != 0)
//     {
//       auto begin_offset {lex_symbol_bucket_offsets_[lex_symbol]};
//       if (IsNotEmptyRange(temp_pattern_range_s))
//       {
//         temp_pattern_range_s =
//         {
//           begin_offset + lex_bwt_.rank(std::get<0>(temp_pattern_range_s), lex_symbol),
//           begin_offset + lex_bwt_.rank(std::get<1>(temp_pattern_range_s) + 1, lex_symbol) - 1
//         };
//       }
//       if (IsNotEmptyRange(temp_pattern_range_l))
//       {
//         temp_pattern_range_l =
//         {
//           begin_offset + lex_bwt_.rank(std::get<0>(temp_pattern_range_l), lex_symbol),
//           begin_offset + lex_bwt_.rank(std::get<1>(temp_pattern_range_l) + 1, lex_symbol) - 1
//         };
//       }
//     }
//     else
//     {
//       temp_pattern_range_s = temp_pattern_range_l = {1, 0};
//     }
//     auto is_not_counted {true};
//     if (IsNotEmptyRange(temp_pattern_range_s) || IsNotEmptyRange(temp_pattern_range_l))
//     {
//       auto prefix_infix_size {size - k};
//       if (prefix_infix_size >= Index::kMaxFactorSize)
//       {
//         rfirst = rlast;
//         rlast = std::next(rend, prefix_infix_size % Index::kMaxFactorSize);
//         BackwardSearchInfix(rfirst, rlast, temp_pattern_range_s, temp_pattern_range_l);
//       }
//       if ((rlast != rend) && (IsNotEmptyRange(temp_pattern_range_s) || IsNotEmptyRange(temp_pattern_range_l)))
//       {
//         auto colex_symbol_range {CalculateSymbolRange(rlast, rend, -1)};
//         if (IsNotEmptyRange(colex_symbol_range))
//         {
//           if (IsNotEmptyRange(temp_pattern_range_s))
//           {
//             auto pattern_range {temp_pattern_range_s};
//             count_s += RangeCount(pattern_range, colex_symbol_range);
//           }
//           if (IsNotEmptyRange(temp_pattern_range_l))
//           {
//             auto pattern_range {temp_pattern_range_l};
//             count_l += RangeCount(pattern_range, colex_symbol_range);
//           }
//         }
//         is_not_counted = false;
//       }
//       if (is_not_counted)
//       {
//         count_s += RangeSize(temp_pattern_range_s);
//         count_l += RangeSize(temp_pattern_range_l);
//       }
//     }
//   }
//   if (size <= Index::kMaxFactorSize)
//   {
//     auto colex_symbol_range {CalculateSymbolRange(rbegin, rend, -1)};
//     if (IsNotEmptyRange(colex_symbol_range))
//     {
//       if (IsNotEmptyRange(pattern_range_s))
//       {
//         auto pattern_range {pattern_range_s};
//         count_s += RangeCount(pattern_range, colex_symbol_range);
//       }
//       if (IsNotEmptyRange(pattern_range_l))
//       {
//         auto pattern_range {pattern_range_l};
//         count_l += RangeCount(pattern_range, colex_symbol_range);
//       }
//     }
//   }
//   pattern_range_s = {1, count_s};
//   pattern_range_l = {1, count_l};
//   return;
// }
}
