#pragma once

#include <deque>
#include <map>
#include <memory>

#include <sdsl/wavelet_trees.hpp>

#include "utility.h"

namespace project
{
template <typename T>
void f (T)
{
  std::cout << __PRETTY_FUNCTION__ << "\n";
  return;
}
class ByteAlphabet
{
public:

  ByteAlphabet () = default;
  ByteAlphabet (ByteAlphabet const &) = default;
  ByteAlphabet (ByteAlphabet &&) = default;
  ByteAlphabet (sdsl::int_vector<8> const &byte_text);
  ByteAlphabet& operator= (ByteAlphabet const &) = default;
  ByteAlphabet& operator= (ByteAlphabet &&) noexcept;

  inline uint64_t ToSymbol (uint64_t const byte) const
  {
    if ((byte < std::size(alphabet_bits_)) && alphabet_bits_[byte])
    {
      return byte_to_symbol_(byte);
    }
    return 0;
  }

  inline uint64_t ToByte (uint64_t const symbol) const
  {
    if (symbol < effective_alphabet_size_)
    {
      return symbol_to_byte_(symbol + 1);
    }
    return 0;
  }

  friend std::ostream& operator<< (std::ostream &out, ByteAlphabet const &byte_alphabet)
  {
    {
      out << "value:\n";
      out << "effective_alphabet_size_:\n";
      out << static_cast<uint64_t>(byte_alphabet.effective_alphabet_size_) << "\n";
      out << "effective_alphabet_width_:\n";
      out << static_cast<uint64_t>(byte_alphabet.effective_alphabet_width_) << "\n";
      out << "byte alphabet:\n";
      for (uint16_t symbol {}; symbol != byte_alphabet.effective_alphabet_size_; ++symbol)
      {
        out << byte_alphabet.ToByte(symbol);
        out << ((symbol != (byte_alphabet.effective_alphabet_size_ - 1)) ? " " : "\n");
      }
    }
    {
      out << "space:\n";
      out << "effective_alphabet_size_: " << sizeof(byte_alphabet.effective_alphabet_size_) << "B\n";
      out << "effective_alphabet_width_: " << sizeof(byte_alphabet.effective_alphabet_width_) << "B\n";
      out << "alphabet_bits_: " << ProperSizeRepresentation(sdsl::size_in_bytes(byte_alphabet.alphabet_bits_)) << "B\n";
      out << "byte_to_symbol_: " << ProperSizeRepresentation(sdsl::size_in_bytes(byte_alphabet.byte_to_symbol_)) << "B\n";
      out << "symbol_to_byte_: " << ProperSizeRepresentation(sdsl::size_in_bytes(byte_alphabet.symbol_to_byte_)) << "B\n";
    }
    return out;
  }

private:

  uint16_t effective_alphabet_size_;
  uint8_t effective_alphabet_width_;
  sdsl::bit_vector alphabet_bits_;
  sdsl::bit_vector::rank_1_type byte_to_symbol_;
  sdsl::bit_vector::select_1_type symbol_to_byte_;

};

ByteAlphabet::ByteAlphabet (sdsl::int_vector<8> const &byte_text)
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

ByteAlphabet& ByteAlphabet::operator= (ByteAlphabet &&byte_alphabet) noexcept
{
  if (this != &byte_alphabet)
  {
    effective_alphabet_size_ = std::move(byte_alphabet.effective_alphabet_size_);
    effective_alphabet_width_ = std::move(byte_alphabet.effective_alphabet_width_);
    alphabet_bits_ = std::move(byte_alphabet.alphabet_bits_);
    byte_to_symbol_ = std::move(byte_alphabet.byte_to_symbol_);
    byte_to_symbol_.set_vector(&alphabet_bits_);
    symbol_to_byte_ = std::move(byte_alphabet.symbol_to_byte_);
    symbol_to_byte_.set_vector(&alphabet_bits_);
  }
  return *this;
}

// class CompactTrie
// {
// public:
//
//   class Node
//   {
//   public:
//     std::deque<uint8_t> edge_labels;
//     std::map<uint8_t, std::shared_ptr<Node>> children;
//     uint64_t count;
//   };
//
//   std::shared_ptr<Node> root;
//
//   CompactTrie (): root {std::make_shared<Node>()} {}
//
//   template <typename Iterator>
//   void Insert (Iterator it, Iterator last)
//   {
//     auto node {root};
//     while (true)
//     {
//       auto symbol {*it++};
//       auto children_it {node->children.find(symbol)};
//       if (children_it == std::end(node->children))
//       {
//         node = node->children[symbol] = std::make_shared<Node>();
//         while (it != last)
//         {
//           node->edge_labels.emplace_back(*it);
//           ++it;
//         }
//         node->count = 1;
//         return;
//       }
//       else
//       {
//         auto child {std::get<1>(*children_it)};
//         auto labels_it {std::begin(child->edge_labels)};
//         auto labels_end {std::end(child->edge_labels)};
//         while ((it != last) && (labels_it != labels_end) && (*it == *labels_it))
//         {
//           ++it;
//           ++labels_it;
//         }
//         if (labels_it == labels_end)
//         {
//           if (it != last)
//           {
//             node = child;
//           }
//           else
//           {
//             ++(child->count);
//             return;
//           }
//         }
//         else
//         {
//           auto internal_node {std::make_shared<Node>()};
//           {
//             auto length {std::distance(std::begin(child->edge_labels), labels_it)};
//             while (length--)
//             {
//               internal_node->edge_labels.emplace_back(child->edge_labels.front());
//               child->edge_labels.pop_front();
//             }
//             internal_node->count = child->count;
//           }
//           std::get<1>(*children_it) = internal_node;
//           internal_node->children[*labels_it] = child;
//           if (it != last)
//           {
//             node = internal_node;
//           }
//           else
//           {
//             ++(internal_node->count);
//             return;
//           }
//         }
//       }
//     }
//     return;
//   }
//
//   friend std::ostream& operator<< (std::ostream &out, std::pair<CompactTrie, bool> const &pair)
//   {
//     auto &trie {std::get<0>(pair)};
//     auto is_preorder {std::get<1>(pair)};
//     std::deque<std::pair<std::shared_ptr<Node>, uint64_t>> nodes;
//     nodes.emplace_back(trie.root, 0);
//     while (!nodes.empty())
//     {
//       auto node {std::get<0>(nodes.back())};
//       auto depth {std::get<1>(nodes.back())};
//       if (is_preorder)
//       {
//         nodes.pop_back();
//       }
//       else
//       {
//         node = std::get<0>(nodes.front());
//         depth = std::get<1>(nodes.front());
//         nodes.pop_front();
//       }
//       if (depth != 0)
//       {
//
//       }
//     }
//     return out;
//   }
// };
//
// template <typename Text>
// void CalculateTemporaryTinyPatternTrie
// (
//   Text const &text,
//   CompactTrie &trie
// )
// {
//   auto text_it {std::begin(text)};
//   while (std::distance(text_it, std::end(text)) >= Index::max_factor_size)
//   {
//     auto it {text_it};
//     auto end {std::next(it, Index::max_factor_size - 1)};
//     while (it != end)
//     {
//       trie.Insert(it, end);
//       ++it;
//     }
//     ++text_it;
//   }
//   while (text_it != std::end(text))
//   {
//     trie.Insert(text_it, std::end(text));
//     ++text_it;
//   }
//   return;
// }
//
// struct TinyPatternTrie
// {
//   sdsl::bit_vector branch_bits;
//   sdsl::bit_vector::select_1_type branches_select_1;
//   sdsl::int_vector<> branch_labels;
//   sdsl::bit_vector edge_bits;
//   sdsl::bit_vector::select_1_type edge_select_1;
//   sdsl::int_vector<> edge_labels;
//   sdsl::int_vector<> counts;
//
//   void Construct (CompactTrie const &trie)
//   {
//     std::deque<uint8_t> branch_bits_;
//     std::deque<uint8_t> branch_labels_;
//     std::deque<uint8_t> edge_bits_;
//     std::deque<uint8_t> edge_labels_;
//     std::deque<uint64_t> counts_;
//     std::deque<std::shared_ptr<CompactTrie::Node>> nodes;
//     nodes.emplace_back(trie.root);
//     while (!nodes.empty())
//     {
//       auto node {nodes.front()};
//       nodes.pop_front();
//       for (auto const &pair : node->children)
//       {
//         auto label {std::get<0>(pair)};
//         auto child {std::get<1>(pair)};
//         nodes.emplace_back(child);
//         branch_bits_.emplace_back(0);
//         branch_labels_.emplace_back(label);
//         {
//           auto it {std::begin(child->edge_labels)};
//           auto prev_end {std::prev(std::end(child->edge_labels))};
//           edge_bits_.emplace_back(1);
//           while (it != prev_end)
//           {
//             edge_bits_.emplace_back(0);
//             edge_labels_.emplace_back(*it);
//           }
//           edge_labels_.emplace_back(*prev_end);
//         }
//         counts_.emplace_back(child->count);
//       }
//       branch_bits_.emplace_back(1);
//     }
//     {
//       branch_bits.resize(std::size(branch_bits_));
//       auto it {std::begin(branch_bits)};
//       for (auto const bit : branch_bits_)
//       {
//         *it++ = bit;
//       }
//       branches_select_1 = decltype(branches_select_1)(&branch_bits);
//     }
//     {
//       branch_labels.width(8);
//       branch_labels.resize(std::size(branch_labels_));
//       std::copy(std::begin(branch_labels_), std::end(branch_labels_), std::begin(branch_labels));
//       sdsl::util::bit_compress(branch_labels);
//     }
//     {
//       edge_bits.resize(std::size(edge_bits_));
//       auto it {std::begin(edge_bits)};
//       for (auto const bit : edge_bits_)
//       {
//         *it++ = bit;
//       }
//       edge_select_1 = decltype(edge_select_1)(&edge_bits);
//     }
//     {
//       edge_labels.width(8);
//       edge_labels.resize(std::size(edge_labels_));
//       std::copy(std::begin(edge_labels_), std::end(edge_labels_), std::begin(edge_labels));
//       sdsl::util::bit_compress(edge_labels);
//     }
//     {
//       counts.resize(std::size(counts_));
//       std::copy(std::begin(counts_), std::end(counts_), std::begin(counts));
//       sdsl::util::bit_compress(counts);
//     }
//     return;
//   }
//
//   template <typename File>
//   void Print (File &file)
//   {
//     file << "space:\n";
//     file << "branch_bits: " << ProperSizeRepresentation(sdsl::size_in_bytes(branch_bits)) << "B\n";
//     file << "select_1: " << ProperSizeRepresentation(sdsl::size_in_bytes(branches_select_1)) << "B\n";
//     file << "branch_labels: " << ProperSizeRepresentation(sdsl::size_in_bytes(branch_labels)) << "B\n";
//     file << "edge_bits: " << ProperSizeRepresentation(sdsl::size_in_bytes(edge_bits)) << "B\n";
//     file << "edge_select_1: " << ProperSizeRepresentation(sdsl::size_in_bytes(edge_select_1)) << "B\n";
//     file << "edge_labels: " << ProperSizeRepresentation(sdsl::size_in_bytes(edge_labels)) << "B\n";
//     file << "counts: " << ProperSizeRepresentation(sdsl::size_in_bytes(counts)) << "B\n";
//     return;
//   }
// };
//
// struct FactorTable
// {
//   sdsl::sd_vector<> bits;
//   sdsl::sd_vector<>::rank_1_type rank_1;
//
//   template <typename File>
//   void Print (File &file)
//   {
//     file << "|factors|: " << rank_1(std::size(bits)) << "\n";
//     file << "space:\n";
//     file << "bits: " << ProperSizeRepresentation(sdsl::size_in_bytes(bits)) << "B\n";
//   }
// };
//
// template <typename StringIterator>
// uint64_t SymbolsToInteger
// (
//   StringIterator it,
//   StringIterator end,
//   uint8_t const width,
//   int8_t const step = 1
// )
// {
//   uint64_t result {};
//   for (uint64_t i {8}; (it != end) && (i != 0); it += step, --i)
//   {
//     result += *it * (1ULL << (width * (i - 1)));
//   }
//   return result;
// }
//
// template <bool Index::max_factor_size>
// template
// <
//   typename Text,
//   typename LexFactorTable,
//   typename LexText
// >
// void Index::CalculateLexText
// (
//   Text const &text,
//   LexFactorTable const &lex_factor_table,
//   LexText &lex_text
// )
// {
//   auto text_prev_begin {std::prev(std::begin(text))};
//   auto text_it {std::prev(std::end(text), 2)};
//   auto next_symbol {*std::prev(std::end(text))};
//   uint8_t sl_type {};
//   uint8_t next_sl_type {Index::L};
//   auto sl_factor_end {std::end(text)};
//   auto lex_text_it {std::prev(std::end(lex_text), 2)};
//   while (text_it != text_prev_begin)
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
//       sl_type = Index::kS;
//     }
//     if ((sl_type == Index::kS) && (next_sl_type == S))
//     {
//       auto factor_it {std::next(text_it)};
//       while (std::distance(factor_it, sl_factor_end) >= Index::max_factor_size)
//       {
//         auto factor_end {std::next(factor_it, Index::max_factor_size)};
//         *lex_text_it-- = lex_factor_table.rank_1(SymbolsToInteger(factor_it, factor_end, text.width()));
//         factor_it = factor_end;
//       }
//       if (std::distance(factor_it, sl_factor_end) != 0)
//       {
//         *lex_text_it-- = lex_factor_table.rank_1(SymbolsToInteger(factor_it, sl_factor_end, text.width()));
//       }
//       sl_factor_end = std::next(text_it);
//     }
//     next_symbol = *text_it--;
//     next_sl_type = sl_type;
//   }
//   *lex_text_it = lex_factor_table.rank_1(SymbolsToInteger(std::next(text_prev_begin), sl_factor_end, text.width()));
//   return;
// }

template <uint8_t max_factor_size = 4>
class Index
{
public:

  static constexpr uint8_t kMaxFactorSize {Index::max_factor_size};
  static constexpr uint8_t kS {1};
  static constexpr uint8_t kL {0};

  Index (std::filesystem::path const &byte_text_path);
  uint64_t Serialize (std::filesystem::path const &index_path);
  void Load (std::filesystem::path const &index_path);
  template <typename PatternIterator>
  uint64_t Count (PatternIterator begin, PatternIterator end);

private:

  ByteAlphabet byte_alphabet;
  // TinyPatternTrie tiny_pattern_trie;
  // FactorTable lex_factor_table;
  // FactorTable colex_factor_table;
  sdsl::int_vector<> colex_to_lex;
  sdsl::int_vector<> bucket_begin_offsets;
  sdsl::wt_rlmn
  <
    sdsl::sd_vector<>,
    typename sdsl::sd_vector<>::rank_1_type,
    typename sdsl::sd_vector<>::select_1_type,
    sdsl::wt_ap<>
  >
  bwt;

  // CalculateFactorInformation()

};

template <uint8_t max_factor_size>
Index<max_factor_size>::Index (std::filesystem::path const &byte_text_path)
{
  std::cout << "construct index of " << std::filesystem::canonical(byte_text_path) << "\n";
  sdsl::int_vector<> text;
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
    byte_alphabet = decltype(byte_alphabet)(byte_text);
    // std::cout << byte_alphabet;
    return;
    // text.width(index.byte_alphabet.symbol_width);
    // text.resize(std::size(byte_text));
    // auto text_it {std::begin(text)};
    // for (auto const &byte : byte_text)
    // {
    //   *text_it++ = index.byte_alphabet.GetSymbol(byte);
    // }
    // Print(text, std::cout, 1, "");
  }
  // {
  //   CompactTrie trie;
  //   CalculateTemporaryTinyPatternTrie(text, Index::max_factor_size, trie);
  //   // index.tiny_pattern_trie.Construct(trie);
  //   // index.tiny_pattern_trie.Print(std::cout);
  // }
  // return;
  // std::map<uint64_t, uint64_t> lex_factor_counts;
  // std::map<uint64_t, uint64_t> colex_to_lex_factor;
  // {
  //   lex_factor_counts[0] = 1;
  //   colex_to_lex_factor[0] = 0;
  //   CalculateFactorInformation(text, Index::max_factor_size, lex_factor_counts, colex_to_lex_factor);
  // }
  // sdsl::int_vector<> lex_text;
  // uint64_t lex_text_size {};
  // uint64_t lex_alphabet_size {};
  // uint64_t lex_text_width {};
  // {
  //   std::vector<uint64_t> lex_factors;
  //   for (auto const &pair : lex_factor_counts)
  //   {
  //     lex_factors.emplace_back(std::get<0>(pair));
  //     lex_text_size += std::get<1>(pair);
  //   }
  //   {
  //     lex_alphabet_size = std::size(lex_factors);
  //     lex_text_width = sdsl::bits::hi(lex_alphabet_size - 1) + 1;
  //     index.lex_factor_table.bits = decltype(index.lex_factor_table.bits)(std::begin(lex_factors), std::end(lex_factors));
  //     index.lex_factor_table.rank_1.set_vector(&(index.lex_factor_table.bits));
  //     // index.lex_factor_table.Print(std::cout);
  //   }
  //   {
  //     lex_text.width(sdsl::bits::hi(std::size(lex_factors) - 1) + 1);
  //     lex_text.resize(lex_text_size);
  //     CalculateLexText(text, index.lex_factor_table, Index::max_factor_size, lex_text);
  //     *std::prev(std::end(lex_text)) = 0;
  //     // Print(lex_text, std::cout);
  //   }
  // }
  // {
  //   std::vector<uint64_t> colex_factors;
  //   for (auto const &pair : colex_to_lex_factor)
  //   {
  //     colex_factors.emplace_back(std::get<0>(pair));
  //   }
  //   index.colex_factor_table.bits = decltype(index.colex_factor_table.bits)(std::begin(colex_factors), std::end(colex_factors));
  //   index.colex_factor_table.rank_1.set_vector(&(index.colex_factor_table.bits));
  //   // index.colex_factor_table.Print(std::cout);
  // }
  // {
  //   index.colex_to_lex.width(lex_text_width);
  //   index.colex_to_lex.resize(lex_alphabet_size);
  //   auto it {std::begin(index.colex_to_lex)};
  //   for (auto const &pair : colex_to_lex_factor)
  //   {
  //     *it++ = index.lex_factor_table.rank_1(std::get<1>(pair));
  //   }
  //   // Print(index.colex_to_lex, std::cout);
  // }
  // {
  //   index.bucket_begin_offsets.width(sdsl::bits::hi(lex_text_size) + 1);
  //   index.bucket_begin_offsets.resize(lex_alphabet_size + 1);
  //   sdsl::util::set_to_value(index.bucket_begin_offsets, 0);
  //   for (auto const lex_symbol : lex_text)
  //   {
  //     ++(index.bucket_begin_offsets[lex_symbol]);
  //   }
  //   std::partial_sum
  //   (
  //     std::begin(index.bucket_begin_offsets),
  //     std::end(index.bucket_begin_offsets),
  //     std::begin(index.bucket_begin_offsets)
  //   );
  //   auto it {std::prev(std::end(index.bucket_begin_offsets))};
  //   auto begin {std::begin(index.bucket_begin_offsets)};
  //   while (it != begin)
  //   {
  //     *it-- = *std::prev(it);
  //   }
  //   *begin = 0;
  //   // Print(index.bucket_begin_offsets, std::cout);
  // }
  // {
  //   sdsl::int_vector<> buffer;
  //   sdsl::qsufsort::construct_sa(buffer, lex_text);
  //   for (auto it {std::begin(buffer)}; it != std::end(buffer); ++it)
  //   {
  //     if (*it != 0)
  //     {
  //       *it = lex_text[(*it - 1)];
  //     }
  //   }
  //   sdsl::construct_im(index.bwt, buffer);
  //   // Print(index.bwt, std::cout);
  // }
}

// template
// <
//   typename Text,
//   typename LexFactorCounts,
//   typename ColexToLexFactor
// >
// void CalculateFactorInformation
// (
//   Text const &text,
//   uint8_t const Index::max_factor_size,
//   LexFactorCounts &lex_factor_counts,
//   ColexToLexFactor &colex_to_lex_factor
// )
// {
//   auto text_prev_begin {std::prev(std::begin(text))};
//   auto text_it {std::prev(std::end(text), 2)};
//   auto next_symbol {*std::prev(std::end(text))};
//   uint8_t sl_type {};
//   uint8_t next_sl_type {Index::kL};
//   auto sl_factor_end {std::end(text)};
//   uint64_t lex_factor;
//   uint64_t colex_factor;
//   while (text_it != text_prev_begin)
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
//       sl_type = L;
//     }
//     if ((sl_type == L) && (next_sl_type == Index::kS))
//     {
//       auto factor_it {std::next(text_it)};
//       uint64_t lex_factor {};
//       uint64_t colex_factor {};
//       while (std::distance(factor_it, sl_factor_end) >= Index::max_factor_size)
//       {
//         auto factor_end {std::next(factor_it, Index::max_factor_size)};
//         lex_factor = SymbolsToInteger(factor_it, factor_end, text.width());
//         colex_factor = SymbolsToInteger(std::prev(factor_end), std::prev(factor_it), text.width(), -1);
//         if (lex_factor_counts.find(lex_factor) == lex_factor_counts.end())
//         {
//           lex_factor_counts[lex_factor] = 0;
//         }
//         ++lex_factor_counts[lex_factor];
//         colex_to_lex_factor[colex_factor] = lex_factor;
//         factor_it = factor_end;
//       }
//       if (std::distance(factor_it, sl_factor_end) != 0)
//       {
//         lex_factor = SymbolsToInteger(factor_it, sl_factor_end, text.width());
//         colex_factor = SymbolsToInteger(std::prev(sl_factor_end), std::prev(factor_it), text.width(), -1);
//         if (lex_factor_counts.find(lex_factor) == lex_factor_counts.end())
//         {
//           lex_factor_counts[lex_factor] = 0;
//         }
//         ++lex_factor_counts[lex_factor];
//         colex_to_lex_factor[colex_factor] = lex_factor;
//       }
//       sl_factor_end = std::next(text_it);
//     }
//     next_symbol = *text_it--;
//     next_sl_type = sl_type;
//   }
//   lex_factor = SymbolsToInteger(std::next(text_prev_begin), sl_factor_end, text.width());
//   colex_factor = SymbolsToInteger(std::prev(sl_factor_end), text_prev_begin, text.width(), -1);
//   if (lex_factor_counts.find(lex_factor) == lex_factor_counts.end())
//   {
//     lex_factor_counts[lex_factor] = 0;
//   }
//   ++lex_factor_counts[lex_factor];
//   colex_to_lex_factor[colex_factor] = lex_factor;
//   return;
// }

// template <typename Index, typename Node = InformationNode<std::string, uint64_t>>
// uint64_t SerializeIndex
// (
//   Index const &index,
//   std::filesystem::path const &index_path,
//   std::shared_ptr<Node> root = nullptr
// )
// {
//   std::fstream index_file(index_path, std::ios_base::out | std::ios_base::trunc);
//   std::cout << "serialize index to " << std::filesystem::canonical(index_path) << "\n";
//   if (root == nullptr)
//   {
//     sdsl::serialize(index.grammar_rules, index_file);
//     SerializeStaticGrammarTrie(index.lex_grammar_count_trie, index_file);
//     SerializeStaticGrammarTrie(index.lex_grammar_rank_trie, index_file);
//     SerializeStaticGrammarTrie(index.colex_grammar_rank_trie, index_file);
//     sdsl::serialize(index.colex_to_lex, index_file);
//     sdsl::serialize(index.lex_rank_bucket_begin_offsets, index_file);
//     sdsl::serialize(index.colex_bwt, index_file);
//   }
//   else
//   {
//     {
//       auto node {std::make_shared<Node>("grammar_rules")};
//       node->value = sdsl::serialize(index.grammar_rules, index_file);
//       root->value += node->value;
//       root->children.emplace_back(node);
//     }
//     {
//       auto node {std::make_shared<Node>("lex_grammar_count_trie")};
//       SerializeStaticGrammarTrie(index.lex_grammar_count_trie, index_file, node);
//       root->value += node->value;
//       root->children.emplace_back(node);
//     }
//     {
//       auto node {std::make_shared<Node>("lex_grammar_rank_trie")};
//       SerializeStaticGrammarTrie(index.lex_grammar_rank_trie, index_file, node);
//       root->value += node->value;
//       root->children.emplace_back(node);
//     }
//     {
//       auto node {std::make_shared<Node>("colex_grammar_rank_trie")};
//       SerializeStaticGrammarTrie(index.colex_grammar_rank_trie, index_file, node);
//       root->value += node->value;
//       root->children.emplace_back(node);
//     }
//     {
//       auto node {std::make_shared<Node>("colex_to_lex")};
//       node->value = sdsl::serialize(index.colex_to_lex, index_file);
//       root->value += node->value;
//       root->children.emplace_back(node);
//     }
//     {
//       auto node {std::make_shared<Node>("lex_rank_bucket_begin_offsets")};
//       node->value = sdsl::serialize(index.lex_rank_bucket_begin_offsets, index_file);
//       root->value += node->value;
//       root->children.emplace_back(node);
//     }
//     {
//       auto node {std::make_shared<Node>("colex_bwt")};
//       node->value = sdsl::serialize(index.colex_bwt, index_file);
//       root->value += node->value;
//       root->children.emplace_back(node);
//     }
//     return root->value;
//   }
//   return 0;
// }
//
// template <typename Index>
// void LoadIndex (Index &index, std::filesystem::path const &index_path)
// {
//   std::ifstream index_file {index_path};
//   std::cout << "load index from " << std::filesystem::canonical(index_path) << "\n";
//   index.grammar_rules.load(index_file);
//   LoadStaticGrammarTrie(index.lex_grammar_count_trie, index_file);
//   SetLabels(index.grammar_rules, index.lex_grammar_count_trie);
//   LoadStaticGrammarTrie(index.lex_grammar_rank_trie, index_file);
//   SetLabels(index.grammar_rules, index.lex_grammar_rank_trie);
//   LoadStaticGrammarTrie(index.colex_grammar_rank_trie, index_file);
//   SetLabels(index.grammar_rules, index.colex_grammar_rank_trie);
//   index.colex_to_lex.load(index_file);
//   index.lex_rank_bucket_begin_offsets.load(index_file);
//   index.colex_bwt.load(index_file);
//   return;
// }
//
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
