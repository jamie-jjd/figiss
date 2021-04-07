#pragma once

#include <map>
#include <memory>
#include <deque>

#include "utility.h"

namespace project
{
constexpr uint64_t S {1};
constexpr uint64_t L {0};

template
<
  typename ByteText,
  typename CondensedText,
  typename ByteToCondensedAlphabet
>
void ConvertByteToCondensedText
(
  ByteText const &byte_text,
  CondensedText &condensed_text,
  ByteToCondensedAlphabet &byte_to_condensed_alphabet
)
{
  byte_to_condensed_alphabet.resize(256);
  sdsl::util::set_to_value(byte_to_condensed_alphabet, 0);
  for (auto const &byte : byte_text)
  {
    byte_to_condensed_alphabet[byte] = 1;
  }
  byte_to_condensed_alphabet[0] = 0;
  std::partial_sum
  (
    std::begin(byte_to_condensed_alphabet),
    std::end(byte_to_condensed_alphabet),
    std::begin(byte_to_condensed_alphabet)
  );
  uint64_t text_size {std::size(byte_text)};
  if (*std::prev(std::end(byte_text)) != 0)
  {
    ++text_size;
  }
  auto max_condensed_character {*std::prev(std::end(byte_to_condensed_alphabet))};
  condensed_text.width(sdsl::bits::hi(max_condensed_character) + 1);
  condensed_text.resize(text_size);
  auto condensed_text_iterator {std::begin(condensed_text)};
  auto byte_text_iterator {std::begin(byte_text)};
  while (condensed_text_iterator != std::end(condensed_text))
  {
    *condensed_text_iterator++ = byte_to_condensed_alphabet[*byte_text_iterator++];
  }
  return;
}

template <typename TextIterator>
void CalculateNextSlFactor
(
  TextIterator reverse_text_end,
  TextIterator &reverse_first_iterator,
  TextIterator &reverse_last_iterator
)
{
  uint64_t previous_sl_type {L};
  reverse_first_iterator = reverse_last_iterator--;
  while
  (
    (reverse_last_iterator != reverse_text_end)
    &&
    ! (
        (previous_sl_type == S)
        &&
        (*reverse_last_iterator > *std::next(reverse_last_iterator))
      )
  )
  {
    if
    (
      (previous_sl_type == L)
      &&
      (*reverse_last_iterator < *std::next(reverse_last_iterator))
    )
    {
      previous_sl_type = S;
    }
    --reverse_last_iterator;
  }
  return;
}

struct SlFactorTrieNode
{
  std::map<uint64_t, std::shared_ptr<SlFactorTrieNode>> branches;
  std::pair<uint64_t, uint64_t> offset_range;
  std::pair<uint64_t, uint64_t> rank_range;
  uint64_t count;

  SlFactorTrieNode() = default;
  SlFactorTrieNode
  (
    std::pair<uint64_t, uint64_t> const offset_range_,
    std::pair<uint64_t, uint64_t> const rank_range_ = {0, 0},
    uint64_t const count_ = 0
  )
  : offset_range {offset_range_},
    rank_range {rank_range_},
    count {count_}
  {
  }
};

template
<
  typename File,
  typename TextIterator
>
void PrintSlFactorTrie
(
  File &file,
  std::shared_ptr<SlFactorTrieNode> root,
  TextIterator text_begin,
  int64_t const offset = 1
)
{
  std::deque<std::pair<std::shared_ptr<SlFactorTrieNode>, uint64_t>> nodes_depth;
  nodes_depth.emplace_back(root, 0);
  while (!nodes_depth.empty())
  {
    auto current_node {std::get<0>(nodes_depth.back())};
    auto depth {std::get<1>(nodes_depth.back())};
    nodes_depth.pop_back();
    file << depth << ":";
    auto edge_iterator {std::next(text_begin, std::get<0>(current_node->offset_range))};
    auto edge_end {std::next(text_begin, std::get<1>(current_node->offset_range))};
    while (edge_iterator != edge_end)
    {
      file << *edge_iterator;
      edge_iterator += offset;
    }
    file
    << "(" << std::get<0>(current_node->rank_range) << "," << std::get<1>(current_node->rank_range) << ")"
    << "(" << current_node->count << ")"
    << "\n";
    if (!current_node->branches.empty())
    {
      auto reverse_branches_iterator {std::rbegin(current_node->branches)};
      auto reverse_branches_end {std::rend(current_node->branches)};
      while (reverse_branches_iterator != reverse_branches_end)
      {
        nodes_depth.emplace_back(std::get<1>(*reverse_branches_iterator), (depth + 1));
        ++reverse_branches_iterator;
      }
    }
  }
  return;
}

template <typename TextIterator>
void InsertSlFactor
(
  std::shared_ptr<SlFactorTrieNode> current_node,
  TextIterator text_begin,
  TextIterator sl_factor_iterator,
  TextIterator sl_factor_end,
  int64_t const offset = 1
)
{
  if (sl_factor_iterator != sl_factor_end)
  {
    while (true)
    {
      auto character {*sl_factor_iterator};
      auto branches_iterator {current_node->branches.find(character)};
      if (branches_iterator == std::end(current_node->branches))
      {
        current_node->branches[character] = std::make_shared<SlFactorTrieNode>
        (
          std::pair<uint64_t, uint64_t>
          {
            std::distance(text_begin, sl_factor_iterator),
            std::distance(text_begin, sl_factor_end)
          },
          std::pair<uint64_t, uint64_t>{1, 0},
          1
        );
        return;
      }
      else
      {
        auto child_node {std::get<1>(*branches_iterator)};
        auto edge_iterator {std::next(text_begin, std::get<0>(child_node->offset_range))};
        auto edge_end {std::next(text_begin, std::get<1>(child_node->offset_range))};
        while
        (
          (sl_factor_iterator != sl_factor_end)
          &&
          (edge_iterator != edge_end)
          &&
          (*sl_factor_iterator == *edge_iterator)
        )
        {
          sl_factor_iterator += offset;
          edge_iterator += offset;
        }
        if (edge_iterator == edge_end)
        {
          if (sl_factor_iterator != sl_factor_end)
          {
            current_node = child_node;
          }
          else
          {
            std::get<0>(child_node->rank_range) = 1;
            ++(child_node->count);
            return;
          }
        }
        else
        {
          auto edge_iterator_offset {std::distance(text_begin, edge_iterator)};
          auto internal_node
          {
            std::make_shared<SlFactorTrieNode>
            (
              std::pair<uint64_t, uint64_t>
              {
                std::get<0>(child_node->offset_range),
                edge_iterator_offset
              }
            )
          };
          std::get<1>(*branches_iterator) = internal_node;
          internal_node->branches[*edge_iterator] = child_node;
          std::get<0>(child_node->offset_range) = edge_iterator_offset;
          if (sl_factor_iterator != sl_factor_end)
          {
            current_node = internal_node;
          }
          else
          {
            std::get<0>(internal_node->rank_range) = 1;
            ++(internal_node->count);
            return;
          }
        }
      }
    }
  }
  return;
}

template <typename Text>
void InsertSlFactors
(
  std::shared_ptr<SlFactorTrieNode> lex_sl_factor_trie_root,
  std::shared_ptr<SlFactorTrieNode> colex_sl_factor_trie_root,
  Text const &text,
  uint64_t &grammar_compressed_text_size
)
{
  auto reverse_text_begin {std::prev(std::end(text), 2)};
  auto reverse_text_end {std::prev(std::begin(text))};
  auto reverse_sl_factor_begin {reverse_text_begin};
  auto reverse_sl_factor_end {reverse_text_begin};
  while (reverse_sl_factor_end != reverse_text_end)
  {
    ++grammar_compressed_text_size;
    CalculateNextSlFactor
    (
      reverse_text_end,
      reverse_sl_factor_begin,
      reverse_sl_factor_end
    );
    // Print
    // (
    //   std::cout,
    //   std::next(reverse_sl_factor_end),
    //   std::next(reverse_sl_factor_begin)
    // );
    InsertSlFactor
    (
      lex_sl_factor_trie_root,
      std::begin(text),
      std::next(reverse_sl_factor_end),
      std::next(reverse_sl_factor_begin)
    );
    // PrintSlFactorTrie
    // (
    //   std::cout,
    //   lex_sl_factor_trie_root,
    //   std::begin(text)
    // );
    // Print
    // (
    //   std::cout,
    //   reverse_sl_factor_begin,
    //   reverse_sl_factor_end,
    //   -1
    // );
    InsertSlFactor
    (
      colex_sl_factor_trie_root,
      std::begin(text),
      reverse_sl_factor_begin,
      reverse_sl_factor_end,
      -1
    );
    // PrintSlFactorTrie
    // (
    //   std::cout,
    //   colex_sl_factor_trie_root,
    //   std::begin(text),
    //   -1
    // );
  }
}

template <typename TextIterator>
auto LookupSlFactorRankRange
(
  std::shared_ptr<SlFactorTrieNode> current_node,
  TextIterator text_begin,
  TextIterator sl_factor_iterator,
  TextIterator sl_factor_end,
  int64_t const offset = 1
)
{
  if
  (
    (current_node != nullptr)
    &&
    (sl_factor_iterator != sl_factor_end)
  )
  {
    auto branches_iterator {std::end(current_node->branches)};
    while ((branches_iterator = current_node->branches.find(*sl_factor_iterator)) != std::end(current_node->branches))
    {
      current_node = std::get<1>(*branches_iterator);
      auto edge_iterator {std::next(text_begin, std::get<0>(current_node->offset_range))};
      auto edge_end {std::next(text_begin, std::get<1>(current_node->offset_range))};
      while
      (
        (sl_factor_iterator != sl_factor_end)
        &&
        (edge_iterator != edge_end)
        &&
        (*sl_factor_iterator == *edge_iterator)
      )
      {
        sl_factor_iterator += offset;
        edge_iterator += offset;
      }
      if (sl_factor_iterator != sl_factor_end)
      {
        if (edge_iterator != edge_end)
        {
          return decltype(current_node->rank_range){};
        }
      }
      else
      {
        return current_node->rank_range;
      }
    }
  }
  return decltype(current_node->rank_range){};
}

void CalculateSlFactorTrieRankRanges (std::shared_ptr<SlFactorTrieNode> root)
{
  uint64_t rank {};
  std::deque<std::pair<std::shared_ptr<SlFactorTrieNode>, bool>> nodes_is_forward;
  nodes_is_forward.emplace_back(root, true);
  while (!nodes_is_forward.empty())
  {
    auto current_node {std::get<0>(nodes_is_forward.back())};
    auto &is_forward {std::get<1>(nodes_is_forward.back())};
    if (is_forward)
    {
      is_forward = false;
      if (std::get<0>(current_node->rank_range) != 0)
      {
        ++rank;
        current_node->rank_range = {rank, rank};
      }
      if (!current_node->branches.empty())
      {
        auto reverse_branches_iterator {std::rbegin(current_node->branches)};
        auto reverse_branches_end {std::rend(current_node->branches)};
        while (reverse_branches_iterator != reverse_branches_end)
        {
          nodes_is_forward.emplace_back(std::get<1>(*reverse_branches_iterator), true);
          ++reverse_branches_iterator;
        }
      }
    }
    else
    {
      if (!current_node->branches.empty())
      {
        auto first_child_node {std::get<1>(*std::begin(current_node->branches))};
        auto last_child_node {std::get<1>(*std::prev(std::end(current_node->branches)))};
        if (std::get<0>(current_node->rank_range) == 0)
        {
          std::get<0>(current_node->rank_range) = std::get<0>(first_child_node->rank_range);
        }
        std::get<1>(current_node->rank_range) = std::get<1>(last_child_node->rank_range);
      }
      nodes_is_forward.pop_back();
    }
  }
  return;
}

template
<
  typename Text,
  typename GrammarCompressedText
>
void CalculateGrammarCompressedText
(
  Text const &text,
  std::shared_ptr<SlFactorTrieNode> lex_sl_factor_trie_root,
  GrammarCompressedText &grammar_compressed_text
)
{
  auto reverse_text_begin {std::prev(std::end(text), 2)};
  auto reverse_text_end {std::prev(std::begin(text))};
  auto reverse_grammar_compressed_text_iterator {std::prev(std::end(grammar_compressed_text), 2)};
  auto reverse_sl_factor_begin {reverse_text_begin};
  auto reverse_sl_factor_end {reverse_text_begin};
  while (reverse_sl_factor_end != reverse_text_end)
  {
    CalculateNextSlFactor
    (
      reverse_text_end,
      reverse_sl_factor_begin,
      reverse_sl_factor_end
    );
    // Print
    // (
    //   std::cout,
    //   std::next(reverse_sl_factor_end),
    //   std::next(reverse_sl_factor_begin)
    // );
    auto rank_range
    {
      LookupSlFactorRankRange
      (
        lex_sl_factor_trie_root,
        std::begin(text),
        std::next(reverse_sl_factor_end),
        std::next(reverse_sl_factor_begin)
      )
    };
    *reverse_grammar_compressed_text_iterator-- = std::get<0>(rank_range);
  }
  return;
}

void CalculateSlFactorTrieCounts (std::shared_ptr<SlFactorTrieNode> root)
{
  std::deque<std::pair<std::shared_ptr<SlFactorTrieNode>, bool>> nodes_is_forward;
  nodes_is_forward.emplace_back(root, true);
  while (!nodes_is_forward.empty())
  {
    auto current_node {std::get<0>(nodes_is_forward.back())};
    auto &is_forward {std::get<1>(nodes_is_forward.back())};
    if (is_forward)
    {
      is_forward = false;
      if (!current_node->branches.empty())
      {
        auto branches_iterator {std::begin(current_node->branches)};
        auto branches_end {std::end(current_node->branches)};
        while (branches_iterator != branches_end)
        {
          nodes_is_forward.emplace_back(std::get<1>(*branches_iterator), true);
          ++branches_iterator;
        }
      }
    }
    else
    {
      if (!current_node->branches.empty())
      {
        auto branches_iterator {std::begin(current_node->branches)};
        auto branches_end {std::end(current_node->branches)};
        while (branches_iterator != branches_end)
        {
          auto child_node {std::get<1>(*branches_iterator)};
          current_node->count += child_node->count;
          ++branches_iterator;
        }
      }
      nodes_is_forward.pop_back();
    }
  }
  return;
}
// template <typename Text>
// uint64_t CalculateMaxSlFactorSize (Text const &text)
// {
//   uint64_t max_size {0};
//   auto reverse_first_iterator {std::prev(std::end(text))};
//   auto reverse_last_iterator {std::prev(std::end(text))};
//   auto reverse_end {std::prev(std::begin(text))};
//   while (reverse_last_iterator != reverse_end)
//   {
//     CalculateSlFactor(reverse_first_iterator, reverse_last_iterator, reverse_end);
//     auto size {static_cast<uint64_t>(std::distance(reverse_last_iterator, reverse_first_iterator))};
//     if (max_size < size)
//     {
//       max_size = size;
//     }
//   }
//   return max_size;
// }

// template
// <
//   typename Text,
//   typename CharacterBucketOffsets
// >
// void CalculateCharacterBucketEndOffsets
// (
//   Text const &text,
//   CharacterBucketOffsets &character_bucket_offsets
// )
// {
//   auto begin {std::begin(character_bucket_offsets)};
//   auto end {std::end(character_bucket_offsets)};
//   sdsl::util::set_to_value(character_bucket_offsets, 0);
//   for (auto const &character : text)
//   {
//     ++character_bucket_offsets[character];
//   }
//   std::partial_sum(begin, end, begin);
//   return;
// }
//
// template
// <
//   typename Text,
//   typename CharacterBucketOffsets
// >
// void CalculateCharacterBucketBeginOffsets
// (
//   Text const &text,
//   CharacterBucketOffsets &character_bucket_offsets
// )
// {
//   CalculateCharacterBucketEndOffsets(text, character_bucket_offsets);
//   auto iterator {std::prev(std::end(character_bucket_offsets))};
//   auto last_iterator {std::begin(character_bucket_offsets)};
//   while (iterator != last_iterator)
//   {
//     *iterator = *std::prev(iterator);
//     --iterator;
//   }
//   *last_iterator = 0;
//   return;
// }

// template <typename LexToColexOrderMapping>
// void CalculateColexTrieRankRangesAndLexToColexOrderMapping
// (
//   TrieNode *current_node,
//   LexToColexOrderMapping &lex_to_colex_order_mapping,
//   uint64_t &colex_rank
// )
// {
//   if (current_node->leftmost_rank != 0)
//   {
//     lex_to_colex_order_mapping[current_node->leftmost_rank] = colex_rank;
//     current_node->leftmost_rank = current_node->rightmost_rank = colex_rank;
//     ++colex_rank;
//   }
//   if (!current_node->branches.empty())
//   {
//     auto branches_begin {std::begin(current_node->branches)};
//     auto branches_iterator {branches_begin};
//     auto branches_end {std::end(current_node->branches)};
//     while (branches_iterator != branches_end)
//     {
//       CalculateColexTrieRankRangesAndLexToColexOrderMapping
//       (
//         std::get<1>(*branches_iterator),
//         lex_to_colex_order_mapping,
//         colex_rank
//       );
//       ++branches_iterator;
//     }
//     auto first_child_node {std::get<1>(*branches_begin)};
//     auto last_child_node {std::get<1>(*std::prev(branches_end))};
//     if (current_node->leftmost_rank == 0)
//     {
//       current_node->leftmost_rank = first_child_node->leftmost_rank;
//     }
//     current_node->rightmost_rank = last_child_node->rightmost_rank;
//   }
//   return;
// }
//
// template
// <
//   typename LexToColexOrderMapping,
//   typename ColexToLexOrderMapping
// >
// void CalculateColexToLexOrderMapping
// (
//   LexToColexOrderMapping const &lex_to_colex_order_mapping,
//   ColexToLexOrderMapping &colex_to_lex_order_mapping
// )
// {
//   for (uint64_t lex_rank {0}; lex_rank != std::size(lex_to_colex_order_mapping); ++lex_rank)
//   {
//     colex_to_lex_order_mapping[lex_to_colex_order_mapping[lex_rank]] = lex_rank;
//   }
//   return;
// }
//
// template
// <
//   typename GrammarCompressedText,
//   typename LexToColexOrderMapping,
//   typename ColexGrammarCompressedBwt
// >
// void CalculateColexGrammarCompressedBwt
// (
//   GrammarCompressedText const &grammar_compressed_text,
//   LexToColexOrderMapping const &lex_to_colex_order_mapping,
//   ColexGrammarCompressedBwt &colex_grammar_compressed_bwt
// )
// {
//   sdsl::int_vector<> buffer;
//   sdsl::qsufsort::construct_sa(buffer, grammar_compressed_text);
//   auto buffer_iterator {std::begin(buffer)};
//   auto buffer_end {std::end(buffer)};
//   while (buffer_iterator != buffer_end)
//   {
//     if (*buffer_iterator != 0)
//     {
//       *buffer_iterator = lex_to_colex_order_mapping[grammar_compressed_text[(*buffer_iterator - 1)]];
//     }
//     ++buffer_iterator;
//   }
//   sdsl::construct_im(colex_grammar_compressed_bwt, buffer);
//   return;
// }

struct Index
{
  sdsl::int_vector<8> byte_to_condensed_alphabet;
};

template <typename File>
void PrintIndex
(
  File &file,
  Index &index
)
{
  Print(file, std::begin(index.byte_to_condensed_alphabet), std::end(index.byte_to_condensed_alphabet));
  return;
}

void Construct
(
  Index &index,
  std::filesystem::path const &byte_text_path
)
{
  sdsl::int_vector<> text;
  {
    sdsl::int_vector<8> byte_text;
    sdsl::load_vector_from_file(byte_text, byte_text_path);
    ConvertByteToCondensedText(byte_text, text, index.byte_to_condensed_alphabet);
  }
  auto lex_sl_factor_trie_root {std::make_shared<SlFactorTrieNode>()};
  auto colex_sl_factor_trie_root {std::make_shared<SlFactorTrieNode>()};
  uint64_t grammar_compressed_text_size {1};
  InsertSlFactors
  (
    lex_sl_factor_trie_root,
    colex_sl_factor_trie_root,
    text,
    grammar_compressed_text_size
  );
  CalculateSlFactorTrieRankRanges(lex_sl_factor_trie_root);
  // PrintSlFactorTrie
  // (
  //   std::cout,
  //   lex_sl_factor_trie_root,
  //   std::begin(text)
  // );
  CalculateSlFactorTrieRankRanges(colex_sl_factor_trie_root);
  // PrintSlFactorTrie
  // (
  //   std::cout,
  //   colex_sl_factor_trie_root,
  //   std::begin(text),
  //   -1
  // );
  auto grammar_compressed_alphabet_size {std::get<1>(lex_sl_factor_trie_root->rank_range) + 1};
  auto grammar_compressed_text_width {sdsl::bits::hi(grammar_compressed_alphabet_size - 1) + 1};
  sdsl::int_vector<> grammar_compressed_text
  (
    grammar_compressed_text_size,
    0,
    grammar_compressed_text_width
  );
  CalculateGrammarCompressedText
  (
    text,
    lex_sl_factor_trie_root,
    grammar_compressed_text
  );
  // auto grammar_compressed_alphabet_size {std::size(index.grammar_rule_sizes) + 1};
  // auto grammar_compressed_character_width {sdsl::bits::hi(std::size(index.grammar_rule_sizes)) + 1};
  // sdsl::int_vector<> lex_to_colex_order_mapping;
  // lex_to_colex_order_mapping.width(grammar_compressed_character_width);
  // lex_to_colex_order_mapping.resize(grammar_compressed_alphabet_size);
  // lex_to_colex_order_mapping[0] = 0;
  // {
  //   uint64_t colex_rank {1};
  //   CalculateColexTrieRankRangesAndLexToColexOrderMapping
  //   (
  //     index.colex_grammar_trie_root,
  //     lex_to_colex_order_mapping,
  //     colex_rank
  //   );
  // }
  // index.colex_to_lex_order_mapping.width(grammar_compressed_character_width);
  // index.colex_to_lex_order_mapping.resize(grammar_compressed_alphabet_size);
  // CalculateColexToLexOrderMapping
  // (
  //   lex_to_colex_order_mapping,
  //   index.colex_to_lex_order_mapping
  // );
  // auto grammar_compressed_text_size
  // {
  //   std::distance
  //   (
  //     temporary_grammar_compressed_text_begin,
  //     temporary_grammar_compressed_text_end
  //   ) + 1
  // };
  // sdsl::int_vector<> grammar_compressed_text;
  // grammar_compressed_text.width(grammar_compressed_character_width);
  // grammar_compressed_text.resize(grammar_compressed_text_size);
  // *std::prev(std::end(grammar_compressed_text)) = 0;
  // index.lex_grammar_compressed_character_bucket_end_offsets.width(sdsl::bits::hi(grammar_compressed_text_size) + 1);
  // index.lex_grammar_compressed_character_bucket_end_offsets.resize(grammar_compressed_alphabet_size);
  // CalculateCharacterBucketEndOffsets
  // (
  //   grammar_compressed_text,
  //   index.lex_grammar_compressed_character_bucket_end_offsets
  // );
  // CalculateColexGrammarCompressedBwt
  // (
  //   grammar_compressed_text,
  //   lex_to_colex_order_mapping,
  //   index.colex_grammar_compressed_bwt
  // );
  return;
}

void Serialize
(
  Index &index,
  std::filesystem::path index_path
)
{
  std::ofstream index_file {index_path};
  sdsl::serialize(index.byte_to_condensed_alphabet, index_file);
  return;
}

void Load
(
  Index &index,
  std::filesystem::path const &index_path
)
{
  std::ifstream index_file {index_path};
  index.byte_to_condensed_alphabet.load(index_file);
  return;
}

// template <typename PatternIterator>
// void BackwardSearchPatternPrefix
// (
//   Index const &index,
//   uint64_t &pattern_range_begin_offset_L,
//   uint64_t &pattern_range_end_offset_L,
//   uint64_t &pattern_range_begin_offset_S,
//   uint64_t &pattern_range_end_offset_S,
//   PatternIterator reverse_pattern_first_iterator,
//   PatternIterator reverse_pattern_last_iterator
// )
// {
//   uint64_t leftmost_colex_rank {0};
//   uint64_t rightmost_colex_rank {0};
//   LookupGrammarRule
//   (
//     std::begin(index.text),
//     index.colex_grammar_trie_root,
//     reverse_pattern_first_iterator,
//     reverse_pattern_last_iterator,
//     -1,
//     leftmost_colex_rank,
//     rightmost_colex_rank
//   );
//   if (leftmost_colex_rank != 0)
//   {
//     if (pattern_range_begin_offset_L != pattern_range_end_offset_L)
//     {
//       pattern_range_end_offset_L = std::get<0>
//       (
//         index.colex_grammar_compressed_bwt.range_search_2d
//         (
//           pattern_range_begin_offset_L,
//           (pattern_range_end_offset_L - 1),
//           leftmost_colex_rank,
//           rightmost_colex_rank,
//           false
//         )
//       );
//       pattern_range_begin_offset_L = 0;
//     }
//     if (pattern_range_begin_offset_S != pattern_range_end_offset_S)
//     {
//       pattern_range_end_offset_S = std::get<0>
//       (
//         index.colex_grammar_compressed_bwt.range_search_2d
//         (
//           pattern_range_begin_offset_S,
//           (pattern_range_end_offset_S - 1),
//           leftmost_colex_rank,
//           rightmost_colex_rank,
//           false
//         )
//       );
//       pattern_range_begin_offset_S = 0;
//     }
//   }
//   // {
//   //   auto iterator {std::next(reverse_pattern_last_iterator)};
//   //   auto end {std::next(reverse_pattern_first_iterator)};
//   //   while ( iterator != end)
//   //   {
//   //     std::cout << *iterator;
//   //     ++iterator;
//   //   }
//   //   std::cout
//   //   << "->(" << leftmost_colex_rank << "," << rightmost_colex_rank << ")"
//   //   << "->L:(" << pattern_range_begin_offset_L << "," << pattern_range_end_offset_L << ")"
//   //   << "->S:(" << pattern_range_begin_offset_S << "," << pattern_range_end_offset_S << ")\n";
//   // }
//   return;
// }
//
// template <typename PatternIterator>
// void BackwardSearchExactSlFactor
// (
//   Index const &index,
//   uint64_t &pattern_range_begin_offset_L,
//   uint64_t &pattern_range_end_offset_L,
//   uint64_t &pattern_range_begin_offset_S,
//   uint64_t &pattern_range_end_offset_S,
//   PatternIterator reverse_pattern_first_iterator,
//   PatternIterator reverse_pattern_last_iterator
// )
// {
//   uint64_t colex_rank {0};
//   LookupGrammarRule
//   (
//     std::begin(index.text),
//     index.colex_grammar_trie_root,
//     reverse_pattern_first_iterator,
//     reverse_pattern_last_iterator,
//     -1,
//     colex_rank,
//     colex_rank
//   );
//   if (colex_rank != 0)
//   {
//     auto character_begin_offset
//     {
//       index.lex_grammar_compressed_character_bucket_end_offsets
//       [
//         index.colex_to_lex_order_mapping[colex_rank] - 1
//       ]
//     };
//     if (pattern_range_begin_offset_L != pattern_range_end_offset_L)
//     {
//       pattern_range_begin_offset_L =
//       (
//         character_begin_offset
//         + index.colex_grammar_compressed_bwt.rank
//         (
//           pattern_range_begin_offset_L,
//           colex_rank
//         )
//       );
//       pattern_range_end_offset_L =
//       (
//         character_begin_offset
//         + index.colex_grammar_compressed_bwt.rank
//         (
//           pattern_range_end_offset_L,
//           colex_rank
//         )
//       );
//     }
//     if (pattern_range_begin_offset_S != pattern_range_end_offset_S)
//     {
//       pattern_range_begin_offset_S =
//       (
//         character_begin_offset
//         + index.colex_grammar_compressed_bwt.rank
//         (
//           pattern_range_begin_offset_S,
//           colex_rank
//         )
//       );
//       pattern_range_end_offset_S =
//       (
//         character_begin_offset
//         + index.colex_grammar_compressed_bwt.rank
//         (
//           pattern_range_end_offset_S,
//           colex_rank
//         )
//       );
//     }
//   }
//   else
//   {
//     pattern_range_begin_offset_L = pattern_range_end_offset_L;
//     pattern_range_begin_offset_S = pattern_range_end_offset_S;
//   }
//   // {
//   //   auto iterator {std::next(reverse_pattern_last_iterator)};
//   //   auto end {std::next(reverse_pattern_first_iterator)};
//   //   while (iterator != end)
//   //   {
//   //     std::cout << *iterator;
//   //     ++iterator;
//   //   }
//   //   std::cout
//   //   << "->(" << colex_rank << ":" << index.colex_to_lex_order_mapping[colex_rank] << ")"
//   //   << "->L:(" << pattern_range_begin_offset_L << "," << pattern_range_end_offset_L << ")"
//   //   << "->S:(" << pattern_range_begin_offset_S << "," << pattern_range_end_offset_S << ")\n";
//   // }
//   return;
// }
//
// template <typename PatternIterator>
// auto CalculateReversePatternLastIteratorOfPatternSuffixS
// (
//   PatternIterator reverse_begin,
//   PatternIterator reverse_end
// )
// {
//   auto iterator {reverse_begin};
//   while
//   (
//     (std::prev(iterator) != reverse_end)
//     &&
//     (*std::prev(iterator) == *iterator)
//   )
//   {
//     --iterator;
//   }
//   if
//   (
//     (std::prev(iterator) != reverse_end)
//     &&
//     (*std::prev(iterator) < *iterator)
//   )
//   {
//     return reverse_begin;
//   }
//   return std::prev(iterator);
// }
//
// template <typename PatternIterator>
// auto BackwardSearchPatternSuffix
// (
//   Index const &index,
//   uint64_t &pattern_range_begin_offset_L,
//   uint64_t &pattern_range_end_offset_L,
//   uint64_t &pattern_range_begin_offset_S,
//   uint64_t &pattern_range_end_offset_S,
//   PatternIterator reverse_pattern_begin,
//   PatternIterator reverse_pattern_end
// )
// {
//   auto reverse_pattern_first_iterator {reverse_pattern_begin};
//   auto reverse_pattern_last_iterator
//   {
//     CalculateReversePatternLastIteratorOfPatternSuffixS
//     (
//       reverse_pattern_begin,
//       reverse_pattern_end
//     )
//   };
//   uint64_t leftmost_lex_rank {0};
//   uint64_t rightmost_lex_rank {0};
//   if
//   (
//     (reverse_pattern_last_iterator != reverse_pattern_begin)
//     &&
//     (reverse_pattern_last_iterator != reverse_pattern_end)
//   )
//   {
//     LookupGrammarRule
//     (
//       std::begin(index.text),
//       index.lex_grammar_trie_root,
//       std::next(reverse_pattern_last_iterator),
//       std::next(reverse_pattern_begin),
//       1,
//       leftmost_lex_rank,
//       rightmost_lex_rank
//     );
//     if (leftmost_lex_rank != 0)
//     {
//       pattern_range_begin_offset_S = index.lex_grammar_compressed_character_bucket_end_offsets[leftmost_lex_rank - 1];
//       pattern_range_end_offset_S = index.lex_grammar_compressed_character_bucket_end_offsets[rightmost_lex_rank];
//     }
//     // {
//     //   auto iterator {std::next(reverse_pattern_last_iterator)};
//     //   auto end {std::next(reverse_pattern_begin)};
//     //   while (iterator != end)
//     //   {
//     //     std::cout << *iterator;
//     //     ++iterator;
//     //   }
//     //   std::cout << "->(" << leftmost_lex_rank << "," << rightmost_lex_rank << ")";
//     //   std::cout << "->S:(" << pattern_range_begin_offset_S << "," << pattern_range_end_offset_S << ")\n";
//     // }
//     CalculateSlFactor(reverse_pattern_first_iterator, reverse_pattern_last_iterator, reverse_pattern_end);
//     if (pattern_range_begin_offset_S != pattern_range_end_offset_S)
//     {
//       uint64_t leftmost_colex_rank {0};
//       uint64_t rightmost_colex_rank {0};
//       LookupGrammarRule
//       (
//         std::begin(index.text),
//         index.colex_grammar_trie_root,
//         reverse_pattern_first_iterator,
//         reverse_pattern_last_iterator,
//         -1,
//         leftmost_colex_rank,
//         rightmost_colex_rank
//       );
//       if (leftmost_colex_rank != 0)
//       {
//         if (reverse_pattern_last_iterator != reverse_pattern_end)
//         {
//           auto character_begin_offset
//           {
//             index.lex_grammar_compressed_character_bucket_end_offsets
//             [
//               index.colex_to_lex_order_mapping[leftmost_colex_rank] - 1
//             ]
//           };
//           pattern_range_begin_offset_S =
//           (
//             character_begin_offset
//             + index.colex_grammar_compressed_bwt.rank
//             (
//               pattern_range_begin_offset_S,
//               leftmost_colex_rank
//             )
//           );
//           pattern_range_end_offset_S =
//           (
//             character_begin_offset
//             + index.colex_grammar_compressed_bwt.rank
//             (
//               pattern_range_end_offset_S,
//               leftmost_colex_rank
//             )
//           );
//         }
//         else
//         {
//           pattern_range_end_offset_S = std::get<0>
//           (
//             index.colex_grammar_compressed_bwt.range_search_2d
//             (
//               pattern_range_begin_offset_S,
//               (pattern_range_end_offset_S - 1),
//               leftmost_colex_rank,
//               rightmost_colex_rank,
//               false
//             )
//           );
//           pattern_range_begin_offset_S = 0;
//         }
//       }
//       else
//       {
//         pattern_range_begin_offset_S = pattern_range_end_offset_S;
//       }
//       // {
//       //   auto iterator {std::next(reverse_pattern_last_iterator)};
//       //   auto end {std::next(reverse_pattern_first_iterator)};
//       //   while (iterator != end)
//       //   {
//       //     std::cout << *iterator;
//       //     ++iterator;
//       //   }
//       //   std::cout
//       //   << "->(" << leftmost_colex_rank << ":" << index.colex_to_lex_order_mapping[leftmost_colex_rank]
//       //   << "," << rightmost_colex_rank << ")"
//       //   << "->S:(" << pattern_range_begin_offset_S << "," << pattern_range_end_offset_S << ")\n";
//       // }
//     }
//   }
//   else if (reverse_pattern_last_iterator == reverse_pattern_begin)
//   {
//     CalculateSlFactor
//     (
//       reverse_pattern_first_iterator,
//       reverse_pattern_last_iterator,
//       reverse_pattern_end
//     );
//   }
//   LookupGrammarRule
//   (
//     std::begin(index.text),
//     index.lex_grammar_trie_root,
//     std::next(reverse_pattern_last_iterator),
//     std::next(reverse_pattern_begin),
//     1,
//     leftmost_lex_rank,
//     rightmost_lex_rank
//   );
//   if (leftmost_lex_rank != 0)
//   {
//     pattern_range_begin_offset_L = index.lex_grammar_compressed_character_bucket_end_offsets[leftmost_lex_rank - 1];
//     pattern_range_end_offset_L = index.lex_grammar_compressed_character_bucket_end_offsets[rightmost_lex_rank];
//   }
//   // {
//   //   auto iterator {std::next(reverse_pattern_last_iterator)};
//   //   auto end {std::next(reverse_pattern_begin)};
//   //   while (iterator != end)
//   //   {
//   //     std::cout << *iterator;
//   //     ++iterator;
//   //   }
//   //   std::cout << "->(" << leftmost_lex_rank << "," << rightmost_lex_rank << ")";
//   //   std::cout << "->L:(" << pattern_range_begin_offset_L << "," << pattern_range_end_offset_L << ")\n";
//   // }
//   return reverse_pattern_last_iterator;
// }
//
// template <typename PatternIterator>
// uint64_t Count
// (
//   Index const &index,
//   PatternIterator pattern_begin,
//   PatternIterator pattern_end
// )
// {
//   uint64_t pattern_range_begin_offset_L {0};
//   uint64_t pattern_range_end_offset_L {0};
//   uint64_t pattern_range_begin_offset_S {0};
//   uint64_t pattern_range_end_offset_S {0};
//   auto reverse_pattern_begin {std::prev(pattern_end)};
//   auto reverse_pattern_end {std::prev(pattern_begin)};
//   auto reverse_pattern_first_iterator {reverse_pattern_begin};
//   auto reverse_pattern_last_iterator {reverse_pattern_begin};
//   // {
//   //   for (auto iterator {pattern_begin}; iterator != pattern_end; ++iterator)
//   //   {
//   //     std::cout << *iterator;
//   //   }
//   //   std::cout << "\n";
//   // }
//   reverse_pattern_last_iterator = BackwardSearchPatternSuffix
//   (
//     index,
//     pattern_range_begin_offset_L,
//     pattern_range_end_offset_L,
//     pattern_range_begin_offset_S,
//     pattern_range_end_offset_S,
//     reverse_pattern_begin,
//     reverse_pattern_end
//   );
//   if (reverse_pattern_last_iterator != reverse_pattern_end)
//   {
//     while
//     (
//       (pattern_range_begin_offset_L != pattern_range_end_offset_L)
//       ||
//       (pattern_range_begin_offset_S != pattern_range_end_offset_S)
//     )
//     {
//       CalculateSlFactor(reverse_pattern_first_iterator, reverse_pattern_last_iterator, reverse_pattern_end);
//       if (reverse_pattern_last_iterator!= reverse_pattern_end)
//       {
//         BackwardSearchExactSlFactor
//         (
//           index,
//           pattern_range_begin_offset_L,
//           pattern_range_end_offset_L,
//           pattern_range_begin_offset_S,
//           pattern_range_end_offset_S,
//           reverse_pattern_first_iterator,
//           reverse_pattern_last_iterator
//         );
//       }
//       else
//       {
//         BackwardSearchPatternPrefix
//         (
//           index,
//           pattern_range_begin_offset_L,
//           pattern_range_end_offset_L,
//           pattern_range_begin_offset_S,
//           pattern_range_end_offset_S,
//           reverse_pattern_first_iterator,
//           reverse_pattern_last_iterator
//         );
//         break;
//       }
//     }
//   }
//   return
//   (
//     (pattern_range_end_offset_L - pattern_range_begin_offset_L)
//     +
//     (pattern_range_end_offset_S - pattern_range_begin_offset_S)
//   );
// }
}
