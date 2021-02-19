#pragma once

#include <map>

#include <sdsl/suffix_trees.hpp>
#include "utility.h"

namespace project
{
constexpr uint64_t S {1};
constexpr uint64_t L {0};

template <typename SlTypesIterator>
constexpr bool IsRightmostLType (SlTypesIterator iterator)
{
  return
  (
    (*iterator == L)
    &&
    (*std::next(iterator) == S)
  );
}

template <typename SlTypesIterator>
constexpr bool IsLeftmostSType (SlTypesIterator iterator)
{
  return
  (
    (*iterator == S)
    &&
    (*std::prev(iterator) == L)
  );
}

template
<
  typename Text,
  typename SlTypes
>
void CalculateSlTypes
(
  Text const &text,
  SlTypes &sl_types
)
{
  sl_types.resize(std::size(text));
  sdsl::util::set_to_value(sl_types, S);
  auto text_iterator {std::prev(std::end(text))};
  auto text_last_iterator {std::begin(text)};
  auto sl_types_iterator {std::prev(std::end(sl_types))};
  while (text_iterator != text_last_iterator)
  {
    if
    (
      (*std::prev(text_iterator) > *text_iterator)
      ||
      (
        (*std::prev(text_iterator) == *text_iterator)
        &&
        (*sl_types_iterator == L)
      )
    )
    {
      *std::prev(sl_types_iterator) = L;
    }
    --sl_types_iterator;
    --text_iterator;
  }
  return;
}

template
<
  typename Text,
  typename CharacterBucketOffsets
>
void CalculateCharacterBucketEndOffsets
(
  Text const &text,
  CharacterBucketOffsets &character_bucket_offsets
)
{
  auto begin {std::begin(character_bucket_offsets)};
  auto end {std::end(character_bucket_offsets)};
  sdsl::util::_set_zero_bits(character_bucket_offsets);
  for (auto const &character : text)
  {
    ++character_bucket_offsets[character];
  }
  std::partial_sum(begin, end, begin);
  return;
}

template
<
  typename Text,
  typename CharacterBucketOffsets
>
void CalculateCharacterBucketBeginOffsets
(
  Text const &text,
  CharacterBucketOffsets &character_bucket_offsets
)
{
  CalculateCharacterBucketEndOffsets(text, character_bucket_offsets);
  auto iterator {std::prev(std::end(character_bucket_offsets))};
  auto last_iterator {std::begin(character_bucket_offsets)};
  while (iterator != last_iterator)
  {
    *iterator = *std::prev(iterator);
    --iterator;
  }
  *last_iterator = 0;
  return;
}

template
<
  typename Text,
  typename SlTypes,
  typename CharacterBucketOffsets,
  typename TextOffsets
>
void BucketSortRightmostLTypeCharacters
(
  Text const &text,
  SlTypes const &sl_types,
  CharacterBucketOffsets &character_bucket_offsets,
  TextOffsets &text_offsets
)
{
  auto text_iterator {std::begin(text)};
  auto text_end {std::end(text)};
  while (text_iterator != text_end)
  {
    auto text_offset {std::distance(std::begin(text), text_iterator)};
    if (IsRightmostLType(std::next(std::begin(sl_types), text_offset)))
    {
      text_offsets[character_bucket_offsets[*text_iterator]++] = text_offset;
    }
    ++text_iterator;
  }
  return;
}

template
<
  typename Text,
  typename SlTypes,
  typename CharacterBucketOffsets,
  typename TextOffsets
>
void InduceSortLTypeCharacters
(
  Text const &text,
  SlTypes const &sl_types,
  CharacterBucketOffsets &character_bucket_offsets,
  TextOffsets &text_offsets
)
{
  auto text_offsets_iterator {std::next(std::begin(text_offsets))};
  auto text_offsets_last_iterator {std::end(text_offsets)};
  auto invalid_text_offset {std::size(text)};
  while (text_offsets_iterator != text_offsets_last_iterator)
  {
    auto text_offset {*text_offsets_iterator};
    if
    (
      (text_offset != invalid_text_offset)
      &&
      (text_offset != 0)
      &&
      (sl_types[text_offset - 1] == L)
    )
    {
      text_offsets[character_bucket_offsets[text[text_offset - 1]]++] = (text_offset - 1);
      *text_offsets_iterator = invalid_text_offset;
    }
    ++text_offsets_iterator;
  }
  return;
}

template
<
  typename Text,
  typename SlTypes,
  typename CharacterBucketOffsets,
  typename TextOffsets
>
void InduceSortSTypeCharacters
(
  Text const &text,
  SlTypes const &sl_types,
  CharacterBucketOffsets &character_bucket_offsets,
  TextOffsets &text_offsets
)
{
  auto text_offsets_iterator {std::prev(std::end(text_offsets))};
  auto text_offsets_last_iterator {std::begin(text_offsets)};
  auto invalid_text_offset {std::size(text)};
  while (text_offsets_iterator != text_offsets_last_iterator)
  {
    auto text_offset {*text_offsets_iterator};
    if
    (
      (text_offset != invalid_text_offset)
      &&
      (text_offset != 0)
      &&
      (sl_types[text_offset - 1] == S)
    )
    {
      text_offsets[--character_bucket_offsets[text[text_offset - 1]]] = (text_offset - 1);
      *text_offsets_iterator = invalid_text_offset;
    }
    --text_offsets_iterator;
  }
  return;
}

template <typename RandomAccessIterator>
auto MoveVaildEntriesToFront
(
  RandomAccessIterator first_iterator,
  RandomAccessIterator last_iterator,
  uint64_t const invalid_value
)
{
  auto iterator {first_iterator};
  auto new_last_iterator {first_iterator};
  while (iterator != last_iterator)
  {
    if (*iterator != invalid_value)
    {
      uint64_t temporary_value {*iterator};
      *iterator = *new_last_iterator;
      *new_last_iterator = temporary_value;
      ++new_last_iterator;
    }
    ++iterator;
  }
  return new_last_iterator;
}

template
<
  typename Text,
  typename SlTypes,
  typename GrammarRuleBeginOffsetsIterator,
  typename TemporaryGrammarCompressedTextIterator
>
void CalculateGrammarRuleBeginOffsetsAndTemporaryGrammarCompressedText
(
  Text const &text,
  SlTypes const &sl_types,
  GrammarRuleBeginOffsetsIterator begin_offsets_begin,
  GrammarRuleBeginOffsetsIterator &begin_offsets_end,
  TemporaryGrammarCompressedTextIterator temporary_grammar_compressed_text_begin,
  TemporaryGrammarCompressedTextIterator &temporary_grammar_compressed_text_end
)
{
  auto invalid_text_offset {std::size(text)};
  uint64_t lex_rank {0};
  auto previous_grammar_rule_iterator {std::prev(std::end(text))};
  auto previous_sl_types_iterator {std::prev(std::end(sl_types))};
  auto begin_offsets_iterator {begin_offsets_begin};
  while (begin_offsets_iterator != begin_offsets_end)
  {
    uint64_t begin_offset {*begin_offsets_iterator};
    auto grammar_rule_iterator {std::next(std::begin(text), begin_offset)};
    auto sl_types_iterator {std::next(std::begin(sl_types), begin_offset)};
    if
    (
      (*previous_grammar_rule_iterator == *grammar_rule_iterator)
      &&
      (*previous_sl_types_iterator == *sl_types_iterator)
    )
    {
      do
      {
        ++previous_grammar_rule_iterator;
        ++grammar_rule_iterator;
        ++previous_sl_types_iterator;
        ++sl_types_iterator;
      }
      while
      (
        !IsLeftmostSType(previous_sl_types_iterator)
        &&
        !IsLeftmostSType(sl_types_iterator)
        &&
        (*previous_grammar_rule_iterator == *grammar_rule_iterator)
        &&
        (*previous_sl_types_iterator == *sl_types_iterator)
      );
      if
      (
        IsLeftmostSType(previous_sl_types_iterator)
        &&
        IsLeftmostSType(sl_types_iterator)
      )
      {
        --lex_rank;
        *begin_offsets_iterator = invalid_text_offset;
      }
    }
    *std::next(temporary_grammar_compressed_text_begin, (begin_offset + 1) / 2) = ++lex_rank;
    previous_grammar_rule_iterator = std::next(std::begin(text), begin_offset);
    previous_sl_types_iterator = std::next(std::begin(sl_types), begin_offset);
    ++begin_offsets_iterator;
  }
  begin_offsets_end = MoveVaildEntriesToFront
  (
    begin_offsets_begin,
    begin_offsets_end,
    invalid_text_offset
  );
  temporary_grammar_compressed_text_end = MoveVaildEntriesToFront
  (
    temporary_grammar_compressed_text_begin,
    temporary_grammar_compressed_text_end,
    invalid_text_offset
  );
  return;
}

template
<
  typename GrammarRuleSizes,
  typename SlTypes,
  typename GrammarRuleBeginOffsetsIterator
>
void CalculateGrammarRuleSizes
(
  GrammarRuleSizes &sizes,
  SlTypes const &sl_types,
  GrammarRuleBeginOffsetsIterator begin_offsets_begin,
  GrammarRuleBeginOffsetsIterator begin_offsets_end
)
{
  sizes.resize(std::distance(begin_offsets_begin, begin_offsets_end));
  auto sizes_iterator {std::begin(sizes)};
  auto begin_offsets_iterator {begin_offsets_begin};
  while (begin_offsets_iterator != begin_offsets_end)
  {
    auto sl_types_first_iterator {std::next(std::begin(sl_types), *begin_offsets_iterator)};
    auto sl_types_last_iterator {std::next(sl_types_first_iterator)};
    while (!IsLeftmostSType(sl_types_last_iterator))
    {
      ++sl_types_last_iterator;
    }
    *sizes_iterator = std::distance(sl_types_first_iterator, sl_types_last_iterator);
    ++sizes_iterator;
    ++begin_offsets_iterator;
  }
  return;
}

struct TrieNode
{
  int64_t edge_begin_offset;
  int64_t edge_end_offset;
  std::map<uint64_t, TrieNode*> branches;
  uint64_t leftmost_rank;
  uint64_t rightmost_rank;

  TrieNode () = default;

  TrieNode
  (
    int64_t begin_offset,
    int64_t end_offset,
    uint64_t rank
  )
  : edge_begin_offset {begin_offset},
    edge_end_offset {end_offset},
    leftmost_rank {rank},
    rightmost_rank {rank}
  {
  }

  ~TrieNode ()
  {
    auto branches_iterator {std::begin(branches)};
    auto branches_end {std::end(branches)};
    while (branches_iterator != branches_end)
    {
      delete std::get<1>(*branches_iterator);
      ++branches_iterator;
    }
  }
};

template
<
  typename Text,
  typename GrammarRuleSizes,
  typename GrammarRules,
  typename GrammarRuleBeginOffsetsIterator
>
void CalculateGrammarRules
(
  Text const &text,
  GrammarRuleSizes const &sizes,
  GrammarRules &rules,
  GrammarRuleBeginOffsetsIterator begin_offsets_begin
)
{
  rules.resize(std::accumulate(std::begin(sizes), std::end(sizes), 0));
  auto rules_iterator {std::begin(rules)};
  auto begin_offsets_iterator {begin_offsets_begin};
  auto sizes_iterator {std::begin(sizes)};
  auto sizes_end {std::end(sizes)};
  while (sizes_iterator != sizes_end)
  {
    auto rule_iterator {std::next(std::begin(text), *begin_offsets_iterator)};
    auto rule_end {std::next(std::begin(text), *begin_offsets_iterator + *sizes_iterator)};
    while (rule_iterator != rule_end)
    {
      *rules_iterator = *rule_iterator;
      ++rules_iterator;
      ++rule_iterator;
    }
    ++begin_offsets_iterator;
    ++sizes_iterator;
  }
  return;
}

template <typename GrammarRulesIterator>
void PrintGrammarTrie
(
  GrammarRulesIterator rules_begin,
  TrieNode *current_node,
  int64_t const edge_offset,
  uint64_t depth = 0
)
{
  if (current_node != nullptr)
  {
    std::cout << depth << ":";
    if (current_node->edge_begin_offset != current_node->edge_end_offset)
    {
      auto edge_iterator {std::next(rules_begin, current_node->edge_begin_offset)};
      auto edge_end {std::next(rules_begin, current_node->edge_end_offset)};
      while (edge_iterator != edge_end)
      {
        std::cout << *edge_iterator;
        edge_iterator += edge_offset;
      }
    }
    std::cout << "(" << current_node->leftmost_rank << "," << current_node->rightmost_rank << ")\n";
    auto branches_iterator {std::begin(current_node->branches)};
    auto branches_end {std::end(current_node->branches)};
    while (branches_iterator != branches_end)
    {
      PrintGrammarTrie(rules_begin, std::get<1>(*branches_iterator), edge_offset, depth + 1);
      ++branches_iterator;
    }
  }
  return;
}

template <typename GrammarRulesIterator>
void InsertGrammarRule
(
  TrieNode *current_node,
  GrammarRulesIterator rules_begin,
  GrammarRulesIterator rule_iterator,
  GrammarRulesIterator rule_end,
  int64_t const offset,
  uint64_t const lex_rank
)
{
  while (rule_iterator != rule_end)
  {
    auto character {*rule_iterator};
    auto branches_iterator {current_node->branches.find(character)};
    if (branches_iterator == std::end(current_node->branches))
    {
      current_node->branches[character] = new TrieNode
      {
        std::distance(rules_begin, rule_iterator),
        std::distance(rules_begin, rule_end),
        lex_rank
      };
      return;
    }
    else
    {
      auto child_node {std::get<1>(*branches_iterator)};
      auto edge_iterator {std::next(rules_begin, child_node->edge_begin_offset)};
      auto edge_end {std::next(rules_begin, child_node->edge_end_offset)};
      while
      (
        (rule_iterator != rule_end)
        &&
        (edge_iterator != edge_end)
        &&
        (*rule_iterator == *edge_iterator)
      )
      {
        rule_iterator += offset;
        edge_iterator += offset;
      }
      if (edge_iterator == edge_end)
      {
        if (rule_iterator != rule_end)
        {
          current_node = child_node;
        }
        else
        {
          child_node->leftmost_rank = lex_rank;
          child_node->rightmost_rank = lex_rank;
        }
      }
      else
      {
        auto edge_it_offset {std::distance(rules_begin, edge_iterator)};
        auto internal_node {new TrieNode {child_node->edge_begin_offset, edge_it_offset, 0}};
        child_node->edge_begin_offset = edge_it_offset;
        internal_node->branches[*edge_iterator] = child_node;
        if (rule_iterator != rule_end)
        {
          internal_node->branches[*rule_iterator] = new TrieNode
          {
            std::distance(rules_begin, rule_iterator),
            std::distance(rules_begin, rule_end),
            lex_rank
          };
        }
        else
        {
          internal_node->leftmost_rank = lex_rank;
          internal_node->rightmost_rank = lex_rank;
        }
        std::get<1>(*branches_iterator) = internal_node;
        return;
      }
    }
  }
  return;
}

template
<
  typename GrammarRuleSizes,
  typename GrammarRules,
  typename TrieNodePointer
>
void InsertGrammarRules
(
  GrammarRuleSizes const &grammar_rule_sizes,
  GrammarRules const &grammar_rules,
  TrieNodePointer &lex_grammar_trie_root,
  TrieNodePointer &colex_grammar_trie_root,
  uint64_t const rank_offset
)
{
  lex_grammar_trie_root = new TrieNode {};
  colex_grammar_trie_root = new TrieNode {};
  uint64_t rank {1};
  auto grammar_rule_sizes_iterator {std::begin(grammar_rule_sizes)};
  auto grammar_rules_iterator {std::begin(grammar_rules)};
  while (grammar_rules_iterator != std::end(grammar_rules))
  {
    auto grammar_rule_begin {grammar_rules_iterator};
    auto grammar_rule_end {std::next(grammar_rule_begin, *grammar_rule_sizes_iterator)};
    InsertGrammarRule
    (
      lex_grammar_trie_root,
      std::begin(grammar_rules),
      grammar_rule_begin,
      grammar_rule_end,
      1,
      rank
    );
    InsertGrammarRule
    (
      colex_grammar_trie_root,
      std::begin(grammar_rules),
      std::prev(grammar_rule_end),
      std::prev(grammar_rule_begin),
      -1,
      rank
    );
    rank += rank_offset;
    ++grammar_rule_sizes_iterator;
    grammar_rules_iterator = grammar_rule_end;
  }
  return;
}

void CalculateLexicographicTrieRankRanges (TrieNode *current_node)
{
  if (!current_node->branches.empty())
  {
    auto branches_begin {std::begin(current_node->branches)};
    auto branches_iterator {branches_begin};
    auto branches_end {std::end(current_node->branches)};
    while (branches_iterator != branches_end)
    {
      CalculateLexicographicTrieRankRanges(std::get<1>(*branches_iterator));
      ++branches_iterator;
    }
    auto first_child_node_node {std::get<1>(*branches_begin)};
    auto last_child_node_node {std::get<1>(*std::prev(branches_end))};
    if (current_node->leftmost_rank == 0)
    {
      current_node->leftmost_rank = first_child_node_node->leftmost_rank;
    }
    current_node->rightmost_rank = last_child_node_node->rightmost_rank;
  }
  return;
}

template <typename LexToColexOrderMapping>
void CalculateColexTrieRankRangesAndLexToColexOrderMapping
(
  TrieNode *current_node,
  LexToColexOrderMapping &lex_to_colex_order_mapping,
  uint64_t &colex_rank
)
{
  if (current_node->leftmost_rank != 0)
  {
    lex_to_colex_order_mapping[current_node->leftmost_rank] = colex_rank;
    current_node->leftmost_rank = current_node->rightmost_rank = colex_rank;
    ++colex_rank;
  }
  if (!current_node->branches.empty())
  {
    auto branches_begin {std::begin(current_node->branches)};
    auto branches_iterator {branches_begin};
    auto branches_end {std::end(current_node->branches)};
    while (branches_iterator != branches_end)
    {
      CalculateColexTrieRankRangesAndLexToColexOrderMapping
      (
        std::get<1>(*branches_iterator),
        lex_to_colex_order_mapping,
        colex_rank
      );
      ++branches_iterator;
    }
    auto first_child_node {std::get<1>(*branches_begin)};
    auto last_child_node {std::get<1>(*std::prev(branches_end))};
    if (current_node->leftmost_rank == 0)
    {
      current_node->leftmost_rank = first_child_node->leftmost_rank;
    }
    current_node->rightmost_rank = last_child_node->rightmost_rank;
  }
  return;
}

template
<
  typename LexToColexOrderMapping,
  typename ColexToLexOrderMapping
>
void CalculateColexToLexOrderMapping
(
  LexToColexOrderMapping const &lex_to_colex_order_mapping,
  ColexToLexOrderMapping &colex_to_lex_order_mapping
)
{
  for (uint64_t lex_rank {0}; lex_rank != std::size(lex_to_colex_order_mapping); ++lex_rank)
  {
    colex_to_lex_order_mapping[lex_to_colex_order_mapping[lex_rank]] = lex_rank;
  }
  return;
}

template
<
  typename GrammarCompressedText,
  typename LexToColexOrderMapping,
  typename ColexGrammarCompressedBwt
>
void CalculateColexGrammarCompressedBurrowWheelerTransform
(
  GrammarCompressedText const &grammar_compressed_text,
  LexToColexOrderMapping const &lex_to_colex_order_mapping,
  ColexGrammarCompressedBwt &colex_grammar_compressed_bwt
)
{
  sdsl::int_vector<> buffer;
  buffer.width(sdsl::bits::hi(std::size(grammar_compressed_text)) + 1);
  buffer.resize(std::size(grammar_compressed_text));
  sdsl::qsufsort::construct_sa(buffer, grammar_compressed_text);
  auto buffer_iterator {std::begin(buffer)};
  auto buffer_end {std::end(buffer)};
  while (buffer_iterator != buffer_end)
  {
    if (*buffer_iterator != 0)
    {
      *buffer_iterator =
      lex_to_colex_order_mapping
      [
        grammar_compressed_text[(*buffer_iterator - 1)]
      ];
    }
    ++buffer_iterator;
  }
  sdsl::construct_im(colex_grammar_compressed_bwt, buffer);
  return;
}

void CalculateTrieRankRanges
(
  TrieNode *current_node,
  uint64_t &rank
)
{
  if (current_node->leftmost_rank != 0)
  {
    current_node->leftmost_rank = current_node->rightmost_rank = ++rank;
  }
  if (!current_node->branches.empty())
  {
    auto branches_begin {std::begin(current_node->branches)};
    auto branches_iterator {branches_begin};
    auto branches_end {std::end(current_node->branches)};
    while (branches_iterator != branches_end)
    {
      CalculateTrieRankRanges(std::get<1>(*branches_iterator), rank);
      ++branches_iterator;
    }
    auto first_child_node {std::get<1>(*branches_begin)};
    auto last_child_node {std::get<1>(*std::prev(branches_end))};
    if (current_node->leftmost_rank == 0)
    {
      current_node->leftmost_rank = first_child_node->leftmost_rank;
    }
    current_node->rightmost_rank = last_child_node->rightmost_rank;
  }
  return;
}

struct Index
{
  sdsl::int_vector<> grammar_rule_sizes;
  sdsl::int_vector<> grammar_rules;
  TrieNode *lex_grammar_trie_root;
  TrieNode *colex_grammar_trie_root;
  sdsl::int_vector<> colex_to_lex_order_mapping;
  sdsl::int_vector<> lex_grammar_compressed_character_bucket_end_offsets;
  sdsl::wm_int<> colex_grammar_compressed_bwt;

  ~Index ()
  {
    delete lex_grammar_trie_root;
    delete colex_grammar_trie_root;
  }
};

void PrintIndex (Index &index)
{
  std::cout << index.grammar_rule_sizes << "\n";
  std::cout << index.grammar_rules << "\n";
  PrintGrammarTrie(std::begin(index.grammar_rules), index.lex_grammar_trie_root, 1);
  PrintGrammarTrie(std::begin(index.grammar_rules), index.colex_grammar_trie_root, -1);
  std::cout << index.colex_to_lex_order_mapping << "\n";
  std::cout << index.lex_grammar_compressed_character_bucket_end_offsets << "\n";
  std::cout << index.colex_grammar_compressed_bwt << "\n";
  return;
}

void Construct
(
  Index &index,
  std::filesystem::path const &text_path
)
{
  sdsl::int_vector<8> text;
  {
    std::ifstream text_file {text_path};
    sdsl::load_vector_from_file(text, text_path);
    sdsl::append_zero_symbol(text);
  }
  sdsl::bit_vector sl_types;
  {
    CalculateSlTypes(text, sl_types);
  }
  auto invalid_text_offset {std::size(text)};
  auto text_size_width {sdsl::bits::hi(std::size(text)) + 1};
  sdsl::int_vector<> text_offsets;
  {
    text_offsets.width(text_size_width);
    text_offsets.resize(std::size(text));
    sdsl::util::set_to_value(text_offsets, invalid_text_offset);
  }
  sdsl::int_vector<> character_bucket_offsets;
  {
    character_bucket_offsets.width(text_size_width);
    character_bucket_offsets.resize(256);
  }
  {
    CalculateCharacterBucketBeginOffsets(text, character_bucket_offsets);
    BucketSortRightmostLTypeCharacters(text, sl_types, character_bucket_offsets, text_offsets);
    InduceSortLTypeCharacters(text, sl_types, character_bucket_offsets, text_offsets);
    CalculateCharacterBucketEndOffsets(text, character_bucket_offsets);
    InduceSortSTypeCharacters(text, sl_types, character_bucket_offsets, text_offsets);
  }
  auto text_offsets_boundary {std::begin(text_offsets)};
  text_offsets_boundary = MoveVaildEntriesToFront
  (
    std::begin(text_offsets),
    std::end(text_offsets),
    invalid_text_offset
  );
  auto grammar_rule_begin_offsets_begin {std::begin(text_offsets)};
  auto grammar_rule_begin_offsets_end {text_offsets_boundary};
  auto temporary_grammar_compressed_text_begin {text_offsets_boundary};
  auto temporary_grammar_compressed_text_end {std::end(text_offsets)};
  CalculateGrammarRuleBeginOffsetsAndTemporaryGrammarCompressedText
  (
    text,
    sl_types,
    grammar_rule_begin_offsets_begin,
    grammar_rule_begin_offsets_end,
    temporary_grammar_compressed_text_begin,
    temporary_grammar_compressed_text_end
  );
  CalculateGrammarRuleSizes
  (
    index.grammar_rule_sizes,
    sl_types,
    grammar_rule_begin_offsets_begin,
    grammar_rule_begin_offsets_end
  );
  CalculateGrammarRules
  (
    text,
    index.grammar_rule_sizes,
    index.grammar_rules,
    grammar_rule_begin_offsets_begin
  );
  sdsl::util::clear(sl_types);
  sdsl::util::clear(text);
  InsertGrammarRules
  (
    index.grammar_rule_sizes,
    index.grammar_rules,
    index.lex_grammar_trie_root,
    index.colex_grammar_trie_root,
    1
  );
  CalculateLexicographicTrieRankRanges(index.lex_grammar_trie_root);
  auto grammar_compressed_alphabet_size {std::size(index.grammar_rule_sizes) + 1};
  auto grammar_compressed_character_width {sdsl::bits::hi(grammar_compressed_alphabet_size) + 1};
  sdsl::int_vector<> lex_to_colex_order_mapping;
  lex_to_colex_order_mapping.width(grammar_compressed_character_width);
  lex_to_colex_order_mapping.resize(grammar_compressed_alphabet_size);
  lex_to_colex_order_mapping[0] = 0;
  {
    uint64_t colex_rank {1};
    CalculateColexTrieRankRangesAndLexToColexOrderMapping
    (
      index.colex_grammar_trie_root,
      lex_to_colex_order_mapping,
      colex_rank
    );
  }
  index.colex_to_lex_order_mapping.width(grammar_compressed_character_width);
  index.colex_to_lex_order_mapping.resize(grammar_compressed_alphabet_size);
  CalculateColexToLexOrderMapping
  (
    lex_to_colex_order_mapping,
    index.colex_to_lex_order_mapping
  );
  auto grammar_compressed_text_size
  {
    1 + std::distance
    (
      temporary_grammar_compressed_text_begin,
      temporary_grammar_compressed_text_end
    )
  };
  sdsl::int_vector<> grammar_compressed_text;
  grammar_compressed_text.width(grammar_compressed_character_width);
  grammar_compressed_text.resize(grammar_compressed_text_size);
  std::copy
  (
    temporary_grammar_compressed_text_begin,
    temporary_grammar_compressed_text_end,
    std::begin(grammar_compressed_text)
  );
  sdsl::util::clear(text_offsets);
  grammar_compressed_text[std::size(grammar_compressed_text) - 1] = 0;
  index.lex_grammar_compressed_character_bucket_end_offsets.width(sdsl::bits::hi(grammar_compressed_text_size) + 1);
  index.lex_grammar_compressed_character_bucket_end_offsets.resize(grammar_compressed_alphabet_size);
  CalculateCharacterBucketEndOffsets
  (
    grammar_compressed_text,
    index.lex_grammar_compressed_character_bucket_end_offsets
  );
  CalculateColexGrammarCompressedBurrowWheelerTransform
  (
    grammar_compressed_text,
    lex_to_colex_order_mapping,
    index.colex_grammar_compressed_bwt
  );
  return;
}

void Serialize
(
  Index &index,
  std::filesystem::path index_path
)
{
  std::ofstream index_file {index_path};
  sdsl::serialize(index.grammar_rule_sizes, index_file);
  sdsl::serialize(index.grammar_rules, index_file);
  sdsl::serialize(index.colex_to_lex_order_mapping, index_file);
  sdsl::serialize(index.lex_grammar_compressed_character_bucket_end_offsets, index_file);
  sdsl::serialize(index.colex_grammar_compressed_bwt, index_file);
  return;
}

void Load
(
  Index &index,
  std::filesystem::path index_path
)
{
  sdsl::util::clear(index.grammar_rule_sizes);
  sdsl::util::clear(index.grammar_rules);
  sdsl::util::clear(index.colex_to_lex_order_mapping);
  sdsl::util::clear(index.lex_grammar_compressed_character_bucket_end_offsets);
  sdsl::util::clear(index.colex_grammar_compressed_bwt);
  load_from_file(index.grammar_rule_sizes, index_path);
  load_from_file(index.grammar_rules, index_path);
  load_from_file(index.colex_to_lex_order_mapping, index_path);
  load_from_file(index.lex_grammar_compressed_character_bucket_end_offsets, index_path);
  load_from_file(index.colex_grammar_compressed_bwt, index_path);
  InsertGrammarRules
  (
    index.grammar_rule_sizes,
    index.grammar_rules,
    index.lex_grammar_trie_root,
    index.colex_grammar_trie_root,
    0
  );
  {
    uint64_t lex_rank {0};
    CalculateTrieRankRanges(index.lex_grammar_trie_root, lex_rank);
  }
  {
    uint64_t colex_rank {0};
    CalculateTrieRankRanges(index.colex_grammar_trie_root, colex_rank);
  }
  return;
}

template <typename TextIterator>
void CalculateSlFactor
(
  TextIterator &reverse_first_iterator,
  TextIterator &reverse_last_iterator,
  TextIterator reverse_end
)
{
  uint64_t previous_sl_type {L};
  reverse_first_iterator = reverse_last_iterator--;
  while
  (
    (reverse_last_iterator != reverse_end)
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

template
<
  typename GrammarRulesIterator,
  typename SlFactorIterator
>
void LookupGrammarRule
(
  GrammarRulesIterator grammar_rules_begin,
  TrieNode *current_node,
  SlFactorIterator sl_factor_iterator,
  SlFactorIterator sl_factor_end,
  int64_t const offset,
  uint64_t &leftmost_rank,
  uint64_t &rightmost_rank
)
{
  leftmost_rank = rightmost_rank = 0;
  auto branches_iterator {std::end(current_node->branches)};
  while
  (
    (branches_iterator = current_node->branches.find(*sl_factor_iterator))
    !=
    std::end(current_node->branches)
  )
  {
    current_node = std::get<1>(*branches_iterator);
    auto edge_iterator {std::next(grammar_rules_begin, current_node->edge_begin_offset)};
    auto edge_end {std::next(grammar_rules_begin, current_node->edge_end_offset)};
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
        return;
      }
    }
    else
    {
      rightmost_rank = current_node->rightmost_rank;
      leftmost_rank = current_node->leftmost_rank;
      return;
    }
  }
  return;
}

template <typename PatternIterator>
void BackwardSearchPatternPrefix
(
  Index const &index,
  uint64_t &pattern_range_begin_offset_L,
  uint64_t &pattern_range_end_offset_L,
  uint64_t &pattern_range_begin_offset_S,
  uint64_t &pattern_range_end_offset_S,
  PatternIterator reverse_pattern_first_iterator,
  PatternIterator reverse_pattern_last_iterator
)
{
  uint64_t leftmost_colex_rank {0};
  uint64_t rightmost_colex_rank {0};
  LookupGrammarRule
  (
    std::begin(index.grammar_rules),
    index.colex_grammar_trie_root,
    reverse_pattern_first_iterator,
    reverse_pattern_last_iterator,
    -1,
    leftmost_colex_rank,
    rightmost_colex_rank
  );
  if (leftmost_colex_rank != 0)
  {
    if (pattern_range_begin_offset_L != pattern_range_end_offset_L)
    {
      pattern_range_end_offset_L = std::get<0>
      (
        index.colex_grammar_compressed_bwt.range_search_2d
        (
          pattern_range_begin_offset_L,
          (pattern_range_end_offset_L - 1),
          leftmost_colex_rank,
          rightmost_colex_rank,
          false
        )
      );
      pattern_range_begin_offset_L = 0;
    }
    if (pattern_range_begin_offset_S != pattern_range_end_offset_S)
    {
      pattern_range_end_offset_S = std::get<0>
      (
        index.colex_grammar_compressed_bwt.range_search_2d
        (
          pattern_range_begin_offset_S,
          (pattern_range_end_offset_S - 1),
          leftmost_colex_rank,
          rightmost_colex_rank,
          false
        )
      );
      pattern_range_begin_offset_S = 0;
    }
  }
  // {
  //   auto iterator {std::next(reverse_pattern_last_iterator)};
  //   auto end {std::next(reverse_pattern_first_iterator)};
  //   while ( iterator != end)
  //   {
  //     std::cout << *iterator;
  //     ++iterator;
  //   }
  //   std::cout
  //   << "->(" << leftmost_colex_rank << "," << rightmost_colex_rank << ")"
  //   << "->L:(" << pattern_range_begin_offset_L << "," << pattern_range_end_offset_L << ")"
  //   << "->S:(" << pattern_range_begin_offset_S << "," << pattern_range_end_offset_S << ")\n";
  // }
  return;
}

template <typename PatternIterator>
void BackwardSearchExactSlFactor
(
  Index const &index,
  uint64_t &pattern_range_begin_offset_L,
  uint64_t &pattern_range_end_offset_L,
  uint64_t &pattern_range_begin_offset_S,
  uint64_t &pattern_range_end_offset_S,
  PatternIterator reverse_pattern_first_iterator,
  PatternIterator reverse_pattern_last_iterator
)
{
  uint64_t colex_rank {0};
  LookupGrammarRule
  (
    std::begin(index.grammar_rules),
    index.colex_grammar_trie_root,
    reverse_pattern_first_iterator,
    reverse_pattern_last_iterator,
    -1,
    colex_rank,
    colex_rank
  );
  if (colex_rank != 0)
  {
    auto character_begin_offset
    {
      index.lex_grammar_compressed_character_bucket_end_offsets
      [
        index.colex_to_lex_order_mapping[colex_rank] - 1
      ]
    };
    if (pattern_range_begin_offset_L != pattern_range_end_offset_L)
    {
      pattern_range_begin_offset_L =
      (
        character_begin_offset
        + index.colex_grammar_compressed_bwt.rank
        (
          pattern_range_begin_offset_L,
          colex_rank
        )
      );
      pattern_range_end_offset_L =
      (
        character_begin_offset
        + index.colex_grammar_compressed_bwt.rank
        (
          pattern_range_end_offset_L,
          colex_rank
        )
      );
    }
    if (pattern_range_begin_offset_S != pattern_range_end_offset_S)
    {
      pattern_range_begin_offset_S =
      (
        character_begin_offset
        + index.colex_grammar_compressed_bwt.rank
        (
          pattern_range_begin_offset_S,
          colex_rank
        )
      );
      pattern_range_end_offset_S =
      (
        character_begin_offset
        + index.colex_grammar_compressed_bwt.rank
        (
          pattern_range_end_offset_S,
          colex_rank
        )
      );
    }
  }
  else
  {
    pattern_range_begin_offset_L = pattern_range_end_offset_L;
    pattern_range_begin_offset_S = pattern_range_end_offset_S;
  }
  // {
  //   auto iterator {std::next(reverse_pattern_last_iterator)};
  //   auto end {std::next(reverse_pattern_first_iterator)};
  //   while (iterator != end)
  //   {
  //     std::cout << *iterator;
  //     ++iterator;
  //   }
  //   std::cout
  //   << "->(" << colex_rank << ":" << index.colex_to_lex_order_mapping[colex_rank] << ")"
  //   << "->L:(" << pattern_range_begin_offset_L << "," << pattern_range_end_offset_L << ")"
  //   << "->S:(" << pattern_range_begin_offset_S << "," << pattern_range_end_offset_S << ")\n";
  // }
  return;
}

template <typename PatternIterator>
auto CalculateReversePatternLastIteratorOfPatternSuffixS
(
  PatternIterator reverse_begin,
  PatternIterator reverse_end
)
{
  auto iterator {reverse_begin};
  while
  (
    (std::prev(iterator) != reverse_end)
    &&
    (*std::prev(iterator) == *iterator)
  )
  {
    --iterator;
  }
  if
  (
    (std::prev(iterator) != reverse_end)
    &&
    (*std::prev(iterator) < *iterator)
  )
  {
    return reverse_begin;
  }
  return std::prev(iterator);
}

template <typename PatternIterator>
auto BackwardSearchPatternSuffix
(
  Index const &index,
  uint64_t &pattern_range_begin_offset_L,
  uint64_t &pattern_range_end_offset_L,
  uint64_t &pattern_range_begin_offset_S,
  uint64_t &pattern_range_end_offset_S,
  PatternIterator reverse_pattern_begin,
  PatternIterator reverse_pattern_end
)
{
  auto reverse_pattern_first_iterator {reverse_pattern_begin};
  auto reverse_pattern_last_iterator
  {
    CalculateReversePatternLastIteratorOfPatternSuffixS
    (
      reverse_pattern_begin,
      reverse_pattern_end
    )
  };
  uint64_t leftmost_lex_rank {0};
  uint64_t rightmost_lex_rank {0};
  if
  (
    (reverse_pattern_last_iterator != reverse_pattern_begin)
    &&
    (reverse_pattern_last_iterator != reverse_pattern_end)
  )
  {
    LookupGrammarRule
    (
      std::begin(index.grammar_rules),
      index.lex_grammar_trie_root,
      std::next(reverse_pattern_last_iterator),
      std::next(reverse_pattern_begin),
      1,
      leftmost_lex_rank,
      rightmost_lex_rank
    );
    if (leftmost_lex_rank != 0)
    {
      pattern_range_begin_offset_S = index.lex_grammar_compressed_character_bucket_end_offsets[leftmost_lex_rank - 1];
      pattern_range_end_offset_S = index.lex_grammar_compressed_character_bucket_end_offsets[rightmost_lex_rank];
    }
    // {
    //   auto iterator {std::next(reverse_pattern_last_iterator)};
    //   auto end {std::next(reverse_pattern_begin)};
    //   while (iterator != end)
    //   {
    //     std::cout << *iterator;
    //     ++iterator;
    //   }
    //   std::cout << "->(" << leftmost_lex_rank << "," << rightmost_lex_rank << ")";
    //   std::cout << "->S:(" << pattern_range_begin_offset_S << "," << pattern_range_end_offset_S << ")\n";
    // }
    CalculateSlFactor(reverse_pattern_first_iterator, reverse_pattern_last_iterator, reverse_pattern_end);
    if (pattern_range_begin_offset_S != pattern_range_end_offset_S)
    {
      uint64_t leftmost_colex_rank {0};
      uint64_t rightmost_colex_rank {0};
      LookupGrammarRule
      (
        std::begin(index.grammar_rules),
        index.colex_grammar_trie_root,
        reverse_pattern_first_iterator,
        reverse_pattern_last_iterator,
        -1,
        leftmost_colex_rank,
        rightmost_colex_rank
      );
      if (leftmost_colex_rank != 0)
      {
        if (reverse_pattern_last_iterator != reverse_pattern_end)
        {
          auto character_begin_offset
          {
            index.lex_grammar_compressed_character_bucket_end_offsets
            [
              index.colex_to_lex_order_mapping[leftmost_colex_rank] - 1
            ]
          };
          pattern_range_begin_offset_S =
          (
            character_begin_offset
            + index.colex_grammar_compressed_bwt.rank
            (
              pattern_range_begin_offset_S,
              leftmost_colex_rank
            )
          );
          pattern_range_end_offset_S =
          (
            character_begin_offset
            + index.colex_grammar_compressed_bwt.rank
            (
              pattern_range_end_offset_S,
              leftmost_colex_rank
            )
          );
        }
        else
        {
          pattern_range_end_offset_S = std::get<0>
          (
            index.colex_grammar_compressed_bwt.range_search_2d
            (
              pattern_range_begin_offset_S,
              (pattern_range_end_offset_S - 1),
              leftmost_colex_rank,
              rightmost_colex_rank,
              false
            )
          );
          pattern_range_begin_offset_S = 0;
        }
      }
      else
      {
        pattern_range_begin_offset_S = pattern_range_end_offset_S;
      }
      // {
      //   auto iterator {std::next(reverse_pattern_last_iterator)};
      //   auto end {std::next(reverse_pattern_first_iterator)};
      //   while (iterator != end)
      //   {
      //     std::cout << *iterator;
      //     ++iterator;
      //   }
      //   std::cout
      //   << "->(" << leftmost_colex_rank << ":" << index.colex_to_lex_order_mapping[leftmost_colex_rank]
      //   << "," << rightmost_colex_rank << ")"
      //   << "->S:(" << pattern_range_begin_offset_S << "," << pattern_range_end_offset_S << ")\n";
      // }
    }
  }
  else if (reverse_pattern_last_iterator == reverse_pattern_begin)
  {
    CalculateSlFactor
    (
      reverse_pattern_first_iterator,
      reverse_pattern_last_iterator,
      reverse_pattern_end
    );
  }
  LookupGrammarRule
  (
    std::begin(index.grammar_rules),
    index.lex_grammar_trie_root,
    std::next(reverse_pattern_last_iterator),
    std::next(reverse_pattern_begin),
    1,
    leftmost_lex_rank,
    rightmost_lex_rank
  );
  if (leftmost_lex_rank != 0)
  {
    pattern_range_begin_offset_L = index.lex_grammar_compressed_character_bucket_end_offsets[leftmost_lex_rank - 1];
    pattern_range_end_offset_L = index.lex_grammar_compressed_character_bucket_end_offsets[rightmost_lex_rank];
  }
  // {
  //   auto iterator {std::next(reverse_pattern_last_iterator)};
  //   auto end {std::next(reverse_pattern_begin)};
  //   while (iterator != end)
  //   {
  //     std::cout << *iterator;
  //     ++iterator;
  //   }
  //   std::cout << "->(" << leftmost_lex_rank << "," << rightmost_lex_rank << ")";
  //   std::cout << "->L:(" << pattern_range_begin_offset_L << "," << pattern_range_end_offset_L << ")\n";
  // }
  return reverse_pattern_last_iterator;
}

template <typename PatternIterator>
uint64_t Count
(
  Index const &index,
  PatternIterator pattern_begin,
  PatternIterator pattern_end
)
{
  uint64_t pattern_range_begin_offset_L {0};
  uint64_t pattern_range_end_offset_L {0};
  uint64_t pattern_range_begin_offset_S {0};
  uint64_t pattern_range_end_offset_S {0};
  auto reverse_pattern_begin {std::prev(pattern_end)};
  auto reverse_pattern_end {std::prev(pattern_begin)};
  auto reverse_pattern_first_iterator {reverse_pattern_begin};
  auto reverse_pattern_last_iterator {reverse_pattern_begin};
  // {
  //   for (auto iterator {pattern_begin}; iterator != pattern_end; ++iterator)
  //   {
  //     std::cout << *iterator;
  //   }
  //   std::cout << "\n";
  // }
  reverse_pattern_last_iterator = BackwardSearchPatternSuffix
  (
    index,
    pattern_range_begin_offset_L,
    pattern_range_end_offset_L,
    pattern_range_begin_offset_S,
    pattern_range_end_offset_S,
    reverse_pattern_begin,
    reverse_pattern_end
  );
  if (reverse_pattern_last_iterator != reverse_pattern_end)
  {
    while
    (
      (pattern_range_begin_offset_L != pattern_range_end_offset_L)
      ||
      (pattern_range_begin_offset_S != pattern_range_end_offset_S)
    )
    {
      CalculateSlFactor(reverse_pattern_first_iterator, reverse_pattern_last_iterator, reverse_pattern_end);
      if (reverse_pattern_last_iterator!= reverse_pattern_end)
      {
        BackwardSearchExactSlFactor
        (
          index,
          pattern_range_begin_offset_L,
          pattern_range_end_offset_L,
          pattern_range_begin_offset_S,
          pattern_range_end_offset_S,
          reverse_pattern_first_iterator,
          reverse_pattern_last_iterator
        );
      }
      else
      {
        BackwardSearchPatternPrefix
        (
          index,
          pattern_range_begin_offset_L,
          pattern_range_end_offset_L,
          pattern_range_begin_offset_S,
          pattern_range_end_offset_S,
          reverse_pattern_first_iterator,
          reverse_pattern_last_iterator
        );
        break;
      }
    }
  }
  return
  (
    (pattern_range_end_offset_L - pattern_range_begin_offset_L)
    +
    (pattern_range_end_offset_S - pattern_range_begin_offset_S)
  );
}

template <typename Text>
uint64_t CalculateMaxSlFactorSize (Text const &text)
{
  uint64_t max_size {0};
  auto reverse_first_iterator {std::prev(std::end(text))};
  auto reverse_last_iterator {std::prev(std::end(text))};
  auto reverse_end {std::prev(std::begin(text))};
  while (reverse_last_iterator != reverse_end)
  {
    CalculateSlFactor(reverse_first_iterator, reverse_last_iterator, reverse_end);
    auto size {static_cast<uint64_t>(std::distance(reverse_last_iterator, reverse_first_iterator))};
    if (max_size < size)
    {
      max_size = size;
    }
  }
  return max_size;
}
}
