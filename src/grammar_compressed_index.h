#pragma once

#include <deque>
#include <map>
#include <memory>

#include "dynamic_grammar_trie.h"
#include "static_grammar_trie.h"
#include "utility.h"

namespace project
{
constexpr uint64_t S {1};
constexpr uint64_t L {0};

template <typename SlTypesIterator>
constexpr bool IsRightmostLType (SlTypesIterator it)
{
  return ((*it == L) && (*std::next(it) == S));
}

template <typename SlTypesIterator>
constexpr bool IsLeftmostSType (SlTypesIterator it)
{
  return ((*it == S) && (*std::prev(it) == L));
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
  auto text_it {std::prev(std::end(text))};
  auto text_last_it {std::begin(text)};
  auto sl_types_it {std::prev(std::end(sl_types))};
  while (text_it != text_last_it)
  {
    if ((*std::prev(text_it) > *text_it) || ((*std::prev(text_it) == *text_it) && (*sl_types_it == L)))
    {
      *std::prev(sl_types_it) = L;
    }
    --sl_types_it;
    --text_it;
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
  sdsl::util::set_to_value(character_bucket_offsets, 0);
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
  auto it {std::prev(std::end(character_bucket_offsets))};
  auto last_it {std::begin(character_bucket_offsets)};
  while (it != last_it)
  {
    *it-- = *std::prev(it);
  }
  *last_it = 0;
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
  auto text_it {std::begin(text)};
  auto text_end {std::end(text)};
  while (text_it != text_end)
  {
    auto text_offset {std::distance(std::begin(text), text_it)};
    if (IsRightmostLType(std::next(std::begin(sl_types), text_offset)))
    {
      text_offsets[character_bucket_offsets[*text_it]++] = text_offset;
    }
    ++text_it;
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
  auto text_offsets_it {std::next(std::begin(text_offsets))};
  auto text_offsets_last_it {std::end(text_offsets)};
  auto invalid_text_offset {std::size(text)};
  while (text_offsets_it != text_offsets_last_it)
  {
    auto text_offset {*text_offsets_it};
    if((text_offset != invalid_text_offset) && (text_offset != 0) && (sl_types[text_offset - 1] == L))
    {
      text_offsets[character_bucket_offsets[text[text_offset - 1]]++] = (text_offset - 1);
      *text_offsets_it = invalid_text_offset;
    }
    ++text_offsets_it;
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
  auto text_offsets_it {std::prev(std::end(text_offsets))};
  auto text_offsets_last_it {std::begin(text_offsets)};
  auto invalid_text_offset {std::size(text)};
  while (text_offsets_it != text_offsets_last_it)
  {
    auto text_offset {*text_offsets_it};
    if ((text_offset != invalid_text_offset) && (text_offset != 0) && (sl_types[text_offset - 1] == S))
    {
      text_offsets[--character_bucket_offsets[text[text_offset - 1]]] = (text_offset - 1);
      *text_offsets_it = invalid_text_offset;
    }
    --text_offsets_it;
  }
  return;
}

template <typename Iterator>
auto MoveVaildEntriesToFront
(
  Iterator begin,
  Iterator end,
  uint64_t const invalid_value
)
{
  auto it {begin};
  auto last {begin};
  while (it != end)
  {
    if (*it != invalid_value)
    {
      uint64_t value {*it}; *it = *last; *last = value;
      ++last;
    }
    ++it;
  }
  return last;
}

template
<
  typename Text,
  typename SlTypes,
  typename GrammarRuleBeginOffsetsIterator,
  typename TemporaryLexTextIterator,
  typename GrammarRuleCounts
>
void CalculateGrammarRuleCountsBeginOffsetsAndTemporaryLexText
(
  Text const &text,
  SlTypes const &sl_types,
  GrammarRuleCounts &rule_counts,
  GrammarRuleBeginOffsetsIterator begin_offsets_begin,
  GrammarRuleBeginOffsetsIterator begin_offsets_end,
  TemporaryLexTextIterator temporary_compressed_text_begin
)
{
  std::deque<uint64_t> counts;
  auto invalid_text_offset {std::size(text)};
  uint64_t lex_rank {};
  auto prev_rule_it {std::prev(std::end(text))};
  auto prev_sl_types_it {std::prev(std::end(sl_types))};
  auto begin_offsets_it {begin_offsets_begin};
  while (begin_offsets_it != begin_offsets_end)
  {
    uint64_t begin_offset {*begin_offsets_it};
    auto rule_it {std::next(std::begin(text), begin_offset)};
    auto sl_types_it {std::next(std::begin(sl_types), begin_offset)};
    if ((*prev_rule_it == *rule_it) && (*prev_sl_types_it == *sl_types_it))
    {
      do
      {
        ++prev_rule_it;
        ++rule_it;
        ++prev_sl_types_it;
        ++sl_types_it;
      }
      while
      (
        !IsLeftmostSType(prev_sl_types_it)
        &&
        !IsLeftmostSType(sl_types_it)
        &&
        (*prev_rule_it == *rule_it)
        &&
        (*prev_sl_types_it == *sl_types_it)
      );
      if (IsLeftmostSType(prev_sl_types_it) && IsLeftmostSType(sl_types_it))
      {
        --lex_rank;
        *begin_offsets_it = invalid_text_offset;
        ++counts.back();
      }
      else
      {
        counts.emplace_back(1);
      }
    }
    else
    {
      counts.emplace_back(1);
    }
    *std::next(temporary_compressed_text_begin, (begin_offset + 1) / 2) = ++lex_rank;
    prev_rule_it = std::next(std::begin(text), begin_offset);
    prev_sl_types_it = std::next(std::begin(sl_types), begin_offset);
    ++begin_offsets_it;
  }
  rule_counts.resize(std::size(counts));
  auto rule_counts_it {std::begin(rule_counts)};
  for (auto const &count : counts)
  {
    *rule_counts_it++ = count;
  }
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
  GrammarRuleBeginOffsetsIterator begin_offsets_end,
  uint64_t const invalid_value
)
{
  begin_offsets_end = MoveVaildEntriesToFront(begin_offsets_begin, begin_offsets_end, invalid_value);
  sizes.resize(std::distance(begin_offsets_begin, begin_offsets_end));
  auto sizes_it {std::begin(sizes)};
  auto begin_offsets_it {begin_offsets_begin};
  while (begin_offsets_it != begin_offsets_end)
  {
    auto sl_types_first_it {std::next(std::begin(sl_types), *begin_offsets_it)};
    auto sl_types_last_it {std::next(sl_types_first_it)};
    while (!IsLeftmostSType(sl_types_last_it))
    {
      ++sl_types_last_it;
    }
    *sizes_it = std::distance(sl_types_first_it, sl_types_last_it);
    ++sizes_it;
    ++begin_offsets_it;
  }
  return;
}

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
  auto rules_it {std::begin(rules)};
  auto begin_offsets_it {begin_offsets_begin};
  auto sizes_it {std::begin(sizes)};
  auto sizes_end {std::end(sizes)};
  while (sizes_it != sizes_end)
  {
    auto rule_it {std::next(std::begin(text), *begin_offsets_it)};
    auto rule_end {std::next(std::begin(text), *begin_offsets_it + *sizes_it)};
    while (rule_it != rule_end)
    {
      *rules_it = *rule_it;
      ++rules_it;
      ++rule_it;
    }
    ++begin_offsets_it;
    ++sizes_it;
  }
  return;
}

template
<
  typename LexText,
  typename TemporaryLexTextIterator
>
void CalculateLexText
(
  LexText &text,
  uint64_t const width,
  TemporaryLexTextIterator begin,
  TemporaryLexTextIterator end,
  uint64_t const invalid_value
)
{
  end = MoveVaildEntriesToFront(begin, end, invalid_value);
  auto size {std::distance(begin, end) + 1};
  {
    text.width(width);
    text.resize(size);
    std::copy(begin, end, std::begin(text));
    *std::prev(std::end(text)) = 0;
  }
  return;
}

template
<
  typename LexText,
  typename LexToColex,
  typename ColexBwt
>
void CalculateColexBwt
(
  LexText const &lex_text,
  LexToColex const &lex_to_colex,
  ColexBwt &colex_bwt
)
{
  sdsl::qsufsort::construct_sa(colex_bwt, lex_text);
  auto it {std::begin(colex_bwt)};
  while (it != std::end(colex_bwt))
  {
    if (*it != 0)
    {
      *it = lex_to_colex[lex_text[(*it - 1)]];
    }
    ++it;
  }
  sdsl::util::bit_compress(colex_bwt);
  return;
}

template <typename RunBitVector = sdsl::bit_vector, typename LengthBitVector = sdsl::sd_vector<>>
struct RunLengthWaveletTree
{
  uint8_t level_size;
  uint64_t runs_size;
  uint64_t alphabet_size;
  RunBitVector all_run_bits;
  typename RunBitVector::rank_1_type all_run_bits_rank;
  LengthBitVector first_length_bits;
  typename LengthBitVector::rank_1_type first_length_bits_rank;
  typename LengthBitVector::select_1_type first_length_bits_select;
  LengthBitVector middle_length_bits;
  typename LengthBitVector::select_1_type middle_length_bits_select;
  LengthBitVector last_length_bits;
  typename LengthBitVector::select_1_type last_length_bits_select;
  sdsl::int_vector<> run_bucket_begin_offsets;
};

uint64_t Access (RunLengthWaveletTree<> const &rlwt, uint64_t offset)
{
  uint64_t run {};
  if (offset < (std::size(rlwt.first_length_bits) - 1))
  {
    offset = rlwt.first_length_bits_rank(offset + 1) - 1;
    uint64_t prefix {};
    uint64_t mask {1ULL << rlwt.level_size};
    for (uint8_t level {}; level != rlwt.level_size; ++level)
    {
      auto begin_offset {(level * rlwt.runs_size) + rlwt.run_bucket_begin_offsets[prefix]};
      auto ones_offset {rlwt.all_run_bits_rank(begin_offset + offset) - rlwt.all_run_bits_rank(begin_offset)};
      run <<= 1;
      mask >>= 1;
      if (rlwt.all_run_bits[begin_offset + offset])
      {
        run |= 1ULL;
        prefix |= mask;
        offset = ones_offset;
      }
      else
      {
        offset -= ones_offset;
      }
    }
  }
  return run;
}

uint64_t Rank
(
  RunLengthWaveletTree<> const &rlwt,
  uint64_t offset,
  uint64_t character
)
{
  if ((offset == 0) || (offset >= std::size(rlwt.first_length_bits)) || character >= rlwt.alphabet_size)
  {
    return 0;
  }
  auto runs_offset {rlwt.first_length_bits_rank(offset)};
  auto run_rank {runs_offset};
  uint64_t prefix {};
  uint64_t mask {1ULL << rlwt.level_size};
  for (uint8_t level {}; level != rlwt.level_size; ++level)
  {
    auto begin_offset {(level * rlwt.runs_size) + rlwt.run_bucket_begin_offsets[prefix]};
    auto ones_size {rlwt.all_run_bits_rank(begin_offset + run_rank) - rlwt.all_run_bits_rank(begin_offset)};
    mask >>= 1;
    if (character & mask)
    {
      prefix |= mask;
      run_rank = ones_size;
    }
    else
    {
      run_rank -= ones_size;
    }
  }
  auto run_bucket_offset {rlwt.run_bucket_begin_offsets[character]};
  auto character_bucket_begin_offset {rlwt.last_length_bits_select(run_bucket_offset + 1)};
  if (character != Access(rlwt, offset - 1))
  {
    return
    (
      rlwt.last_length_bits_select(run_bucket_offset + run_rank + 1)
      - character_bucket_begin_offset
    );
  }
  return
  (
    rlwt.last_length_bits_select(run_bucket_offset + run_rank)
    - character_bucket_begin_offset
    + (offset - rlwt.first_length_bits_select(runs_offset))
  );
}

uint64_t RangeCount
(
  RunLengthWaveletTree<> const &rlwt,
  uint64_t lower_offset,
  uint64_t upper_offset,
  uint64_t lower_value,
  uint64_t upper_value
)
{
  uint64_t count {};
  if
  (
    (lower_offset < std::size(rlwt.first_length_bits) - 1)
    && (upper_offset < std::size(rlwt.first_length_bits) - 1)
    && (lower_value < rlwt.alphabet_size)
    && (upper_value < rlwt.alphabet_size)
    && (lower_offset <= upper_offset)
    && (lower_value <= upper_value)
  )
  {
    auto lower_runs_offset {rlwt.first_length_bits_rank(++lower_offset)};
    auto upper_runs_offset {rlwt.first_length_bits_rank(++upper_offset)};
    auto run {Access(rlwt, lower_offset - 1)};
    if (lower_runs_offset == upper_runs_offset)
    {
      if ((lower_value <= run) && (run <= upper_value))
      {
        count = (upper_offset - lower_offset + 1);
      }
    }
    else
    {
      if ((lower_value <= run) && (run <= upper_value))
      {
        count += (rlwt.first_length_bits_select(lower_runs_offset + 1) - lower_offset + 1);
      }
      if (lower_runs_offset != upper_runs_offset)
      {
        run = Access(rlwt, upper_offset - 1);
        if ((lower_value <= run) && (run <= upper_value))
        {
          count += (upper_offset - rlwt.first_length_bits_select(upper_runs_offset));
        }
      }
      std::deque<std::tuple<uint8_t, uint64_t, uint64_t, uint64_t, uint64_t>> lists;
      lists.emplace_back(0, lower_runs_offset + 1, upper_runs_offset - 1, 0, rlwt.alphabet_size - 1);
      while (!lists.empty())
      {
        auto level {std::get<0>(lists.front())};
        lower_runs_offset = std::get<1>(lists.front());
        upper_runs_offset = std::get<2>(lists.front());
        auto lower_prefix {std::get<3>(lists.front())};
        auto upper_prefix {std::get<4>(lists.front())};
        lists.pop_front();
        if ((lower_runs_offset <= upper_runs_offset) && !((upper_prefix < lower_value) || (upper_value < lower_prefix)))
        {
          auto begin_offset {(level * rlwt.runs_size) + rlwt.run_bucket_begin_offsets[lower_prefix]};
          if ((lower_value <= lower_prefix) && (upper_prefix <= upper_value))
          {
            auto *select {&rlwt.first_length_bits_select};
            if (level == rlwt.level_size)
            {
              begin_offset = rlwt.run_bucket_begin_offsets[lower_prefix];
              select = &rlwt.last_length_bits_select;
            }
            else if (level != 0)
            {
              begin_offset = begin_offset - rlwt.runs_size + (level - 1);
              select = &rlwt.middle_length_bits_select;
            }
            count += (*select)(begin_offset + upper_runs_offset + 1) - (*select)(begin_offset + lower_runs_offset);
          }
          else
          {
            auto ones_begin_offset {rlwt.all_run_bits_rank(begin_offset)};
            auto lower_ones_offset {rlwt.all_run_bits_rank(begin_offset + lower_runs_offset - 1) + 1 - ones_begin_offset};
            auto upper_ones_offset {rlwt.all_run_bits_rank(begin_offset + upper_runs_offset) - ones_begin_offset};
            uint64_t middle_prefix {lower_prefix | (1ULL << (rlwt.level_size - 1 - level))};
            lists.emplace_back
            (
              level + 1,
              lower_runs_offset - lower_ones_offset + 1,
              upper_runs_offset - upper_ones_offset,
              lower_prefix,
              middle_prefix - 1
            );
            lists.emplace_back
            (
              level + 1,
              lower_ones_offset,
              upper_ones_offset,
              middle_prefix,
              upper_prefix
            );
          }
        }
      }
    }
  }
  return count;
}

template <typename File>
void PrintRunLengthWaveletTree (File &file, RunLengthWaveletTree<> const &rlwt)
{
  auto length_width {sdsl::bits::hi(std::size(rlwt.first_length_bits) - 1) + 1};
  sdsl::int_vector<> runs(rlwt.runs_size, 0, rlwt.level_size);
  sdsl::int_vector<> lengths(rlwt.runs_size, 0, length_width);
  file << static_cast<uint64_t>(rlwt.level_size) << "\n";
  file << rlwt.runs_size << "\n";
  file << rlwt.alphabet_size << "\n";
  for (uint64_t runs_offset {}; runs_offset != rlwt.runs_size; ++runs_offset)
  {
    runs[runs_offset] = Access(rlwt, rlwt.first_length_bits_select(runs_offset + 1));
  }
  Print(file, runs);
  for (uint64_t runs_offset {}; runs_offset != rlwt.runs_size; ++runs_offset)
  {
    lengths[runs_offset] = rlwt.first_length_bits_select(runs_offset + 2) - rlwt.first_length_bits_select(runs_offset + 1);
  }
  Print(file, lengths);
  {
    sdsl::int_vector<> one_runs(rlwt.runs_size, 0, rlwt.level_size);
    sdsl::int_vector<> one_lengths(rlwt.runs_size, 0, length_width);
    uint64_t prefix_mask {};
    uint64_t mask {1ULL << (rlwt.level_size - 1)};
    auto run_bits_it {std::begin(rlwt.all_run_bits)};
    for (uint8_t level {}; level != rlwt.level_size; ++level)
    {
      auto runs_it {std::begin(runs)};
      auto runs_end {std::end(runs)};
      auto new_runs_it {std::begin(runs)};
      auto one_runs_it {std::begin(one_runs)};
      auto one_runs_end {std::end(one_runs)};
      auto lengths_it {std::begin(lengths)};
      auto new_lengths_it {std::begin(lengths)};
      auto one_lengths_it {std::begin(one_lengths)};
      auto prefix {*runs_it & prefix_mask};
      while (runs_it != runs_end)
      {
        while ((runs_it != runs_end) && ((*runs_it & prefix_mask) == prefix))
        {
          if (*run_bits_it++)
          {
            *one_runs_it++ = *runs_it++;
            *one_lengths_it++ = *lengths_it++;
          }
          else
          {
            *new_runs_it++ = *runs_it++;
            *new_lengths_it++ = *lengths_it++;
          }
        }
        {
          one_runs_end = one_runs_it;
          one_runs_it = std::begin(one_runs);
          one_lengths_it = std::begin(one_lengths);
          while (one_runs_it != one_runs_end)
          {
            *new_runs_it++ = *one_runs_it++;
            *new_lengths_it++ = *one_lengths_it++;
          }
          one_runs_it = std::begin(one_runs);
          one_lengths_it = std::begin(one_lengths);
        }
        prefix = *runs_it & prefix_mask;
      }
      Print
      (
        file,
        std::next(std::begin(rlwt.all_run_bits), level * rlwt.runs_size),
        std::next(std::begin(rlwt.all_run_bits), (level + 1) * rlwt.runs_size)
      );
      if (level == 0)
      {
        Print(file, rlwt.first_length_bits);
      }
      else
      {
        Print
        (
          file,
          std::next(std::begin(rlwt.middle_length_bits), (level - 1) * std::size(rlwt.first_length_bits)),
          std::next(std::begin(rlwt.middle_length_bits), level * std::size(rlwt.first_length_bits))
        );
      }
      Print(file, runs);
      Print(file, lengths);
      prefix_mask |= mask;
      mask >>= 1;
    }
    Print(file, rlwt.last_length_bits);
    Print(file, rlwt.run_bucket_begin_offsets);
    for (uint64_t run {}; run != std::size(rlwt.run_bucket_begin_offsets); ++run)
    {
      file << rlwt.last_length_bits_select(rlwt.run_bucket_begin_offsets[run] + 1);
      file << ((run != std::size(rlwt.run_bucket_begin_offsets) - 1) ? " " : "\n");
    }
  }
  return;
}

template
<
  typename Text,
  typename Runs,
  typename Lengths
>
void CalculateRunsAndLengths
(
  Text const &text,
  Runs &runs,
  Lengths &lengths
)
{
  auto it {std::begin(text)};
  auto runs_it {std::prev(std::begin(runs))};
  auto lengths_it {std::prev(std::begin(lengths))};
  auto prev_character {std::size(text)};
  while (it != std::end(text))
  {
    if (*it != prev_character)
    {
      ++runs_it;
      ++lengths_it;
      *runs_it = prev_character = *it;
      *lengths_it = 1;
    }
    else
    {
      ++(*lengths_it);
    }
    ++it;
  }
  return;
}

void ConstructRunLengthWaveletTree
(
  sdsl::int_vector<> const &text,
  RunLengthWaveletTree<> &rlwt
)
{
  // Print(std::cout, text);
  auto length_width {sdsl::bits::hi(std::size(text)) + 1};
  // std::cout << length_width << "\n";
  {
    rlwt.level_size = text.width();
    // std::cout << rlwt.level_size << "\n";
    rlwt.runs_size = CalculateRunsSize(text);
    // std::cout << rlwt.runs_size << "\n";
    rlwt.all_run_bits.resize(rlwt.level_size * rlwt.runs_size);
  }
  sdsl::int_vector<> runs(rlwt.runs_size, 0, rlwt.level_size);
  sdsl::int_vector<> lengths(rlwt.runs_size, 0, length_width);
  CalculateRunsAndLengths(text, runs, lengths);
  rlwt.alphabet_size = (*max_element(std::begin(runs), std::end(runs)) + 1);
  // std::cout << rlwt.alphabet_size << "\n";
  // Print(std::cout, runs);
  // Print(std::cout, lengths);
  {
    sdsl::bit_vector first_length_bits(std::size(text) + 1);
    sdsl::bit_vector middle_length_bits((rlwt.level_size - 1) * (std::size(text) + 1));
    sdsl::int_vector<> one_runs(rlwt.runs_size, 0, rlwt.level_size);
    sdsl::int_vector<> one_lengths(rlwt.runs_size, 0, length_width);
    uint64_t prefix_mask {};
    uint64_t mask {1ULL << (rlwt.level_size - 1)};
    auto run_bits_it {std::begin(rlwt.all_run_bits)};
    auto length_bits_it {std::begin(first_length_bits)};
    for (uint8_t level {}; level != rlwt.level_size; ++level)
    {
      if (level == 1)
      {
        length_bits_it = std::begin(middle_length_bits);
      }
      auto runs_it {std::begin(runs)};
      auto runs_end {std::end(runs)};
      auto new_runs_it {std::begin(runs)};
      auto one_runs_it {std::begin(one_runs)};
      auto one_runs_end {std::end(one_runs)};
      auto lengths_it {std::begin(lengths)};
      auto new_lengths_it {std::begin(lengths)};
      auto one_lengths_it {std::begin(one_lengths)};
      auto prefix {*runs_it & prefix_mask};
      while (runs_it != runs_end)
      {
        while ((runs_it != runs_end) && ((*runs_it & prefix_mask) == prefix))
        {
          *run_bits_it = *runs_it & mask;
          *length_bits_it = 1;
          length_bits_it += *lengths_it;
          if (*run_bits_it++)
          {
            *one_runs_it++ = *runs_it++;
            *one_lengths_it++ = *lengths_it++;
          }
          else
          {
            *new_runs_it++ = *runs_it++;
            *new_lengths_it++ = *lengths_it++;
          }
        }
        {
          one_runs_end = one_runs_it;
          one_runs_it = std::begin(one_runs);
          one_lengths_it = std::begin(one_lengths);
          while (one_runs_it != one_runs_end)
          {
            *new_runs_it++ = *one_runs_it++;
            *new_lengths_it++ = *one_lengths_it++;
          }
          one_runs_it = std::begin(one_runs);
          one_lengths_it = std::begin(one_lengths);
        }
        prefix = *runs_it & prefix_mask;
      }
      *length_bits_it++ = 1;
      // Print
      // (
      //   std::cout,
      //   std::next(std::begin(rlwt.all_run_bits), level * rlwt.runs_size),
      //   std::next(std::begin(rlwt.all_run_bits), (level + 1) * rlwt.runs_size)
      // );
      if (level == 0)
      {
        rlwt.first_length_bits = decltype(rlwt.first_length_bits)(first_length_bits);
        rlwt.first_length_bits_rank = decltype(rlwt.first_length_bits_rank)(&rlwt.first_length_bits);
        rlwt.first_length_bits_select = decltype(rlwt.first_length_bits_select)(&rlwt.first_length_bits);
        // Print(std::cout, first_length_bits);
      }
      // else
      // {
      //   Print
      //   (
      //     std::cout,
      //     std::next(std::begin(middle_length_bits), (level - 1) * (std::size(text) + 1)),
      //     std::next(std::begin(middle_length_bits), level * (std::size(text) + 1))
      //   );
      // }
      // Print(std::cout, runs);
      // Print(std::cout, lengths);
      prefix_mask |= mask;
      mask >>= 1;
    }
    rlwt.all_run_bits_rank = decltype(rlwt.all_run_bits_rank)(&rlwt.all_run_bits);
    rlwt.middle_length_bits = decltype(rlwt.middle_length_bits)(middle_length_bits);
    rlwt.middle_length_bits_select = decltype(rlwt.middle_length_bits_select)(&rlwt.middle_length_bits);
  }
  {
    sdsl::bit_vector last_length_bits(std::size(text) + 1);
    rlwt.run_bucket_begin_offsets.width(sdsl::bits::hi(rlwt.runs_size) + 1);
    rlwt.run_bucket_begin_offsets.resize(rlwt.alphabet_size + 1);
    auto runs_begin {std::begin(runs)};
    auto runs_it {std::begin(runs)};
    auto lengths_it {std::begin(lengths)};
    auto run_bucket_it {std::begin(rlwt.run_bucket_begin_offsets)};
    auto length_bits_it {std::begin(last_length_bits)};
    uint64_t prev_run {std::numeric_limits<uint64_t>::max()};
    while (runs_it != std::end(runs))
    {
      if (*runs_it != prev_run)
      {
        prev_run = *runs_it;
        *run_bucket_it++ = std::distance(runs_begin, runs_it);
      }
      ++runs_it;
      *length_bits_it = 1;
      length_bits_it += *lengths_it++;
    }
    *length_bits_it = 1;
    *(std::prev(std::end(rlwt.run_bucket_begin_offsets))) = rlwt.runs_size;
    rlwt.last_length_bits = decltype(rlwt.last_length_bits)(last_length_bits);
    rlwt.last_length_bits_select = decltype(rlwt.last_length_bits_select)(&rlwt.last_length_bits);
    // Print(std::cout, rlwt.last_length_bits);
    // Print(std::cout, rlwt.run_bucket_begin_offsets);
    // for (uint64_t run {}; run != (rlwt.alphabet_size + 1); ++run)
    // {
    //   std::cout << rlwt.last_length_bits_select(rlwt.run_bucket_begin_offsets[run] + 1);
    //   std::cout << ((run != rlwt.alphabet_size) ? " " : "\n");
    // }
  }
  return;
}

template <typename File, typename Node = InformationNode<std::string, uint64_t>>
uint64_t SerializeRunLengthWaveletTree
(
  RunLengthWaveletTree<> const &rlwt,
  File &index_file,
  std::shared_ptr<Node> root = nullptr
)
{
  if (root == nullptr)
  {
    sdsl::write_member(rlwt.level_size, index_file);
    sdsl::write_member(rlwt.runs_size, index_file);
    sdsl::write_member(rlwt.alphabet_size, index_file);
    sdsl::serialize(rlwt.all_run_bits, index_file);
    sdsl::serialize(rlwt.all_run_bits_rank, index_file);
    sdsl::serialize(rlwt.first_length_bits, index_file);
    sdsl::serialize(rlwt.first_length_bits_rank, index_file);
    sdsl::serialize(rlwt.first_length_bits_select, index_file);
    sdsl::serialize(rlwt.middle_length_bits, index_file);
    sdsl::serialize(rlwt.middle_length_bits_select, index_file);
    sdsl::serialize(rlwt.last_length_bits, index_file);
    sdsl::serialize(rlwt.last_length_bits_select, index_file);
    sdsl::serialize(rlwt.run_bucket_begin_offsets, index_file);
  }
  else
  {
    {
      auto node {std::make_shared<Node>("level_size")};
      node->value = sdsl::write_member(rlwt.level_size, index_file);
      root->value += node->value;
      root->children.emplace_back(node);
    }
    {
      auto node {std::make_shared<Node>("runs_size")};
      node->value = sdsl::write_member(rlwt.runs_size, index_file);
      root->value += node->value;
      root->children.emplace_back(node);
    }
    {
      auto node {std::make_shared<Node>("alphabet_size")};
      node->value = sdsl::write_member(rlwt.alphabet_size, index_file);
      root->value += node->value;
      root->children.emplace_back(node);
    }
    {
      auto node {std::make_shared<Node>("all_run_bits")};
      node->value = sdsl::serialize(rlwt.all_run_bits, index_file);
      root->value += node->value;
      root->children.emplace_back(node);
    }
    {
      auto node {std::make_shared<Node>("all_run_bits_rank")};
      node->value = sdsl::serialize(rlwt.all_run_bits_rank, index_file);
      root->value += node->value;
      root->children.emplace_back(node);
    }
    {
      auto node {std::make_shared<Node>("first_length_bits")};
      node->value = sdsl::serialize(rlwt.first_length_bits, index_file);
      root->value += node->value;
      root->children.emplace_back(node);
    }
    {
      auto node {std::make_shared<Node>("first_length_bits_rank")};
      node->value = sdsl::serialize(rlwt.first_length_bits_rank, index_file);
      root->value += node->value;
      root->children.emplace_back(node);
    }
    {
      auto node {std::make_shared<Node>("first_length_bits_select")};
      node->value = sdsl::serialize(rlwt.first_length_bits_select, index_file);
      root->value += node->value;
      root->children.emplace_back(node);
    }
    {
      auto node {std::make_shared<Node>("middle_length_bits")};
      node->value = sdsl::serialize(rlwt.middle_length_bits, index_file);
      // root->value += node->value;
      // root->children.emplace_back(node);
    }
    {
      auto node {std::make_shared<Node>("middle_length_bits_select")};
      node->value = sdsl::serialize(rlwt.middle_length_bits_select, index_file);
      root->value += node->value;
      root->children.emplace_back(node);
    }
    {
      auto node {std::make_shared<Node>("last_length_bits")};
      node->value = sdsl::serialize(rlwt.last_length_bits, index_file);
      root->value += node->value;
      root->children.emplace_back(node);
    }
    {
      auto node {std::make_shared<Node>("last_length_bits_select")};
      node->value = sdsl::serialize(rlwt.last_length_bits_select, index_file);
      root->value += node->value;
      root->children.emplace_back(node);
    }
    {
      auto node {std::make_shared<Node>("run_bucket_begin_offsets")};
      node->value = sdsl::serialize(rlwt.run_bucket_begin_offsets, index_file);
      root->value += node->value;
      root->children.emplace_back(node);
    }
    return root->value;
  }
  return 0;
}

template <typename File>
void LoadRunLengthWaveletTree
(
  RunLengthWaveletTree<> &rlwt,
  File &index_file
)
{
  sdsl::read_member(rlwt.level_size, index_file);
  sdsl::read_member(rlwt.runs_size, index_file);
  sdsl::read_member(rlwt.alphabet_size, index_file);
  sdsl::load(rlwt.all_run_bits, index_file);
  sdsl::load(rlwt.all_run_bits_rank, index_file);
  rlwt.all_run_bits_rank.set_vector(&rlwt.all_run_bits);
  sdsl::load(rlwt.first_length_bits, index_file);
  sdsl::load(rlwt.first_length_bits_rank, index_file);
  sdsl::load(rlwt.first_length_bits_select, index_file);
  rlwt.first_length_bits_rank.set_vector(&rlwt.first_length_bits);
  rlwt.first_length_bits_select.set_vector(&rlwt.first_length_bits);
  sdsl::load(rlwt.middle_length_bits, index_file);
  rlwt.middle_length_bits_select.set_vector(&rlwt.middle_length_bits);
  sdsl::load(rlwt.last_length_bits, index_file);
  rlwt.last_length_bits_select.set_vector(&rlwt.last_length_bits);
  sdsl::load(rlwt.run_bucket_begin_offsets, index_file);
  return;
}

struct Index
{
  sdsl::int_vector<8> grammar_rules;
  StaticGrammarTrie lex_grammar_count_trie;
  StaticGrammarTrie lex_grammar_rank_trie;
  StaticGrammarTrie colex_grammar_rank_trie;
  sdsl::int_vector<> colex_to_lex;
  sdsl::int_vector<> lex_rank_bucket_begin_offsets;
  RunLengthWaveletTree<> colex_bwt_rlwt;
};

template <typename File>
void PrintIndex
(
  File &file,
  Index &index
)
{
  Print(file, index.grammar_rules);
  PrintStaticGrammarTrie(file, index.grammar_rules, index.lex_grammar_count_trie, false);
  PrintStaticGrammarTrie(file, index.grammar_rules, index.lex_grammar_rank_trie);
  PrintStaticGrammarTrie(file, index.grammar_rules, index.colex_grammar_rank_trie);
  Print(file, index.colex_to_lex);
  Print(file, index.lex_rank_bucket_begin_offsets);
  PrintRunLengthWaveletTree(file, index.colex_bwt_rlwt);
  return;
}

void ConstructIndex
(
  Index &index,
  std::filesystem::path const &text_path
)
{
  sdsl::int_vector<8> text;
  {
    std::ifstream text_file {text_path};
    sdsl::load_vector_from_file(text, text_path);
    if (*std::prev(std::end(text)) != 0)
    {
      sdsl::append_zero_symbol(text);
      // std::cout << "zero symbol appended\n";
    }
    // Print(std::cout, text);
  }
  sdsl::bit_vector sl_types;
  {
    CalculateSlTypes(text, sl_types);
    // Print(std::cout, sl_types);
  }
  auto invalid_text_offset {std::size(text)};
  auto text_size_width {sdsl::bits::hi(std::size(text)) + 1};
  sdsl::int_vector<> text_offsets;
  {
    text_offsets.width(text_size_width);
    text_offsets.resize(std::size(text));
    sdsl::util::set_to_value(text_offsets, invalid_text_offset);
    // Print(std::cout, text_offsets);
  }
  sdsl::int_vector<> character_bucket_offsets;
  {
    character_bucket_offsets.width(text_size_width);
    character_bucket_offsets.resize(256);
  }
  {
    CalculateCharacterBucketBeginOffsets(text, character_bucket_offsets);
    // Print(std::cout, character_bucket_offsets);
    BucketSortRightmostLTypeCharacters(text, sl_types, character_bucket_offsets, text_offsets);
    // Print(std::cout, text_offsets);
    InduceSortLTypeCharacters(text, sl_types, character_bucket_offsets, text_offsets);
    // Print(std::cout, text_offsets);
    CalculateCharacterBucketEndOffsets(text, character_bucket_offsets);
    // Print(std::cout, character_bucket_offsets);
    InduceSortSTypeCharacters(text, sl_types, character_bucket_offsets, text_offsets);
    // Print(std::cout, text_offsets);
  }
  auto text_offsets_boundary {std::begin(text_offsets)};
  {
    text_offsets_boundary = MoveVaildEntriesToFront
    (
      std::begin(text_offsets),
      std::end(text_offsets),
      invalid_text_offset
    );
    // Print(std::cout, text_offsets);
  }
  auto grammar_rule_begin_offsets_begin {std::begin(text_offsets)};
  auto grammar_rule_begin_offsets_end {text_offsets_boundary};
  auto temporary_lex_text_begin {text_offsets_boundary};
  auto temporary_lex_text_end {std::end(text_offsets)};
  sdsl::int_vector<> grammar_rule_counts;
  {
    CalculateGrammarRuleCountsBeginOffsetsAndTemporaryLexText
    (
      text,
      sl_types,
      grammar_rule_counts,
      grammar_rule_begin_offsets_begin,
      grammar_rule_begin_offsets_end,
      temporary_lex_text_begin
    );
    // Print(std::cout, grammar_rule_counts);
    // Print(std::cout, grammar_rule_begin_offsets_begin, grammar_rule_begin_offsets_end);
    // Print(std::cout, temporary_lex_text_begin, temporary_lex_text_end);
  }
  sdsl::int_vector<> grammar_rule_sizes;
  {
    CalculateGrammarRuleSizes
    (
      grammar_rule_sizes,
      sl_types,
      grammar_rule_begin_offsets_begin,
      grammar_rule_begin_offsets_end,
      invalid_text_offset
    );
    // Print(std::cout, grammar_rule_sizes);
  }
  {
    CalculateGrammarRules
    (
      text,
      grammar_rule_sizes,
      index.grammar_rules,
      grammar_rule_begin_offsets_begin
    );
    // Print(std::cout, index.grammar_rules);
  }
  sdsl::int_vector<> lex_text;
  sdsl::int_vector<> lex_to_colex;
  auto grammar_ranks_size {std::size(grammar_rule_sizes) + 1};
  auto lex_text_width {sdsl::bits::hi(std::size(grammar_rule_sizes)) + 1};
  {
    CalculateLexText
    (
      lex_text,
      lex_text_width,
      temporary_lex_text_begin,
      temporary_lex_text_end,
      invalid_text_offset
    );
    // Print(std::cout, lex_text);
  }
  {
    sdsl::util::clear(sl_types);
    sdsl::util::clear(text);
    sdsl::util::clear(text_offsets);
  }
  {
    DynamicGrammarTrie lex_grammar_count_trie;
    InsertGrammarRuleSuffixesAndCounts
    (
      grammar_rule_sizes,
      index.grammar_rules,
      grammar_rule_counts,
      lex_grammar_count_trie
    );
    // PrintDynamicGrammarTrie(std::cout, index.grammar_rules, lex_grammar_count_trie, false);
    CalculateCumulativeGrammarCount(lex_grammar_count_trie);
    // PrintDynamicGrammarTrie(std::cout, index.grammar_rules, lex_grammar_count_trie, false);
    ConstructStaticGrammarTrie
    (
      lex_grammar_count_trie,
      index.lex_grammar_count_trie,
      false
    );
    // PrintStaticGrammarTrie(std::cout, index.grammar_rules, index.lex_grammar_count_trie, false);
  }
  {
    DynamicGrammarTrie lex_grammar_rank_trie;
    DynamicGrammarTrie colex_grammar_rank_trie(-1);
    InsertGrammarRules
    (
      grammar_rule_sizes,
      index.grammar_rules,
      lex_grammar_rank_trie,
      colex_grammar_rank_trie
    );
    // PrintDynamicGrammarTrie(std::cout, index.grammar_rules, lex_grammar_rank_trie);
    // PrintDynamicGrammarTrie(std::cout, index.grammar_rules, colex_grammar_rank_trie);
    CalculateCumulativeLexRankRanges(lex_grammar_rank_trie);
    // PrintDynamicGrammarTrie(std::cout, index.grammar_rules, lex_grammar_rank_trie);
    ConstructStaticGrammarTrie
    (
      lex_grammar_rank_trie,
      index.lex_grammar_rank_trie
    );
    // PrintStaticGrammarTrie(std::cout, index.grammar_rules, index.lex_grammar_rank_trie);
    lex_to_colex.width(lex_text_width);
    lex_to_colex.resize(grammar_ranks_size);
    lex_to_colex[0] = 0;
    CalculateCumulativeColexRankRangesAndLexToColex(colex_grammar_rank_trie, lex_to_colex);
    // Print(std::cout, lex_to_colex);
    // PrintDynamicGrammarTrie(std::cout, index.grammar_rules, colex_grammar_rank_trie);
    ConstructStaticGrammarTrie
    (
      colex_grammar_rank_trie,
      index.colex_grammar_rank_trie
    );
    // PrintStaticGrammarTrie(std::cout, index.grammar_rules, index.colex_grammar_rank_trie);
  }
  {
    index.colex_to_lex.width(lex_text_width);
    index.colex_to_lex.resize(grammar_ranks_size);
    for (uint64_t rank {}; rank != std::size(lex_to_colex); ++rank)
    {
      index.colex_to_lex[lex_to_colex[rank]] = rank;
    }
    // Print(std::cout, index.colex_to_lex);
  }
  {
    index.lex_rank_bucket_begin_offsets.width(sdsl::bits::hi(std::size(lex_text)) + 1);
    index.lex_rank_bucket_begin_offsets.resize(grammar_ranks_size + 1);
    CalculateCharacterBucketBeginOffsets(lex_text, index.lex_rank_bucket_begin_offsets);
    // Print(std::cout, index.lex_rank_bucket_begin_offsets);
  }
  {
    sdsl::int_vector<> colex_bwt;
    CalculateColexBwt(lex_text, lex_to_colex, colex_bwt);
    // Print(std::cout, colex_bwt);
    ConstructRunLengthWaveletTree(colex_bwt, index.colex_bwt_rlwt);
    // PrintRunLengthWaveletTree(std::cout, index.colex_bwt_rlwt);
  }
  return;
}

template <typename Node = InformationNode<std::string, uint64_t>>
uint64_t SerializeIndex
(
  Index const &index,
  std::filesystem::path const &index_path,
  std::shared_ptr<Node> root = nullptr
)
{
  std::fstream index_file(index_path, std::ios_base::out | std::ios_base::trunc);
  if (root == nullptr)
  {
    sdsl::serialize(index.grammar_rules, index_file);
    SerializeStaticGrammarTrie(index.lex_grammar_count_trie, index_file);
    SerializeStaticGrammarTrie(index.lex_grammar_rank_trie, index_file);
    SerializeStaticGrammarTrie(index.colex_grammar_rank_trie, index_file);
    sdsl::serialize(index.colex_to_lex, index_file);
    sdsl::serialize(index.lex_rank_bucket_begin_offsets, index_file);
    SerializeRunLengthWaveletTree(index.colex_bwt_rlwt, index_file);
  }
  else
  {
    {
      auto node {std::make_shared<Node>("grammar_rules")};
      node->value = sdsl::serialize(index.grammar_rules, index_file);
      root->value += node->value;
      root->children.emplace_back(node);
    }
    {
      auto node {std::make_shared<Node>("lex_grammar_count_trie")};
      SerializeStaticGrammarTrie(index.lex_grammar_count_trie, index_file, node);
      root->value += node->value;
      root->children.emplace_back(node);
    }
    {
      auto node {std::make_shared<Node>("lex_grammar_rank_trie")};
      SerializeStaticGrammarTrie(index.lex_grammar_rank_trie, index_file, node);
      root->value += node->value;
      root->children.emplace_back(node);
    }
    {
      auto node {std::make_shared<Node>("colex_grammar_rank_trie")};
      SerializeStaticGrammarTrie(index.colex_grammar_rank_trie, index_file, node);
      root->value += node->value;
      root->children.emplace_back(node);
    }
    {
      auto node {std::make_shared<Node>("colex_to_lex")};
      node->value = sdsl::serialize(index.colex_to_lex, index_file);
      root->value += node->value;
      root->children.emplace_back(node);
    }
    {
      auto node {std::make_shared<Node>("lex_rank_bucket_begin_offsets")};
      node->value = sdsl::serialize(index.lex_rank_bucket_begin_offsets, index_file);
      root->value += node->value;
      root->children.emplace_back(node);
    }
    {
      auto node {std::make_shared<Node>("colex_bwt_rlwt")};
      SerializeRunLengthWaveletTree(index.colex_bwt_rlwt, index_file, node);
      root->value += node->value;
      root->children.emplace_back(node);
    }
    return root->value;
  }
  return 0;
}

void LoadIndex
(
  Index &index,
  std::filesystem::path const &index_path
)
{
  std::ifstream index_file {index_path};
  index.grammar_rules.load(index_file);
  LoadStaticGrammarTrie(index.lex_grammar_count_trie, index_file);
  LoadStaticGrammarTrie(index.lex_grammar_rank_trie, index_file);
  LoadStaticGrammarTrie(index.colex_grammar_rank_trie, index_file);
  index.colex_to_lex.load(index_file);
  index.lex_rank_bucket_begin_offsets.load(index_file);
  LoadRunLengthWaveletTree(index.colex_bwt_rlwt, index_file);
  return;
}

template <typename Range>
constexpr bool IsEmptyRange (Range const &range)
{
  return (std::get<0>(range) > std::get<1>(range));
}

template <typename Range>
uint64_t CalculateRangeSize (Range const &range)
{
  if (std::get<0>(range) <= std::get<1>(range))
  {
    return (std::get<1>(range) - std::get<0>(range) + 1);
  }
  return 0;
}

template <typename TextIterator>
void CalculateSlFactor
(
  TextIterator const rend,
  TextIterator &rfirst,
  TextIterator &rlast
)
{
  uint64_t prev_sl_type {L};
  rfirst = rlast--;
  while ((rlast != rend) && !((prev_sl_type == S) && (*rlast > *std::next(rlast))))
  {
    if((prev_sl_type == L) && (*rlast < *std::next(rlast)))
    {
      prev_sl_type = S;
    }
    --rlast;
  }
  return;
}

template
<
  typename Labels,
  typename SlFactorIterator
>
uint64_t LookUpSlFactorCountInStaticGrammarTrie
(
  Labels const &labels,
  StaticGrammarTrie const &trie,
  SlFactorIterator it,
  SlFactorIterator last
)
{
  uint64_t begin_offset {};
  uint64_t end_offset {trie.level_order_select(1)};
  uint64_t offset {};
  while (begin_offset != end_offset)
  {
    while (begin_offset != end_offset)
    {
      offset = begin_offset + (end_offset - begin_offset) / 2;
      auto branch_character {labels[trie.edge_begin_offsets[offset]]};
      if (*it == branch_character)
      {
        break;
      }
      else if (*it < branch_character)
      {
        end_offset = offset;
      }
      else
      {
        begin_offset = offset + 1;
      }
    }
    if (begin_offset != end_offset)
    {
      auto edge_it {std::next(std::begin(labels), trie.edge_begin_offsets[offset])};
      auto edge_end {std::next(std::begin(labels), trie.edge_prev_end_offsets[offset] + trie.step)};
      while ((it != last) && (edge_it != edge_end) && (*it == *edge_it))
      {
        it += trie.step;
        edge_it += trie.step;
      }
      if (it != last)
      {
        if (edge_it != edge_end)
        {
          break;
        }
        else
        {
          begin_offset = trie.level_order_select(offset + 1) - (offset + 1) + 1;
          end_offset = trie.level_order_select(offset + 2) - (offset + 2) + 1;
        }
      }
      else
      {
        return trie.counts[offset];
      }
    }
  }
  return 0;
}

template
<
  typename Labels,
  typename SlFactorIterator
>
std::pair<uint64_t, uint64_t> LookUpSlFactorRankRangeInStaticGrammarTrie
(
  Labels const &labels,
  StaticGrammarTrie const &trie,
  SlFactorIterator it,
  SlFactorIterator last
)
{
  uint64_t begin_offset {};
  uint64_t end_offset {trie.level_order_select(1)};
  uint64_t offset {};
  while (begin_offset != end_offset)
  {
    while (begin_offset != end_offset)
    {
      offset = begin_offset + (end_offset - begin_offset) / 2;
      auto branch_character {labels[trie.edge_begin_offsets[offset]]};
      if (*it == branch_character)
      {
        break;
      }
      else if (*it < branch_character)
      {
        end_offset = offset;
      }
      else
      {
        begin_offset = offset + 1;
      }
    }
    if (begin_offset != end_offset)
    {
      auto edge_it {std::next(std::begin(labels), trie.edge_begin_offsets[offset])};
      auto edge_end {std::next(std::begin(labels), trie.edge_prev_end_offsets[offset] + trie.step)};
      while ((it != last) && (edge_it != edge_end) && (*it == *edge_it))
      {
        it += trie.step;
        edge_it += trie.step;
      }
      if (it != last)
      {
        if (edge_it != edge_end)
        {
          break;
        }
        else
        {
          begin_offset = trie.level_order_select(offset + 1) - (offset + 1) + 1;
          end_offset = trie.level_order_select(offset + 2) - (offset + 2) + 1;
        }
      }
      else
      {
        return
        {
          trie.leftmost_ranks[offset],
          trie.rightmost_ranks[offset]
        };
      }
    }
  }
  return {1, 0};
}

template
<
  typename PatternRange,
  typename PatternIterator
>
void BackwardSearchPatternPrefix
(
  Index const &index,
  PatternRange &pattern_range_l,
  PatternRange &pattern_range_s,
  PatternIterator rfirst,
  PatternIterator rlast
)
{
  auto colex_rank_range
  {
    LookUpSlFactorRankRangeInStaticGrammarTrie
    (
      index.grammar_rules,
      index.colex_grammar_rank_trie,
      rfirst,
      rlast
    )
  };
  if (!IsEmptyRange(colex_rank_range))
  {
    if (!IsEmptyRange(pattern_range_l))
    {
      pattern_range_l =
      {
        1,
        RangeCount
        (
          index.colex_bwt_rlwt,
          std::get<0>(pattern_range_l),
          std::get<1>(pattern_range_l),
          std::get<0>(colex_rank_range),
          std::get<1>(colex_rank_range)
        )
      };
    }
    if (!IsEmptyRange(pattern_range_s))
    {
      pattern_range_s =
      {
        1,
        RangeCount
        (
          index.colex_bwt_rlwt,
          std::get<0>(pattern_range_s),
          std::get<1>(pattern_range_s),
          std::get<0>(colex_rank_range),
          std::get<1>(colex_rank_range)
        )
      };
    }
  }
  // {
  //   Print(std::cout, rfirst, rlast, -1);
  //   std::cout
  //   << "->[" << std::get<0>(colex_rank_range) << "," << std::get<1>(colex_rank_range) << "]"
  //   << "->L:[" << std::get<0>(pattern_range_l) << "," << std::get<1>(pattern_range_l) << "]"
  //   << "->S:[" << std::get<0>(pattern_range_s) << "," << std::get<1>(pattern_range_s) << "]\n";
  // }
  return;
}

template
<
  typename PatternRange,
  typename PatternIterator
>
auto BackwardSearchPatternProperSubstring
(
  Index const &index,
  PatternRange &pattern_range_l,
  PatternRange &pattern_range_s,
  PatternIterator rfirst,
  PatternIterator rlast
)
{
  auto colex_rank_range
  {
    LookUpSlFactorRankRangeInStaticGrammarTrie
    (
      index.grammar_rules,
      index.colex_grammar_rank_trie,
      rfirst,
      rlast
    )
  };
  auto colex_rank {std::get<0>(colex_rank_range)};
  if (!IsEmptyRange(colex_rank_range))
  {
    auto begin_offset {index.lex_rank_bucket_begin_offsets[index.colex_to_lex[colex_rank]]};
    if (!IsEmptyRange(pattern_range_l))
    {
      pattern_range_l =
      {
        begin_offset + Rank(index.colex_bwt_rlwt, std::get<0>(pattern_range_l), colex_rank),
        begin_offset + Rank(index.colex_bwt_rlwt, std::get<1>(pattern_range_l) + 1, colex_rank) - 1
      };
    }
    if (!IsEmptyRange(pattern_range_s))
    {
      pattern_range_s =
      {
        begin_offset + Rank(index.colex_bwt_rlwt, std::get<0>(pattern_range_s), colex_rank),
        begin_offset + Rank(index.colex_bwt_rlwt, std::get<1>(pattern_range_s) + 1, colex_rank) - 1
      };
    }
  }
  else
  {
    pattern_range_l = pattern_range_s = {1, 0};
  }
  // {
  //   Print(std::cout, rfirst, rlast, -1);
  //   std::cout
  //   << "->[" << colex_rank << ":" << index.colex_to_lex[colex_rank] << "]"
  //   << "->L:[" << std::get<0>(pattern_range_l) << "," << std::get<1>(pattern_range_l) << "]"
  //   << "->S:[" << std::get<0>(pattern_range_s) << "," << std::get<1>(pattern_range_s) << "]\n";
  // }
  return rlast;
}

template <typename PatternIterator>
auto CalculatePatternSuffixS
(
  PatternIterator rbegin,
  PatternIterator rend
)
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

template
<
  typename PatternRange,
  typename PatternIterator
>
auto BackwardSearchPatternSuffix
(
  Index const &index,
  PatternRange &pattern_range_l,
  PatternRange &pattern_range_s,
  PatternIterator rbegin,
  PatternIterator rend
)
{
  auto rfirst {rbegin};
  auto rlast {CalculatePatternSuffixS(rbegin, rend)};
  if (rlast == rend)
  {
    std::get<1>(pattern_range_l) = LookUpSlFactorCountInStaticGrammarTrie
    (
      index.grammar_rules,
      index.lex_grammar_count_trie,
      std::next(rend),
      std::next(rbegin)
    );
    // {
    //   Print(std::cout, std::next(rend), std::next(rbegin));
    //   std::cout << "->L:[" << std::get<0>(pattern_range_l) << "," << std::get<1>(pattern_range_l) << "]\n";
    // }
  }
  else
  {
    if (rlast != rbegin)
    {
      auto lex_rank_range
      {
        LookUpSlFactorRankRangeInStaticGrammarTrie
        (
          index.grammar_rules,
          index.lex_grammar_rank_trie,
          std::next(rlast),
          std::next(rbegin)
        )
      };
      if (!IsEmptyRange(lex_rank_range))
      {
        pattern_range_s =
        {
          index.lex_rank_bucket_begin_offsets[std::get<0>(lex_rank_range)],
          (index.lex_rank_bucket_begin_offsets[std::get<1>(lex_rank_range) + 1] - 1)
        };
      }
      // {
      //   Print(std::cout, std::next(rlast), std::next(rbegin));
      //   std::cout << "->[" << std::get<0>(lex_rank_range) << "," << std::get<1>(lex_rank_range) << "]";
      //   std::cout << "->S:[" << std::get<0>(pattern_range_s) << "," << std::get<1>(pattern_range_s) << "]\n";
      // }
      CalculateSlFactor(rend, rfirst, rlast);
      if (!IsEmptyRange(pattern_range_s))
      {
        auto colex_rank_range
        {
          LookUpSlFactorRankRangeInStaticGrammarTrie
          (
            index.grammar_rules,
            index.colex_grammar_rank_trie,
            rfirst,
            rlast
          )
        };
        if (!IsEmptyRange(colex_rank_range))
        {
          if (rlast == rend)
          {
            pattern_range_s =
            {
              1,
              RangeCount
              (
                index.colex_bwt_rlwt,
                std::get<0>(pattern_range_s),
                std::get<1>(pattern_range_s),
                std::get<0>(colex_rank_range),
                std::get<1>(colex_rank_range)
              )
            };
          }
          else
          {
            auto colex_rank {std::get<0>(colex_rank_range)};
            auto begin_offset {index.lex_rank_bucket_begin_offsets[index.colex_to_lex[colex_rank]]};
            pattern_range_s =
            {
              begin_offset + Rank(index.colex_bwt_rlwt, std::get<0>(pattern_range_s), colex_rank),
              begin_offset + Rank(index.colex_bwt_rlwt, std::get<1>(pattern_range_s) + 1, colex_rank) - 1
            };
          }
        }
        else
        {
          pattern_range_s = {1, 0};
        }
        // {
        //   Print(std::cout, rfirst, rlast, -1);
        //   std::cout
        //   << "->[" << std::get<0>(colex_rank_range)
        //   << ":" << index.colex_to_lex[std::get<0>(colex_rank_range)]
        //   << "," << std::get<1>(colex_rank_range) << "]"
        //   << "->S:[" << std::get<0>(pattern_range_s) << "," << std::get<1>(pattern_range_s) << "]\n";
        // }
      }
    }
    else
    {
      CalculateSlFactor(rend, rfirst, rlast);
    }
    if (rlast == rend)
    {
      std::get<1>(pattern_range_l) = LookUpSlFactorCountInStaticGrammarTrie
      (
        index.grammar_rules,
        index.lex_grammar_count_trie,
        std::next(rend),
        std::next(rbegin)
      );
      // {
      //   Print(std::cout, std::next(rend), std::next(rbegin));
      //   std::cout << "->L:[" << std::get<0>(pattern_range_l) << "," << std::get<1>(pattern_range_l) << "]\n";
      // }
    }
    else
    {
      auto lex_rank_range
      {
        LookUpSlFactorRankRangeInStaticGrammarTrie
        (
          index.grammar_rules,
          index.lex_grammar_rank_trie,
          std::next(rlast),
          std::next(rbegin)
        )
      };
      if (!IsEmptyRange(lex_rank_range))
      {
        pattern_range_l =
        {
          index.lex_rank_bucket_begin_offsets[std::get<0>(lex_rank_range)],
          (index.lex_rank_bucket_begin_offsets[std::get<1>(lex_rank_range) + 1] - 1)
        };
      }
      // {
      //   Print(std::cout, std::next(rlast), std::next(rbegin));
      //   std::cout << "->[" << std::get<0>(lex_rank_range) << "," << std::get<1>(lex_rank_range) << "]";
      //   std::cout << "->L:[" << std::get<0>(pattern_range_l) << "," << std::get<1>(pattern_range_l) << "]\n";
      // }
    }
  }
  return rlast;
}

template <typename PatternIterator>
uint64_t Count
(
  Index const &index,
  PatternIterator begin,
  PatternIterator end
)
{
  std::pair<uint64_t, uint64_t> pattern_range_l {1, 0};
  std::pair<uint64_t, uint64_t> pattern_range_s {1, 0};
  auto rbegin {std::prev(end)};
  auto rend {std::prev(begin)};
  auto rfirst {rbegin};
  auto rlast {rbegin};
  // Print(std::cout, begin, end);
  rlast = BackwardSearchPatternSuffix(index, pattern_range_l, pattern_range_s, rbegin, rend);
  if (rlast != rend)
  {
    while (!IsEmptyRange(pattern_range_l) || !IsEmptyRange(pattern_range_s))
    {
      CalculateSlFactor(rend, rfirst, rlast);
      if (rlast != rend)
      {
        rlast = BackwardSearchPatternProperSubstring(index, pattern_range_l, pattern_range_s, rfirst, rlast);
      }
      else
      {
        BackwardSearchPatternPrefix(index, pattern_range_l, pattern_range_s, rfirst, rlast);
        break;
      }
    }
  }
  return (CalculateRangeSize(pattern_range_l) + CalculateRangeSize(pattern_range_s));
}
}
