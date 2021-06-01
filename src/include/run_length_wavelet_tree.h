#pragma once

#include <sdsl/bit_vectors.hpp>

#include "utility.h"

namespace project
{
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
}
