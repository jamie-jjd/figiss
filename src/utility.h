#pragma once

#include <cmath>

#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <map>
#include <memory>
#include <random>
#include <string>
#include <unordered_map>
#include <unordered_set>

#include <sdsl/suffix_trees.hpp>

namespace project
{
std::filesystem::path CreateParentDirectoryByCategory
(
  std::string const &category,
  std::filesystem::path const &path
)
{
  std::filesystem::path parent_category_path
  {
    std::string{"../data/"}
    + category + "/"
    + path.parent_path().filename().string()
  };
  if (!std::filesystem::exists(parent_category_path))
  {
    std::filesystem::create_directories(parent_category_path);
  }
  return parent_category_path;
}

std::filesystem::path CreatePath
(
  std::filesystem::path const &parent_path,
  std::string const &filename,
  std::string const &extensions = ""
)
{
  return (parent_path / (filename + extensions));
}

void CalculatePrefix
(
  std::filesystem::path const &text_path,
  uint64_t const size_in_megabytes
)
{
  sdsl::int_vector<8> text;
  sdsl::load_vector_from_file(text, text_path);
  auto parent_prefix_path
  {
    std::string{"../data/corpus/"}
    + std::to_string(size_in_megabytes) + "mb"
  };
  if (!std::filesystem::exists(parent_prefix_path))
  {
    std::filesystem::create_directories(parent_prefix_path);
  }
  auto prefix_path
  {
    CreatePath
    (
      parent_prefix_path,
      text_path.filename().string(),
      std::string{"."} + std::to_string(size_in_megabytes) + "mb"
    )
  };
  std::ofstream prefix_file {prefix_path};
  for (uint64_t i {}; i != (size_in_megabytes * 1024 * 1024); ++i)
  {
    prefix_file << text[i];
  }
  return;
}

template <typename Patterns>
void GeneratePatterns
(
  std::filesystem::path const &text_path,
  uint64_t const pattern_amount,
  uint64_t const pattern_size,
  Patterns &patterns
)
{
  sdsl::int_vector<8> text;
  sdsl::load_vector_from_file(text, text_path);
  if (pattern_size < std::size(text))
  {
    std::mt19937 engine {std::random_device{}()};
    std::uniform_int_distribution<uint64_t> distribution(0, std::size(text) - pattern_size);
    auto random_begin_offset {std::bind(distribution, engine)};
    patterns.resize(pattern_amount * pattern_size);
    auto patterns_iterator {std::begin(patterns)};
    for (uint64_t i {0}; i != pattern_amount; ++i)
    {
      auto text_iterator {std::next(std::begin(text), random_begin_offset())};
      for (uint64_t j {0}; j != pattern_size; ++j)
      {
        *patterns_iterator = *text_iterator;
        ++patterns_iterator;
        ++text_iterator;
      }
    }
  }
  return;
}

auto GenerateAndSerializePatterns
(
  std::filesystem::path const &text_path,
  uint64_t const pattern_amount,
  uint64_t const pattern_size
)
{
  sdsl::int_vector<8> patterns;
  GeneratePatterns(text_path, pattern_amount, pattern_size, patterns);
  auto parent_patterns_path {CreateParentDirectoryByCategory("patterns", text_path)};
  auto patterns_path
  {
    CreatePath
    (
      parent_patterns_path,
      text_path.filename().string(),
      std::string{"."} + std::to_string(pattern_amount) + "." + std::to_string(pattern_size)
    )
  };
  if (std::size(patterns) == (pattern_amount * pattern_size))
  {
    std::ofstream patterns_file {patterns_path};
    sdsl::int_vector<> pattern_amount_size(2, 0, 64);
    pattern_amount_size[0] = pattern_amount;
    pattern_amount_size[1] = pattern_size;
    sdsl::serialize(pattern_amount_size, patterns_file);
    sdsl::serialize(patterns, patterns_file);
  }
  return patterns_path;
}

template <typename Patterns>
void LoadPatterns
(
  std::filesystem::path const &patterns_path,
  uint64_t &pattern_amount,
  uint64_t &pattern_size,
  Patterns &patterns
)
{
  std::ifstream patterns_file {patterns_file};
  sdsl::int_vector<> pattern_amount_size;
  pattern_amount_size.load(patterns_file);
  pattern_amount = pattern_amount_size[0];
  pattern_size = pattern_amount_size[1];
  patterns.load(patterns_file);
  return;
}

void ConvertByteToIntText
(
  std::string const &category,
  std::filesystem::path const &byte_text_path
)
{
  sdsl::int_vector<8> byte_text;
  sdsl::load_vector_from_file(byte_text, byte_text_path);
  auto int_text_size {std::size(byte_text)};
  if (*std::prev(std::end(byte_text)) != 0)
  {
    ++int_text_size;
  }
  sdsl::int_vector<> int_text(int_text_size, 0, 8);
  std::copy(std::begin(byte_text), std::end(byte_text), std::begin(int_text));
  auto parent_int_text_path {CreateParentDirectoryByCategory(category, byte_text_path)};
  auto int_text_path {CreatePath(parent_int_text_path, byte_text_path.filename().string())};
  std::ofstream int_text_file {int_text_path};
  sdsl::serialize(int_text, int_text_file);
  return;
}

void CalculateAndSerializeCompactText
(
  std::string const &category,
  std::filesystem::path const &text_path
)
{
  sdsl::int_vector<> text;
  sdsl::load_from_file(text, text_path);
  if (*std::prev(std::end(text)) != 0)
  {
    sdsl::append_zero_symbol(text);
  }
  uint64_t codebook_size {*std::max_element(std::begin(text), std::end(text)) + 1};
  sdsl::int_vector<> codebook(codebook_size, 0);
  {
    auto text_iterator {std::begin(text)};
    auto text_end {std::end(text)};
    while (text_iterator != text_end)
    {
      if (*text_iterator != 0)
      {
        codebook[*text_iterator] = 1;
      }
      ++text_iterator;
    }
  }
  std::partial_sum(std::begin(codebook), std::end(codebook), std::begin(codebook));
  {
    auto text_iterator {std::begin(text)};
    auto text_end {std::end(text)};
    while (text_iterator != text_end)
    {
      *text_iterator = codebook[*text_iterator];
      ++text_iterator;
    }
  }
  sdsl::util::bit_compress(text);
  auto parent_compact_text_path {CreateParentDirectoryByCategory(category, text_path)};
  auto compact_text_path {CreatePath(parent_compact_text_path, text_path.filename().string())};
  std::ofstream compact_text_file {compact_text_path};
  sdsl::serialize(text, compact_text_file);
  return;
}

void CalculateAndSerializeBwt
(
  std::string const &category,
  std::filesystem::path const &text_path
)
{
  sdsl::int_vector<> text;
  sdsl::load_from_file(text, text_path);
  sdsl::int_vector<> buffer;
  sdsl::qsufsort::construct_sa(buffer, text);
  auto buffer_iterator {std::begin(buffer)};
  auto buffer_end {std::end(buffer)};
  while (buffer_iterator != buffer_end)
  {
    if (*buffer_iterator != 0)
    {
      *buffer_iterator = text[*buffer_iterator - 1];
    }
    ++buffer_iterator;
  }
  sdsl::util::bit_compress(buffer);
  auto parent_bwt_path {CreateParentDirectoryByCategory(category, text_path)};
  auto bwt_path {CreatePath(parent_bwt_path, text_path.filename().string())};
  std::ofstream bwt_file {bwt_path};
  sdsl::serialize(buffer, bwt_file);
  return;
}

template <typename Text>
uint64_t CalculateRuns (Text const &text)
{
  uint64_t runs {};
  uint64_t previous_character {};
  auto iterator {std::begin(text)};
  while (iterator != std::end(text))
  {
    if (previous_character != *iterator)
    {
      ++runs;
      previous_character = *iterator;
    }
    ++iterator;
  }
  return runs;
}

template <typename Text>
double CalculateZerothEmpiricalEntropy (Text const &text)
{
  std::map<uint64_t, uint64_t> alphabet_count;
  auto text_iterator {std::begin(text)};
  while (text_iterator != std::end(text))
  {
    auto alphabet_count_iterator {alphabet_count.find(*text_iterator)};
    if (alphabet_count_iterator != std::end(alphabet_count))
    {
      ++std::get<1>(*alphabet_count_iterator);
    }
    else
    {
      alphabet_count[*text_iterator] = 1;
    }
    ++text_iterator;
  }
  double zeroth_empirical_entropy {};
  auto lg_size {std::log2(std::size(text))};
  for (auto const &character_count : alphabet_count)
  {
    auto count {static_cast<double>(std::get<1>(character_count))};
    zeroth_empirical_entropy += (count * (lg_size - std::log2(count)));
  }
  return (zeroth_empirical_entropy / std::size(text));
}

template <typename Text>
struct KmerTrie
{
  struct Node
  {
    uint64_t edge_begin_offset;
    uint64_t edge_end_offset;
    std::map<uint64_t, std::shared_ptr<Node>> branches;
    std::vector<uint64_t> context;

    Node () = default;
    Node
    (
      uint64_t const begin_offset,
      uint64_t const end_offset
    )
    : edge_begin_offset {begin_offset},
      edge_end_offset {end_offset},
      branches {},
      context {}
    {
    }
  };

  std::shared_ptr<Node> root;
  double kth_empirical_entropy;

  KmerTrie
  (
    Text const &text,
    uint64_t const k
  )
  : root {std::make_shared<Node>()},
    kth_empirical_entropy {}
  {
    auto text_begin {std::begin(text)};
    auto text_iterator {text_begin};
    auto text_last_iterator {std::prev(std::end(text), k)};
    while (text_iterator != text_last_iterator)
    {
      auto k_mer_begin {text_iterator};
      auto k_mer_end {std::next(k_mer_begin, k)};
      Insert(text_begin, k_mer_begin, k_mer_end);
      ++text_iterator;
    }
    kth_empirical_entropy = CalculateCumulative0thEmpiricalEntropy(root) / std::size(text);
  }

  template <typename TextIterator>
  void Insert
  (
    TextIterator text_begin,
    TextIterator k_mer_iterator,
    TextIterator k_mer_end
  )
  {
    auto current_node {root};
    while (k_mer_iterator != k_mer_end)
    {
      auto character {*k_mer_iterator};
      auto branches_iterator {current_node->branches.find(character)};
      if (branches_iterator == std::end(current_node->branches))
      {
        current_node = current_node->branches[character] = std::make_shared<Node>
        (
          std::distance(text_begin, k_mer_iterator),
          std::distance(text_begin, k_mer_end)
        );
        k_mer_iterator = k_mer_end;
      }
      else
      {
        auto child_node {std::get<1>(*branches_iterator)};
        auto edge_iterator {std::next(text_begin, child_node->edge_begin_offset)};
        auto edge_end {std::next(text_begin, child_node->edge_end_offset)};
        while
        (
          (k_mer_iterator != k_mer_end)
          &&
          (edge_iterator != edge_end)
          &&
          (*k_mer_iterator == *edge_iterator)
        )
        {
          ++k_mer_iterator;
          ++edge_iterator;
        }
        if
        (
          (k_mer_iterator != k_mer_end)
          &&
          (edge_iterator != edge_end)
        )
        {
          auto edge_it_offset {std::distance(text_begin, edge_iterator)};
          auto internal_node {std::make_shared<Node>(child_node->edge_begin_offset, edge_it_offset)};
          child_node->edge_begin_offset = edge_it_offset;
          internal_node->branches[*edge_iterator] = child_node;
          current_node = internal_node->branches[*k_mer_iterator] = std::make_shared<Node>
          (
            std::distance(text_begin, k_mer_iterator),
            std::distance(text_begin, k_mer_end)
          );
          std::get<1>(*branches_iterator) = internal_node;
        }
        else
        {
          current_node = child_node;
        }
      }
    }
    current_node->context.push_back(*k_mer_end);
    return;
  }

  double CalculateCumulative0thEmpiricalEntropy (std::shared_ptr<Node> current_node)
  {
    double cumulative_0th_empirical_entropy {};
    if (std::size(current_node->context) == 0)
    {
      auto branches_iterator {std::begin(current_node->branches)};
      auto branches_end {std::end(current_node->branches)};
      while (branches_iterator != branches_end)
      {
        cumulative_0th_empirical_entropy += CalculateCumulative0thEmpiricalEntropy(std::get<1>(*branches_iterator));
        ++branches_iterator;
      }
    }
    else
    {
      cumulative_0th_empirical_entropy +=
      (CalculateZerothEmpiricalEntropy(current_node->context) * std::size(current_node->context));
    }
    return cumulative_0th_empirical_entropy;
  }
};

uint64_t CalculateCompressedBitSize
(
  std::string const &command_prefix,
  std::filesystem::path const &text_path,
  std::string const &extension
)
{
  auto command {command_prefix + text_path.string()};
  uint64_t bit_size {};
  if (std::system(command.c_str()) == 0)
  {
    auto compressed_text_path {std::filesystem::path{text_path.string() + "." + extension}};
    bit_size = (std::filesystem::file_size(compressed_text_path) * 8);
    std::filesystem::remove(compressed_text_path);
  }
  return bit_size;
}

void PrintTextStatistics
(
  std::string const &category,
  std::filesystem::path const &text_path
)
{
  sdsl::int_vector<> text;
  sdsl::load_from_file(text, text_path);
  auto parent_statistics_path {CreateParentDirectoryByCategory(category, text_path)};
  auto statistics_path {CreatePath(parent_statistics_path, text_path.filename().string())};
  std::ofstream statistics_file {statistics_path};
  auto n {std::size(text)};
  auto r {CalculateRuns(text)};
  auto lg_sigma {static_cast<uint64_t>(text.width())};
  auto H_0 {CalculateZerothEmpiricalEntropy(text)};
  std::vector<double> H_ks;
  for (uint64_t power {0}; power != 4; ++power)
  {
    KmerTrie k_mer_trie(text, (1ULL << power));
    H_ks.push_back(k_mer_trie.kth_empirical_entropy);
  }
  auto bzip2_bit_size {CalculateCompressedBitSize("bzip2 -k -9 -v ", text_path, "bz2")};
  auto p7zip_bit_size {CalculateCompressedBitSize("p7zip -k ", text_path, "7z")};
  auto repair_bit_size {CalculateCompressedBitSize("repair ", text_path, "rp")};
  statistics_file << std::fixed << std::setprecision(2);
  statistics_file << text_path.filename().string() << "\n";
  statistics_file << (n / (1024.0 * 1024.0)) << "M & ";
  statistics_file << (r / (1024.0 * 1024.0)) << "M & ";
  statistics_file << lg_sigma << " & ";
  statistics_file << H_0 << " & ";
  for (auto const &H_k : H_ks)
  {
    statistics_file << H_k << " & ";
  }
  statistics_file << (bzip2_bit_size / (1024.0 * 1024.0 * 8.0)) << "MB & ";
  statistics_file << (p7zip_bit_size / (1024.0 * 1024.0 * 8.0)) << "MB & ";
  statistics_file << (repair_bit_size / (1024.0 * 1024.0 * 8.0)) << "MB\n";
  statistics_file << " & ";
  statistics_file << "(" << (static_cast<double>(r) / n) * 100.0 << "\\%) & ";
  statistics_file << " & ";
  statistics_file << "(" << (H_0 / lg_sigma) * 100.0 << "\\%) & ";
  for (auto const &H_k : H_ks)
  {
    statistics_file << "(" << (H_k / lg_sigma) * 100.0 << "\\%) & ";
  }
  statistics_file << "(" << (static_cast<double>(bzip2_bit_size) / text.bit_size()) * 100.0 << "\\%) & ";
  statistics_file << "(" << (static_cast<double>(p7zip_bit_size) / text.bit_size()) * 100.0 << "\\%) & ";
  statistics_file << "(" << (static_cast<double>(repair_bit_size) / text.bit_size()) * 100.0 << "\\%)\n";
  return;
}
}
