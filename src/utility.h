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

void GeneratePrefix
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
  std::filesystem::path const &byte_text_path,
  std::string const &category
)
{
  sdsl::int_vector<8> byte_text;
  sdsl::load_vector_from_file(byte_text, byte_text_path);
  sdsl::int_vector<> int_text(std::size(byte_text), 0, 8);
  std::copy(std::begin(byte_text), std::end(byte_text), std::begin(int_text));
  auto parent_int_text_path {CreateParentDirectoryByCategory(category, byte_text_path)};
  auto int_text_path {CreatePath(parent_int_text_path, byte_text_path.filename().string())};
  std::ofstream int_text_file {int_text_path};
  sdsl::serialize(int_text, int_text_file);
  return;
}

void GenerateAndSerializeCompactText
(
  std::filesystem::path const &text_path,
  std::string const &category
)
{
  sdsl::int_vector<> text;
  sdsl::load_from_file(text, text_path);
  uint64_t codebook_size {*std::max_element(std::begin(text), std::end(text)) + 1};
  sdsl::int_vector<> codebook(codebook_size, 0);
  {
    auto text_iterator {std::begin(text)};
    auto text_end {std::end(text)};
    while (text_iterator != text_end)
    {
      codebook[*text_iterator] = 1;
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

template <typename CategoriesCountMap>
double CalculateZerothEmpiricalEntropy
(
  CategoriesCountMap const &categories_count,
  uint64_t size
)
{
  double zeroth_empirical_entropy {};
  auto log_size {std::log2(size)};
  for (auto const &category_count : categories_count)
  {
    auto count {static_cast<double>(std::get<1>(category_count))};
    zeroth_empirical_entropy += (count * (log_size - std::log2(count)));
  }
  return (zeroth_empirical_entropy / size);
}

template
<
  typename File,
  typename Text
>
void PrintTextSize
(
  File &file,
  Text const &text
)
{
  file
  << std::fixed
  << std::setprecision(2)
  << "size,"
  << std::size(text)
  << ",bit_width,"
  << static_cast<uint64_t>(text.width())
  << ",bit_size,"
  << (static_cast<uint64_t>(text.bit_size()) / (1024.0 * 1024.0 * 8.0)) << "(MB)\n";
  return;
}

template
<
  typename File,
  typename Text
>
void PrintAlphabetSize
(
  File &file,
  Text const &text
)
{
  std::unordered_set<uint64_t> alphabet;
  for (auto const &character : text)
  {
    alphabet.insert(character);
  }
  auto bit_width {static_cast<double>(sdsl::bits::hi(alphabet.size()) + 1)};
  file
  << std::fixed
  << std::setprecision(2)
  << "alphabet_size,"
  << std::size(alphabet)
  << ",compression_ratio,"
  << (bit_width / text.width()) * 100.0
  << "%\n";
  return;
}

template
<
  typename File,
  typename Text
>
void Print0thEmpiricalEntropy
(
  File &file,
  Text const &text
)
{
  std::map<uint64_t, uint64_t> alphabet_count;
  for (auto const &character : text)
  {
    auto iterator {alphabet_count.find(character)};
    if (iterator != std::end(alphabet_count))
    {
      ++std::get<1>(*iterator);
    }
    else
    {
      alphabet_count[character] = 1;
    }
  }
  auto zeroth_empirical_entropy {CalculateZerothEmpiricalEntropy(alphabet_count, std::size(text))};
  file
  << std::fixed
  << std::setprecision(2)
  << "nH_0,"
  << (zeroth_empirical_entropy * std::size(text) / (1024.0 * 1024.0 * 8.0)) << "(MB),"
  << "compression_ratio,"
  << (zeroth_empirical_entropy / text.width()) * 100.0
  << "%\n";
  return;
}

template <typename Text>
struct KmerTrie
{
  struct Node
  {
    uint64_t edge_begin_offset;
    uint64_t edge_end_offset;
    std::map<uint64_t, std::shared_ptr<Node>> branches;
    uint64_t context_size;
    std::map<uint64_t, uint64_t> alphabet_count;

    Node () = default;
    Node
    (
      uint64_t const begin_offset,
      uint64_t const end_offset
    )
    : edge_begin_offset {begin_offset},
      edge_end_offset {end_offset},
      branches {},
      context_size {},
      alphabet_count {}
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
    ++(current_node->context_size);
    auto next_character {*k_mer_end};
    auto alphabet_count_iterator {current_node->alphabet_count.find(next_character)};
    if (alphabet_count_iterator == std::end(current_node->alphabet_count))
    {
      current_node->alphabet_count[next_character] = 1;
    }
    else
    {
      ++std::get<1>(*alphabet_count_iterator);
    }
    return;
  }

  double CalculateCumulative0thEmpiricalEntropy (std::shared_ptr<Node> current_node)
  {
    double cumulative_0th_empirical_entropy {};
    if (current_node->context_size == 0)
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
      current_node->context_size *
      CalculateZerothEmpiricalEntropy
      (
        current_node->alphabet_count,
        current_node->context_size
      );
    }
    return cumulative_0th_empirical_entropy;
  }
};

template
<
  typename File,
  typename Text
>
void PrintKthEmpiricalEntropy
(
  File &file,
  Text const &text,
  uint64_t const k
)
{
  KmerTrie k_mer_trie(text, k);
  file
  << std::fixed
  << std::setprecision(2)
  << "nH_" << k << ","
  << (k_mer_trie.kth_empirical_entropy * std::size(text) / (1024.0 * 1024.0 * 8.0)) << "(MB),"
  << "compression_ratio,"
  << (k_mer_trie.kth_empirical_entropy / text.width()) * 100.0
  << "%\n";
  return;
}

template <typename File>
void PrintBzip2CompressedSize
(
  std::filesystem::path const &text_path,
  File &file
)
{
  auto command {std::string{"bzip2 -k -9 -v "} + text_path.string()};
  if (std::system(command.c_str()) == 0)
  {
    auto bzip2_path {std::filesystem::path{text_path.string() + ".bz2"}};
    auto bzip2_file_size {static_cast<double>(std::filesystem::file_size(bzip2_path))};
    auto text_file_size {static_cast<double>(std::filesystem::file_size(text_path))};
    file
    << std::fixed
    << std::setprecision(2)
    << "bzip2,"
    << (bzip2_file_size / (1024.0 * 1024.0)) << "(MB),"
    << (bzip2_file_size / text_file_size) * 100.0 << "%\n";
    std::filesystem::remove(bzip2_path);
  }
  return;
}

template <typename File>
void PrintP7zipCompressedSize
(
  std::filesystem::path const &text_path,
  File &file
)
{
  auto command {std::string{"p7zip -k "} + text_path.string()};
  if (std::system(command.c_str()) == 0)
  {
    auto p7zip_path {std::filesystem::path{text_path.string() + ".7z"}};
    auto p7zip_file_size {static_cast<double>(std::filesystem::file_size(p7zip_path))};
    auto text_file_size {static_cast<double>(std::filesystem::file_size(text_path))};
    file
    << std::fixed
    << std::setprecision(2)
    << "p7zip,"
    << (p7zip_file_size / (1024.0 * 1024.0)) << "(MB),"
    << (p7zip_file_size / text_file_size) * 100.0 << "%\n";
    std::filesystem::remove(p7zip_path);
  }
  return;
}

template <typename File>
void PrintRepairCompressedSize
(
  std::filesystem::path const &text_path,
  File &file
)
{
  auto command {std::string{"repair "} + text_path.string()};
  if (std::system(command.c_str()) == 0)
  {
    auto repair_path {std::filesystem::path{text_path.string() + ".rp"}};
    auto repair_file_size {static_cast<double>(std::filesystem::file_size(repair_path))};
    auto text_file_size {static_cast<double>(std::filesystem::file_size(text_path))};
    file
    << std::fixed
    << std::setprecision(2)
    << "repair,"
    << (repair_file_size / (1024.0 * 1024.0)) << "(MB),"
    << (repair_file_size / text_file_size) * 100.0 << "%\n";
    std::filesystem::remove(repair_path);
  }
  return;
}

void PrintTextStatistics (std::filesystem::path const &text_path)
{
  sdsl::int_vector<> text;
  sdsl::load_from_file(text, text_path);
  if (*std::prev(std::end(text)) != 0)
  {
    sdsl::append_zero_symbol(text);
  }
  auto parent_statistics_path {CreateParentDirectoryByCategory("statistics", text_path)};
  auto statistics_path {CreatePath(parent_statistics_path, text_path.filename().string(), ".csv")};
  std::ofstream statistics_file {statistics_path};
  PrintTextSize(statistics_file, text);
  PrintAlphabetSize(statistics_file, text);
  Print0thEmpiricalEntropy(statistics_file, text);
  for (uint64_t power {0}; power != 4; ++power)
  {
    PrintKthEmpiricalEntropy(statistics_file, text, (1ULL << power));
  }
  PrintRunLengthCompressedSize(statistics_file, text);
  PrintBzip2CompressedSize(text_path, statistics_file);
  PrintP7zipCompressedSize(text_path, statistics_file);
  PrintRepairCompressedSize(text_path, statistics_file);
  return;
}
}
