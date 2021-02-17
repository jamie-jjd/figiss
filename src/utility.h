#pragma once

#include <cmath>

#include <filesystem>
#include <fstream>
#include <functional>
#include <map>
#include <memory>
#include <random>
#include <string>
#include <unordered_map>
#include <unordered_set>

#include <sdsl/construct.hpp>
#include <sdsl/int_vector.hpp>

namespace project
{
auto GeneratePattern
(
  std::filesystem::path const text_path,
  uint64_t const pattern_amount,
  uint64_t const pattern_size
)
-> std::filesystem::path
{
  std::ifstream text_file {text_path};
  sdsl::int_vector<> text;
  text.load(text_file);
  auto parent_pattern_path
  {
    std::filesystem::path{"../data/pattern"}
    / text_path.parent_path().filename()
  };
  if (!std::filesystem::exists(parent_pattern_path))
  {
    std::filesystem::create_directories(parent_pattern_path);
  }
  if (pattern_size < std::size(text))
  {
    auto pattern_path
    {
      parent_pattern_path
      /
      std::filesystem::path
      {
        text_path.filename().string()
        + "." + std::to_string(pattern_amount)
        + "." + std::to_string(pattern_size)
      }
    };
    std::ofstream pattern_file {pattern_path};
    std::mt19937 engine {std::random_device{}()};
    std::uniform_int_distribution<uint64_t> distribution(0, std::size(text) - pattern_size);
    auto random_begin_offset {std::bind(distribution, engine)};
    sdsl::int_vector<> pattern_amount_size(2, 0, 64);
    pattern_amount_size[0] = pattern_amount;
    pattern_amount_size[1] = pattern_size;
    pattern_amount_size.serialize(pattern_file);
    sdsl::int_vector<> patterns(pattern_amount * pattern_size, 0, text.width());
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
    patterns.serialize(pattern_file);
    return pattern_path;
  }
  else
  {
    auto pattern_path
    {
      std::filesystem::path
      {
        text_path.filename().string()
        + ".1"
        + "." + std::to_string(std::size(text))
      }
    };
    std::ofstream pattern_file {pattern_path};
    std::cout
    << "pattern size >= text size\n"
    << "generated pattern is the text\n";
    pattern_file << "1\n";
    text.serialize(pattern_file);
    return pattern_path;
  }
  return std::filesystem::path{};
}

void ByteToCompactText (std::filesystem::path const text_path)
{
  sdsl::int_vector<> text;
  sdsl::load_vector_from_file(text, text_path);
  sdsl::int_vector<> codebook(256, 0);
  {
    auto text_iterator {std::begin(text)};
    auto text_end {std::end(text)};
    while (text_iterator != text_end)
    {
      codebook[*text_iterator] = 1;
      ++text_iterator;
    }
  }
  std::partial_sum
  (
    std::begin(codebook),
    std::end(codebook),
    std::begin(codebook)
  );
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
  std::filesystem::path parent_compact_text_path
  {
    std::filesystem::path("../data/compact")
    / text_path.parent_path().filename()
  };
  if (!std::filesystem::exists(parent_compact_text_path))
  {
    std::filesystem::create_directories(parent_compact_text_path);
  }
  std::filesystem::path compact_text_path
  {
    parent_compact_text_path
    / (text_path.filename().string() + ".sdsl")
  };
  std::ofstream compact_text_file {compact_text_path};
  text.serialize(compact_text_file);
  return;
}

template <typename Text>
void PrintRunAmount
(
  std::ofstream &statistics_file,
  Text const &text
)
{
  uint64_t previous_character {0};
  uint64_t run_amount {0};
  for (auto const &character : text)
  {
    if (character != previous_character)
    {
      ++run_amount;
      previous_character = character;
    }
  }
  statistics_file << "run_amount," << run_amount << "\n";
  return;
}

template <typename Text>
void PrintAlphabetSize
(
  std::ofstream &statistics_file,
  Text const &text
)
{
  std::unordered_set<uint64_t> alphabet;
  for (auto const &character : text)
  {
    alphabet.insert(character);
  }
  statistics_file << "alphabet_size," << std::size(alphabet) << "\n";
  return;
}

template <typename Text>
void Print0thEmpiricalEntropy
(
  std::ofstream &statistics_file,
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
  double zeroth_empirical_entropy {};
  auto log_text_size {std::log2(std::size(text))};
  for (auto const &character_count : alphabet_count)
  {
    auto count {static_cast<double>(std::get<1>(character_count))};
    zeroth_empirical_entropy += (count * (log_text_size - std::log2(count)));
  }
  zeroth_empirical_entropy /= std::size(text);
  statistics_file
  << "0th_empirical_entropy,"
  << zeroth_empirical_entropy
  << ",compression_ratio,"
  << ((std::size(text) * zeroth_empirical_entropy) / text.bit_size()) * 100.0
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
      auto log_context_size {std::log2(static_cast<double>(current_node->context_size))};
      for (auto const &character_count : current_node->alphabet_count)
      {
        auto count {static_cast<double>(std::get<1>(character_count))};
        cumulative_0th_empirical_entropy += count * (log_context_size - std::log2(count));
      }
    }
    return cumulative_0th_empirical_entropy;
  }
};

template <typename Text>
void PrintKthEmpiricalEntropy
(
  std::ofstream &statistics_file,
  Text const &text,
  uint64_t const k
)
{
  KmerTrie k_mer_trie(text, k);
  statistics_file
  << k << "th_empirical_entropy,"
  << k_mer_trie.kth_empirical_entropy
  << ",compression_ratio,"
  << ((std::size(text) * k_mer_trie.kth_empirical_entropy) / text.bit_size()) * 100.0
  << "%\n";
  return;
}

void PrintP7zipCompressedSize
(
  std::filesystem::path const &text_path,
  std::ofstream &statistics_file
)
{
  auto command {std::string{"p7zip -k "} + text_path.string()};
  if (std::system(command.c_str()) == 0)
  {
    auto p7zip_path {std::filesystem::path{text_path.string() + ".7z"}};
    auto p7zip_file_size {static_cast<double>(std::filesystem::file_size(p7zip_path))};
    auto text_file_size {static_cast<double>(std::filesystem::file_size(text_path))};
    std::filesystem::remove(p7zip_path);
    statistics_file << "p7zip," << (p7zip_file_size / text_file_size) * 100.0 << "%\n";
  }
  return;
}

void PrintBzip2CompressedSize
(
  std::filesystem::path const &text_path,
  std::ofstream &statistics_file
)
{
  auto command {std::string{"bzip2 -k -9 -v "} + text_path.string()};
  if (std::system(command.c_str()) == 0)
  {
    auto bzip2_path {std::filesystem::path{text_path.string() + ".bz2"}};
    auto bzip2_file_size {static_cast<double>(std::filesystem::file_size(bzip2_path))};
    auto text_file_size {static_cast<double>(std::filesystem::file_size(text_path))};
    std::filesystem::remove(bzip2_path);
    statistics_file << "bzip2," << (bzip2_file_size / text_file_size) * 100.0 << "%\n";
  }
  return;
}

void PrintRepairCompressedSize
(
  std::filesystem::path const &text_path,
  std::ofstream &statistics_file
)
{
  auto command {std::string{"repair "} + text_path.string()};
  if (std::system(command.c_str()) == 0)
  {
    auto repair_path {std::filesystem::path{text_path.string() + ".rp"}};
    auto repair_file_size {static_cast<double>(std::filesystem::file_size(repair_path))};
    auto text_file_size {static_cast<double>(std::filesystem::file_size(text_path))};
    std::filesystem::remove(repair_path);
    statistics_file << "repair," << (repair_file_size / text_file_size) * 100.0 << "%\n";
  }
  return;
}

void PrintTextStatistics (std::filesystem::path const &text_path)
{
  std::ifstream text_file {text_path};
  sdsl::int_vector<> text;
  text.load(text_file);
  auto parent_statistics_path
  {
    std::filesystem::path{"../data/statistics"}
    / text_path.parent_path().filename()
  };
  if (!std::filesystem::exists(parent_statistics_path))
  {
    std::filesystem::create_directories(parent_statistics_path);
  }
  auto statistics_path {parent_statistics_path / (text_path.filename().string() + ".csv")};
  std::ofstream statistics_file {statistics_path};
  statistics_file << "text_size," << std::size(text) << "\n";
  statistics_file << "text_bit_width," << static_cast<uint64_t>(text.width()) << "\n";
  PrintAlphabetSize(statistics_file, text);
  PrintRunAmount(statistics_file, text);
  Print0thEmpiricalEntropy(statistics_file, text);
  for (uint64_t power {0}; power != 4; ++power)
  {
    PrintKthEmpiricalEntropy(statistics_file, text, (1ULL << power));
  }
  PrintP7zipCompressedSize(text_path, statistics_file);
  PrintBzip2CompressedSize(text_path, statistics_file);
  PrintRepairCompressedSize(text_path, statistics_file);
  return;
}
}
