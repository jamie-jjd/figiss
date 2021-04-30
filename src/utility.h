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

#include <sdsl/suffix_trees.hpp>
#include <sdsl/wavelet_trees.hpp>

namespace project
{
template
<
  typename File,
  typename Container
>
void Print
(
  File &file,
  Container const &container,
  int8_t step = 1,
  std::string const &delimiter = " ",
  std::string const &endmarker = "\n"
)
{
  if (std::size(container) != 0)
  {
    auto it {std::begin(container)};
    auto prev_end {std::prev(std::end(container))};
    if (step < 0)
    {
      it = std::prev(std::end(container));
      prev_end = std::begin(container);
    }
    while (it != prev_end)
    {
      file << *it << delimiter;
      std::advance(it, step);
    }
    file << *prev_end << endmarker;
  }
  return;
}

template
<
  typename File,
  typename Iterator
>
void Print
(
  File &file,
  Iterator first,
  Iterator last,
  int8_t const step = 1,
  std::string const &delimiter = " ",
  std::string const &endmarker = "\n"
)
{
  if (std::distance(first, last) * step > 0)
  {
    auto it {first};
    auto prev_last {std::prev(last, step)};
    while (it != prev_last)
    {
      file << *it << delimiter;
      std::advance(it, step);
    }
    file << *prev_last << endmarker;
  }
  return;
}

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

struct ConcatenatedPatterns
{
  uint64_t amount;
  uint64_t unit_size;
  sdsl::int_vector<8> labels;

  ConcatenatedPatterns () = default;

  ConcatenatedPatterns
  (
    uint64_t const amount_,
    uint64_t const unit_size_
  )
  : amount {amount_},
    unit_size {unit_size_}
  {
  }
};

void GenerateConcatnatedPatterns
(
  std::filesystem::path const &text_path,
  ConcatenatedPatterns &patterns
)
{
  sdsl::int_vector<8> text;
  sdsl::load_vector_from_file(text, text_path);
  if (patterns.unit_size <= std::size(text))
  {
    std::mt19937 engine {std::random_device{}()};
    std::uniform_int_distribution<uint64_t> distribution(0, std::size(text) - patterns.unit_size);
    auto random_begin_offset {std::bind(distribution, engine)};
    patterns.labels.resize(patterns.amount * patterns.unit_size);
    auto patterns_it {std::begin(patterns.labels)};
    for (uint64_t i {}; i != patterns.amount; ++i)
    {
      auto text_it {std::next(std::begin(text), random_begin_offset())};
      for (uint64_t j {}; j != patterns.unit_size; ++j)
      {
        *patterns_it = *text_it;
        ++patterns_it;
        ++text_it;
      }
    }
  }
  return;
}

auto GenerateAndSerializeConcatenatedPatterns (std::filesystem::path const &text_path)
{
  ConcatenatedPatterns patterns;
  GenerateConcatnatedPatterns(text_path, patterns);
  auto parent_patterns_path {CreateParentDirectoryByCategory("patterns", text_path)};
  auto patterns_path
  {
    CreatePath
    (
      parent_patterns_path,
      text_path.filename().string(),
      std::string{"."} + std::to_string(patterns.amount) + "." + std::to_string(patterns.unit_size)
    )
  };
  if (std::size(patterns.labels) != 0)
  {
    std::ofstream patterns_file {patterns_path};
    sdsl::write_member(patterns.amount, patterns_file);
    sdsl::write_member(patterns.unit_size, patterns_file);
    sdsl::serialize(patterns.labels, patterns_file);
  }
  return patterns_path;
}

template <typename Patterns>
void LoadConcatenatedPatterns
(
  std::filesystem::path const &patterns_path,
  Patterns &patterns
)
{
  std::ifstream patterns_file {patterns_file};
  sdsl::read_member(patterns.amount, patterns_file);
  sdsl::read_member(patterns.unit_size, patterns_file);
  patterns.labels.load(patterns_file);
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
    auto text_it {std::begin(text)};
    auto text_end {std::end(text)};
    while (text_it != text_end)
    {
      if (*text_it != 0)
      {
        codebook[*text_it] = 1;
      }
      ++text_it;
    }
  }
  std::partial_sum(std::begin(codebook), std::end(codebook), std::begin(codebook));
  {
    auto text_it {std::begin(text)};
    auto text_end {std::end(text)};
    while (text_it != text_end)
    {
      *text_it = codebook[*text_it];
      ++text_it;
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
  auto buffer_it {std::begin(buffer)};
  auto buffer_end {std::end(buffer)};
  while (buffer_it != buffer_end)
  {
    if (*buffer_it != 0)
    {
      *buffer_it = text[*buffer_it - 1];
    }
    ++buffer_it;
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
  auto it {std::begin(text)};
  while (it != std::end(text))
  {
    if (previous_character != *it)
    {
      ++runs;
      previous_character = *it;
    }
    ++it;
  }
  return runs;
}

template <typename Text>
double CalculateZerothEmpiricalEntropy (Text const &text)
{
  std::map<uint64_t, uint64_t> alphabet_count;
  auto text_it {std::begin(text)};
  while (text_it != std::end(text))
  {
    auto alphabet_count_it {alphabet_count.find(*text_it)};
    if (alphabet_count_it != std::end(alphabet_count))
    {
      ++std::get<1>(*alphabet_count_it);
    }
    else
    {
      alphabet_count[*text_it] = 1;
    }
    ++text_it;
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
    auto text_it {text_begin};
    auto text_last_it {std::prev(std::end(text), k)};
    while (text_it != text_last_it)
    {
      auto k_mer_begin {text_it};
      auto k_mer_end {std::next(k_mer_begin, k)};
      Insert(text_begin, k_mer_begin, k_mer_end);
      ++text_it;
    }
    kth_empirical_entropy = CalculateCumulative0thEmpiricalEntropy(root) / std::size(text);
  }

  template <typename TextIterator>
  void Insert
  (
    TextIterator text_begin,
    TextIterator k_mer_it,
    TextIterator k_mer_end
  )
  {
    auto current_node {root};
    while (k_mer_it != k_mer_end)
    {
      auto character {*k_mer_it};
      auto branches_it {current_node->branches.find(character)};
      if (branches_it == std::end(current_node->branches))
      {
        current_node = current_node->branches[character] = std::make_shared<Node>
        (
          std::distance(text_begin, k_mer_it),
          std::distance(text_begin, k_mer_end)
        );
        k_mer_it = k_mer_end;
      }
      else
      {
        auto child_node {std::get<1>(*branches_it)};
        auto edge_it {std::next(text_begin, child_node->edge_begin_offset)};
        auto edge_end {std::next(text_begin, child_node->edge_end_offset)};
        while
        (
          (k_mer_it != k_mer_end)
          &&
          (edge_it != edge_end)
          &&
          (*k_mer_it == *edge_it)
        )
        {
          ++k_mer_it;
          ++edge_it;
        }
        if
        (
          (k_mer_it != k_mer_end)
          &&
          (edge_it != edge_end)
        )
        {
          auto edge_it_offset {std::distance(text_begin, edge_it)};
          auto internal_node {std::make_shared<Node>(child_node->edge_begin_offset, edge_it_offset)};
          child_node->edge_begin_offset = edge_it_offset;
          internal_node->branches[*edge_it] = child_node;
          current_node = internal_node->branches[*k_mer_it] = std::make_shared<Node>
          (
            std::distance(text_begin, k_mer_it),
            std::distance(text_begin, k_mer_end)
          );
          std::get<1>(*branches_it) = internal_node;
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
      auto branches_it {std::begin(current_node->branches)};
      auto branches_end {std::end(current_node->branches)};
      while (branches_it != branches_end)
      {
        cumulative_0th_empirical_entropy += CalculateCumulative0thEmpiricalEntropy(std::get<1>(*branches_it));
        ++branches_it;
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

template
<
  typename WaveletTree,
  typename Text
>
auto CalculateRunLengthWaveletTreeSizeInMegaBytes (Text const &text)
{
  sdsl::wt_rlmn
  <
    sdsl::sd_vector<>,
    typename sdsl::sd_vector<>::rank_1_type,
    typename sdsl::sd_vector<>::select_1_type,
    WaveletTree
  >
  run_length_wavelet_tree;
  sdsl::construct_im(run_length_wavelet_tree, text);
  return sdsl::size_in_mega_bytes(run_length_wavelet_tree);
}

void PrintRunLengthWaveletTreeStatistics
(
  std::string const &category,
  std::filesystem::path const &bwt_path
)
{
  sdsl::int_vector<> bwt;
  sdsl::load_from_file(bwt, bwt_path);
  auto bwt_size_in_mega_bytes {std::filesystem::file_size(bwt_path) / (1024.0 * 1024.0)};
  auto rlwm_size_in_mega_bytes {CalculateRunLengthWaveletTreeSizeInMegaBytes<sdsl::wm_int<>>(bwt)};
  auto rlwt_hutu_size_in_mega_bytes {CalculateRunLengthWaveletTreeSizeInMegaBytes<sdsl::wt_hutu_int<>>(bwt)};
  auto parent_path {CreateParentDirectoryByCategory(category, bwt_path)};
  auto path {CreatePath(parent_path, bwt_path.filename().string())};
  std::ofstream file {path};
  file
  << std::fixed << std::setprecision(2)
  << bwt_path.filename().string() << "\n"
  << "\\multirow{2}{*}{" << bwt_size_in_mega_bytes << "} & "
  << rlwm_size_in_mega_bytes << " & "
  << rlwt_hutu_size_in_mega_bytes << " \\\\\n"
  << "(" << (rlwm_size_in_mega_bytes / bwt_size_in_mega_bytes) * 100.0 << "\\%) & "
  << "(" << (rlwt_hutu_size_in_mega_bytes / bwt_size_in_mega_bytes) * 100.0 << "\\%) \\\\\n";
  return;
}

template <typename Size>
std::string ProperSizeRepresentation (Size const size)
{
  std::stringstream stringstream;
  stringstream << std::fixed << std::setprecision(2);
  if (size > (1024 * 1024))
  {
    stringstream << (size / (1024.0 * 1024.0));
    return (stringstream.str() + "M");
  }
  else if (size > 1024)
  {
    stringstream << (size / 1024.0);
    return (stringstream.str() + "K");
  }
  else
  {
    stringstream << (size * 1.0);
    return (stringstream.str());
  }
  return "0";
}
}
