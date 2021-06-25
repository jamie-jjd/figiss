#pragma once

#include <cmath>

#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <map>
#include <memory>
#include <random>

#include <sdsl/csa_wt.hpp>
#include <sdsl/int_vector.hpp>

namespace project
{
template <typename File, typename Container>
void Print
(
  Container const &container,
  File &file,
  int8_t const step = 1,
  std::string const &separator = " ",
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
      file << static_cast<uint64_t>(*it) << separator;
      std::advance(it, step);
    }
    file << static_cast<uint64_t>(*prev_end) << endmarker;
  }
  return;
}

template <typename File, typename Iterator>
void Print
(
  Iterator first,
  Iterator last,
  File &file,
  int8_t const step = 1,
  std::string const &separator = " ",
  std::string const &endmarker = "\n"
)
{
  if (std::distance(first, last) * step > 0)
  {
    auto it {first};
    auto prev_last {std::prev(last, step)};
    while (it != prev_last)
    {
      file << static_cast<uint64_t>(*it) << separator;
      std::advance(it, step);
    }
    file << static_cast<uint64_t>(*prev_last) << endmarker;
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
  std::fstream prefix_file(prefix_path, std::ios_base::out | std::ios_base::trunc);
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

  ConcatenatedPatterns
  (
    uint64_t const amount_ = 0,
    uint64_t const unit_size_ = 0
  )
  : amount {amount_},
    unit_size {unit_size_}
  {
  }
};

void GenerateConcatenatedPatterns
(
  std::filesystem::path const &text_path,
  ConcatenatedPatterns &patterns
)
{
  sdsl::int_vector<8> text;
  sdsl::load_vector_from_file(text, text_path);
  if (patterns.unit_size > std::size(text))
  {
    std::cout << "warning: pattern size is replaced with text size\n";
    patterns.unit_size = std::size(text);
  }
  std::cout
  << "generate " << patterns.amount << " random patterns of size " << patterns.unit_size
  << " from " << std::filesystem::canonical(text_path) << "\n" ;
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
  return;
}

auto GenerateAndSerializeConcatenatedPatterns
(
  std::filesystem::path const &text_path,
  ConcatenatedPatterns & patterns
)
{
  GenerateConcatenatedPatterns(text_path, patterns);
  auto parent_patterns_path {CreateParentDirectoryByCategory("patterns", text_path)};
  auto patterns_path
  {
    CreatePath
    (
      parent_patterns_path,
      text_path.filename().string(),
      "." + std::to_string(patterns.amount) + "." + std::to_string(patterns.unit_size)
    )
  };
  if (std::size(patterns.labels) != 0)
  {
    std::fstream patterns_file {patterns_path, std::ios_base::out | std::ios_base::trunc};
    std::cout << "serialize patterns into " << std::filesystem::canonical(patterns_path) << "\n";
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
  std::ifstream patterns_file {patterns_path};
  std::cout << "load patterns from " << std::filesystem::canonical(patterns_path) << "\n";
  sdsl::read_member(patterns.amount, patterns_file);
  sdsl::read_member(patterns.unit_size, patterns_file);
  patterns.labels.load(patterns_file);
  return;
}

template <typename Text>
uint64_t CalculateRunsSize (Text const &text)
{
  uint64_t runs_size {};
  uint64_t prev_character {std::numeric_limits<uint64_t>::max()};
  auto it {std::begin(text)};
  while (it != std::end(text))
  {
    if (prev_character != *it)
    {
      ++runs_size;
      prev_character = *it;
    }
    ++it;
  }
  return runs_size;
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

template <typename Size>
std::string ProperSizeRepresentation (Size const size)
{
  std::stringstream stringstream;
  stringstream << std::fixed << std::setprecision(2);
  if (size >= (1ULL << 30))
  {
    stringstream << (static_cast<double>(size) / (1ULL << 30));
    return (stringstream.str() + "G");
  }
  else if (size >= (1ULL << 20))
  {
    stringstream << (static_cast<double>(size) / (1ULL << 20));
    return (stringstream.str() + "M");
  }
  else if (size >= (1ULL << 10))
  {
    stringstream << (static_cast<double>(size) / (1ULL << 10));
    return (stringstream.str() + "K");
  }
  else
  {
    stringstream << static_cast<double>(size);
    return stringstream.str();
  }
  return "0.00";
}

std::string ProperTimeRepresentation (double const nanoseconds)
{
  std::stringstream stringstream;
  stringstream << std::fixed << std::setprecision(2);
  if (nanoseconds >= 1'000'000'000.0)
  {
    stringstream << (nanoseconds / 1'000'000'000.0);
    return (stringstream.str() + "s");
  }
  else if (nanoseconds >= 1'000'000.0)
  {
    stringstream << (nanoseconds / 1'000'000.0);
    return (stringstream.str() + "ms");
  }
  else if (nanoseconds >= 1'000.0)
  {
    stringstream << (nanoseconds / 1'000.0);
    return (stringstream.str() + "us");
  }
  else
  {
    stringstream << nanoseconds;
    return (stringstream.str() + "ns");
  }
  return "0.00";
}

class SpaceNode
{
public:

  SpaceNode () = default;
  SpaceNode (std::string const &name, uint64_t const size_in_bytes = 0) noexcept;

  void AccumalateSizeInBytes (uint64_t const size_in_bytes);
  void AddChild (std::shared_ptr<SpaceNode> child);
  void AddLeaf (std::string const &name, uint64_t const size_in_bytes);

  inline uint64_t GetSizeInBytes () const noexcept
  {
    return size_in_bytes_;
  }

  friend std::ostream& operator<< (std::ostream &out, std::shared_ptr<SpaceNode> root);

private:

  std::string name_;
  uint64_t size_in_bytes_;
  std::deque<std::shared_ptr<SpaceNode>> children_;

};

SpaceNode::SpaceNode (std::string const &name, uint64_t const size_in_bytes) noexcept
: name_ {name},
  size_in_bytes_ {size_in_bytes},
  children_ {}
{
}

void SpaceNode::AccumalateSizeInBytes (uint64_t const size_in_bytes)
{
  size_in_bytes_ += size_in_bytes;
  return;
}

void SpaceNode::AddChild (std::shared_ptr<SpaceNode> child)
{
  children_.emplace_back(child);
  size_in_bytes_ += child->GetSizeInBytes();
  return;
}

void SpaceNode::AddLeaf (std::string const &name, uint64_t const size_in_bytes)
{
  auto node {std::make_shared<SpaceNode>(name, size_in_bytes)};
  children_.emplace_back(node);
  size_in_bytes_ += size_in_bytes;
  return;
}

std::ostream& operator<< (std::ostream &out, std::shared_ptr<SpaceNode> root)
{
  std::deque<std::pair<std::shared_ptr<SpaceNode>, uint64_t>> nodes;
  nodes.emplace_back(root, 0);
  while (!nodes.empty())
  {
    auto const node {std::get<0>(nodes.back())};
    auto const depth {std::get<1>(nodes.back())};
    nodes.pop_back();
    std::string whitespaces(depth * 2, ' ');
    out << whitespaces << node->name_ << ": " << ProperSizeRepresentation(node->size_in_bytes_) << "\n";
    for (auto it {std::rbegin(node->children_)}; it != std::rend(node->children_); ++it)
    {
      nodes.emplace_back(*it, depth + 1);
    }
  }
  return out;
}

template <typename Index>
void PrintSpace (Index &index, std::filesystem::path const &text_path)
{
  auto parent_space_path {CreateParentDirectoryByCategory("space", text_path)};
  auto space_path {CreatePath(parent_space_path, text_path.filename().string())};
  std::fstream space_file {space_path, std::ios_base::out | std::ios_base::trunc};
  std::cout << "write space information to " << std::filesystem::canonical(space_path) << "\n";
  {
    auto parent_index_path {CreateParentDirectoryByCategory("index", text_path)};
    auto index_path {CreatePath(parent_index_path, text_path.filename().string(), ".index")};
    Index index {text_path};
    auto root {std::make_shared<SpaceNode>("index")};
    index.Serialize(index_path, root);
    space_file << root;
  }
  {
    auto parent_index_path {CreateParentDirectoryByCategory("index", text_path)};
    auto index_path {CreatePath(parent_index_path, text_path.filename().string(), ".rlfm")};
    sdsl::int_vector<8> text;
    sdsl::load_vector_from_file(text, text_path);
    sdsl::csa_wt
    <
      sdsl::wt_rlmn<>,
      std::numeric_limits<int32_t>::max(),
      std::numeric_limits<int32_t>::max()
    >
    rlfm;
    sdsl::construct_im(rlfm, text);
    std::fstream index_file {index_path, std::ios_base::out | std::ios_base::trunc};
    space_file << ProperSizeRepresentation(sdsl::serialize(rlfm, index_file)) << "\n";
  }
  space_file << ProperSizeRepresentation(std::filesystem::file_size(text_path)) << "\n";
  return;
}

template <typename Times>
void PrintTime
(
  std::string const &category,
  std::filesystem::path const &text_path,
  Times const &times
)
{
  auto parent_path {CreateParentDirectoryByCategory(category, text_path)};
  auto path {CreatePath(parent_path, text_path.filename().string())};
  std::fstream file {path, std::ios_base::out | std::ios_base::trunc};
  file << std::fixed << std::setprecision(2);
  for (auto pair : times)
  {
    file
    << std::get<0>(pair) << ","
    << ProperTimeRepresentation(std::get<1>(pair)) << "/char"
    << "\n";
  }
  return;
}
}
