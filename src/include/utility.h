#pragma once

#include <cmath>

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <map>
#include <memory>

#include <sdsl/csa_wt.hpp>

namespace figiss
{
template <typename File, typename Container>
void Print
(
  Container const& container,
  File& file,
  int8_t const step = 1,
  std::string const& separator = " ",
  std::string const& endmarker = "\n"
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
  File& file,
  int8_t const step = 1,
  std::string const& separator = " ",
  std::string const& endmarker = "\n"
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

void GeneratePrefix
(
  std::filesystem::path const& text_path,
  uint64_t const size_in_megabytes
)
{
  sdsl::int_vector<8> text;
  sdsl::load_vector_from_file(text, text_path);
  auto path
  {
    std::filesystem::path
    {
      text_path.filename().string()
      + "." + std::to_string(size_in_megabytes) + "mb"
    }
  };
  std::fstream fout {path, std::ios_base::out | std::ios_base::trunc};
  for (uint64_t i {}; i != (size_in_megabytes * (1ULL << 20)); ++i)
  {
    fout << text[i];
  }
  return;
}

template <typename Text>
double CalculateZerothEmpiricalEntropy (Text const& text)
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
  for (auto const& character_count : alphabet_count)
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
    Text const& text,
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
    return (stringstream.str() + "Gi");
  }
  else if (size >= (1ULL << 20))
  {
    stringstream << (static_cast<double>(size) / (1ULL << 20));
    return (stringstream.str() + "Mi");
  }
  else if (size >= (1ULL << 10))
  {
    stringstream << (static_cast<double>(size) / (1ULL << 10));
    return (stringstream.str() + "Ki");
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
  SpaceNode (std::string const& name, uint64_t const size_in_bytes = 0);

  void AccumalateSizeInBytes (uint64_t const size_in_bytes);
  void AddChild (std::shared_ptr<SpaceNode> child);
  void AddLeaf (std::string const& name, uint64_t const size_in_bytes);

  inline uint64_t GetSizeInBytes () const
  {
    return size_in_bytes_;
  }

  friend std::ostream& operator<<
  (
    std::ostream& out,
    std::pair<std::shared_ptr<SpaceNode>, bool> pair
  );

private:

  std::string name_;
  uint64_t size_in_bytes_;
  std::deque<std::shared_ptr<SpaceNode>> children_;

};

SpaceNode::SpaceNode (std::string const& name, uint64_t const size_in_bytes)
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

void SpaceNode::AddLeaf (std::string const& name, uint64_t const size_in_bytes)
{
  auto node {std::make_shared<SpaceNode>(name, size_in_bytes)};
  children_.emplace_back(node);
  size_in_bytes_ += size_in_bytes;
  return;
}

std::ostream& operator<<
(
  std::ostream& out,
  std::pair<std::shared_ptr<SpaceNode>, bool> pair
)
{
  std::deque<std::pair<std::shared_ptr<SpaceNode>, uint64_t>> nodes;
  auto root {std::get<0>(pair)};
  auto is_proper {std::get<1>(pair)};
  nodes.emplace_back(root, 0);
  while (!nodes.empty())
  {
    auto const node {std::get<0>(nodes.back())};
    auto const depth {std::get<1>(nodes.back())};
    nodes.pop_back();
    std::string whitespaces(depth * 2, ' ');
    out << whitespaces << node->name_ << ":";
    if (is_proper)
    {
      out << ProperSizeRepresentation(node->size_in_bytes_);
    }
    else
    {
      out << node->size_in_bytes_;
    }
    out << "B\n";
    for (auto it {std::rbegin(node->children_)}; it != std::rend(node->children_); ++it)
    {
      nodes.emplace_back(*it, depth + 1);
    }
  }
  return out;
}

template <typename Index>
void PrintIndexSpace
(
  std::filesystem::path const& byte_text_path,
  Index& index,
  bool const is_proper_representation = false
)
{
  std::filesystem::path output_path {std::string{"../data/space/"} + byte_text_path.filename().string()};
  if (!std::filesystem::exists(output_path.parent_path()))
  {
    std::filesystem::create_directories(output_path.parent_path());
  }
  std::fstream fout {output_path, std::ios_base::out | std::ios_base::trunc};
  std::cout << "write space information to " << std::filesystem::canonical(output_path) << "\n";
  {
    auto index_path {std::filesystem::path("_.index")};
    auto root {std::make_shared<SpaceNode>("index")};
    index.Serialize(index_path, root);
    fout << std::make_pair(root, is_proper_representation);
    std::cout << "remove " << std::filesystem::canonical(index_path) << "\n";
    std::filesystem::remove(index_path);
  }
  return;
}

void PrintRlfmSpace (std::filesystem::path const& text_path, bool const is_proper_representation = false)
{
  auto output_path {std::filesystem::path{std::string{"../data/space/rlfm/"} + text_path.filename().string()}};
  if (!std::filesystem::exists(output_path.parent_path()))
  {
    std::filesystem::create_directories(output_path.parent_path());
  }
  std::fstream fout {output_path, std::ios_base::out | std::ios_base::trunc};
  std::cout << "write space information to " << std::filesystem::canonical(output_path) << "\n";
  {
    sdsl::csa_wt<sdsl::wt_rlmn<>, 0xFFFF'FFFF, 0xFFFF'FFFF> rlfm;
    {
      sdsl::int_vector<8> text;
      sdsl::load_vector_from_file(text, text_path);
      std::cout << "construct rlfm of " << std::filesystem::canonical(text_path) << "\n";
      sdsl::construct_im(rlfm, text);
      if (is_proper_representation)
      {
        fout << ProperSizeRepresentation(sdsl::size_in_bytes(rlfm)) << "\n";
      }
      else
      {
        fout << sdsl::size_in_bytes(rlfm) << "\n";
      }
    }
  }
  return;
}
}
