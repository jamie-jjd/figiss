#pragma once

#include "utility.h"

namespace figiss
{
class Trie
{
public:

  struct Node
  {
    std::map<uint8_t, std::shared_ptr<Node>> children;
    uint64_t count;
  };

  Trie (): root_ {std::make_shared<Node>()} {}

  template <typename Iterator>
  void Insert (Iterator it, Iterator last);

  inline auto GetRoot () const
  {
    return root_;
  }

  friend std::ostream& operator<< (std::ostream &out, std::pair<Trie, bool> const &pair);

private:

  std::shared_ptr<Node> root_;

};

template <typename Iterator>
void Trie::Insert (Iterator it, Iterator end)
{
  auto node {root_};
  while (it != end)
  {
    if (node->children.find(*it) == node->children.end())
    {
      node->children[*it] = std::make_shared<Node>();
    }
    node = node->children[*it];
    ++(node->count);
    ++it;
  }
  return;
}

std::ostream& operator<< (std::ostream &out, std::pair<Trie, bool> const &pair)
{
  out << "trie:\n";
  auto const &trie {std::get<0>(pair)};
  auto const is_level_order {std::get<1>(pair)};
  std::deque<std::tuple<std::shared_ptr<Trie::Node>, uint64_t, uint64_t>> nodes;
  nodes.emplace_back(trie.root_, 0, 0);
  uint64_t size {};
  while (!nodes.empty())
  {
    ++size;
    auto node {std::get<0>(nodes.front())};
    auto label {std::get<1>(nodes.front())};
    auto depth {std::get<2>(nodes.front())};
    if (is_level_order)
    {
      nodes.pop_front();
    }
    else
    {
      node = std::get<0>(nodes.back());
      label = std::get<1>(nodes.back());
      depth = std::get<2>(nodes.back());
      nodes.pop_back();
    }
    if (depth != 0)
    {
      out << depth << ":" << label << "(" << node->count << ")\n";
    }
    if (is_level_order)
    {
      for (auto it {std::begin(node->children)}; it != std::end(node->children); ++it)
      {
        nodes.emplace_back(std::get<1>(*it), std::get<0>(*it), depth + 1);
      }
    }
    else
    {
      for (auto it {std::rbegin(node->children)}; it != std::rend(node->children); ++it)
      {
        nodes.emplace_back(std::get<1>(*it), std::get<0>(*it), depth + 1);
      }
    }
  }
  out << "|nodes|:\n";
  out << (size - 1) << "\n";
  return out;
}

class SubFactorTrie
{
public:

  SubFactorTrie () = default;
  SubFactorTrie (SubFactorTrie const&);
  SubFactorTrie (SubFactorTrie&&);
  SubFactorTrie (Trie const &trie);
  SubFactorTrie& operator= (SubFactorTrie const&) = delete;
  SubFactorTrie& operator= (SubFactorTrie&&);

  void Swap (SubFactorTrie&);

  uint64_t Serialize
  (
    std::ostream &out,
    std::shared_ptr<SpaceNode> parent = nullptr,
    std::string const name = ""
  );
  void Load (std::istream &in);

  template <typename Iterator>
  uint64_t Count (Iterator it, Iterator end);

  friend std::ostream& operator<< (std::ostream &out, std::pair<SubFactorTrie, bool> const &pair);

private:

  sdsl::bit_vector level_order_bits_;
  sdsl::bit_vector::select_1_type level_order_select_1_;
  sdsl::int_vector<> labels_;
  sdsl::int_vector<> counts_;

};

SubFactorTrie::SubFactorTrie (SubFactorTrie const& sub_factor_trie)
{
  if (this != &sub_factor_trie)
  {
    auto temp {SubFactorTrie()};
    this->Swap(temp);
    {
      level_order_bits_ = decltype(level_order_bits_)(sub_factor_trie.level_order_bits_);
      level_order_select_1_ = decltype(level_order_select_1_)(&sub_factor_trie.level_order_bits_);
      labels_ = decltype(labels_)(sub_factor_trie.labels_);
      counts_ = decltype(counts_)(sub_factor_trie.counts_);
    }
  }
}

SubFactorTrie::SubFactorTrie (SubFactorTrie&& sub_factor_trie)
{
  if (this != &sub_factor_trie)
  {
    this->Swap(sub_factor_trie);
  }
}

SubFactorTrie::SubFactorTrie (Trie const &trie)
{
  std::vector<bool> level_order_bits;
  std::vector<uint8_t> labels;
  std::vector<uint64_t> counts;
  std::deque<std::shared_ptr<Trie::Node>> nodes;
  nodes.emplace_back(trie.GetRoot());
  while (!nodes.empty())
  {
    auto node {nodes.front()};
    nodes.pop_front();
    for (auto const &pair : node->children)
    {
      auto label {std::get<0>(pair)};
      auto child {std::get<1>(pair)};
      nodes.emplace_back(child);
      level_order_bits.emplace_back(0);
      labels.emplace_back(label);
      counts.emplace_back(child->count);
    }
    level_order_bits.emplace_back(1);
  }
  {
    level_order_bits_.resize(std::size(level_order_bits));
    std::copy(std::begin(level_order_bits), std::end(level_order_bits), std::begin(level_order_bits_));
    level_order_select_1_ = decltype(level_order_select_1_)(&level_order_bits_);
  }
  {
    labels_.resize(std::size(labels));
    std::copy(std::begin(labels), std::end(labels), std::begin(labels_));
    sdsl::util::bit_compress(labels_);
  }
  {
    counts_.resize(std::size(counts));
    std::copy(std::begin(counts), std::end(counts), std::begin(counts_));
    sdsl::util::bit_compress(counts_);
  }
  return;
}

SubFactorTrie& SubFactorTrie::operator= (SubFactorTrie &&sub_factor_trie)
{
  if (this != &sub_factor_trie)
  {
    SubFactorTrie temp {std::move(sub_factor_trie)};
    this->Swap(temp);
  }
  return *this;
}

void SubFactorTrie::Swap (SubFactorTrie& sub_factor_trie)
{
  if (this != &sub_factor_trie)
  {
    level_order_bits_.swap(sub_factor_trie.level_order_bits_);
    level_order_select_1_.swap(sub_factor_trie.level_order_select_1_);
    level_order_select_1_.set_vector(&level_order_bits_);
    sub_factor_trie.level_order_select_1_.set_vector(&sub_factor_trie.level_order_bits_);
    labels_.swap(sub_factor_trie.labels_);
    counts_.swap(sub_factor_trie.counts_);
  }
  return;
}

uint64_t SubFactorTrie::Serialize
(
  std::ostream &out,
  std::shared_ptr<SpaceNode> parent,
  std::string const name
)
{
  uint64_t size_in_bytes {};
  if (!parent)
  {
    sdsl::serialize(level_order_bits_, out);
    sdsl::serialize(level_order_select_1_, out);
    sdsl::serialize(labels_, out);
    sdsl::serialize(counts_, out);
  }
  else
  {
    auto node {std::make_shared<SpaceNode>(name)};
    node->AddLeaf("level_order_bits_", sdsl::serialize(level_order_bits_, out));
    node->AddLeaf("level_order_select_1_", sdsl::serialize(level_order_select_1_, out));
    node->AddLeaf("labels_", sdsl::serialize(labels_, out));
    node->AddLeaf("counts_", sdsl::serialize(counts_, out));
    parent->AddChild(node);
    size_in_bytes = node->GetSizeInBytes();
  }
  return size_in_bytes;
}

void SubFactorTrie::Load (std::istream &in)
{
  level_order_bits_.load(in);
  level_order_select_1_.load(in);
  level_order_select_1_.set_vector(&level_order_bits_);
  labels_.load(in);
  counts_.load(in);
  return;
}

template <typename Iterator>
uint64_t SubFactorTrie::Count (Iterator it, Iterator end)
{
  uint64_t offset {};
  uint64_t end_offset {level_order_select_1_(1)};
  while (true)
  {
    while ((offset != end_offset) && (*it != labels_[offset]))
    {
      ++offset;
    }
    if (offset != end_offset)
    {
      if (++it != end)
      {
        end_offset = level_order_select_1_(offset + 2) - (offset + 2) + 1;
        offset = level_order_select_1_(offset + 1) - (offset + 1) + 1;
      }
      else
      {
        return counts_[offset];
      }
    }
    else
    {
      return 0;
    }
  }
  return 0;
}

std::ostream& operator<< (std::ostream &out, std::pair<SubFactorTrie, bool> const &pair)
{
  out << "trie:\n";
  auto const &trie {std::get<0>(pair)};
  auto const is_level_order {std::get<1>(pair)};
  std::deque<std::pair<uint64_t, uint64_t>> offsets;
  offsets.emplace_back(0, 0);
  while (!offsets.empty())
  {
    auto offset {std::get<0>(offsets.front())};
    auto depth {std::get<1>(offsets.front())};
    if (is_level_order)
    {
      offsets.pop_front();
    }
    else
    {
      offset = std::get<0>(offsets.back());
      depth = std::get<1>(offsets.back());
      offsets.pop_back();
    }
    uint64_t end_offset {};
    if (depth != 0)
    {
      out << depth << ":" << trie.labels_[offset] << "(" << trie.counts_[offset] << ")\n";
      end_offset = trie.level_order_select_1_(offset + 2) - (offset + 2) + 1;
      offset = trie.level_order_select_1_(offset + 1) - (offset + 1) + 1;
    }
    else
    {
      end_offset = trie.level_order_select_1_(1);
    }
    if (is_level_order)
    {
      while (offset != end_offset)
      {
        offsets.emplace_back(offset, depth + 1);
        ++offset;
      }
    }
    else
    {
      std::swap(--offset, --end_offset);
      while (offset != end_offset)
      {
        offsets.emplace_back(offset, depth + 1);
        --offset;
      }
    }
  }
  out << "|nodes|:\n";
  out << std::size(trie.labels_) << "\n";
  out << "space:\n";
  out << "level_order_bits_: " << ProperSizeRepresentation(sdsl::size_in_bytes(trie.level_order_bits_)) << "B\n";
  out << "level_order_select_1_: " << ProperSizeRepresentation(sdsl::size_in_bytes(trie.level_order_select_1_)) << "B\n";
  out << "labels_: " << ProperSizeRepresentation(sdsl::size_in_bytes(trie.labels_)) << "B\n";
  out << "counts_: " << ProperSizeRepresentation(sdsl::size_in_bytes(trie.counts_)) << "B\n";
  return out;
}
}
