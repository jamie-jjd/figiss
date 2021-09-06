#pragma once

#include "utility.h"

namespace figiss
{
class SparsePrefixSum
{
public:

  SparsePrefixSum () = default;
  template <typename RandomAccessContainer>
  SparsePrefixSum (RandomAccessContainer const& prefix_sum);
  SparsePrefixSum& operator= (SparsePrefixSum&&);

  uint64_t Serialize
  (
    std::ostream& out,
    std::shared_ptr<SpaceNode> parent = nullptr,
    std::string const name = ""
  );
  void Load (std::istream& in);

  uint64_t operator[] (uint64_t const index) const;

  friend std::ostream& operator<< (std::ostream& out, SparsePrefixSum const& sparse_prefix_sum);

private:

  sdsl::sd_vector<> prefix_sum_bits_;
  sdsl::sd_vector<>::select_1_type prefix_sum_select_1_;

};

template <typename RandomAccessContainer>
SparsePrefixSum::SparsePrefixSum (RandomAccessContainer const& prefix_sum)
{
  prefix_sum_bits_ = decltype(prefix_sum_bits_)(std::begin(prefix_sum), std::end(prefix_sum));
  prefix_sum_select_1_.set_vector(&prefix_sum_bits_);
}

SparsePrefixSum& SparsePrefixSum::operator= (SparsePrefixSum&& sparse_prefix_sum)
{
  if (this != &sparse_prefix_sum)
  {
    prefix_sum_bits_ = std::move(sparse_prefix_sum.prefix_sum_bits_);
    prefix_sum_select_1_.set_vector(&prefix_sum_bits_);
  }
  return *this;
}

uint64_t SparsePrefixSum::operator[] (uint64_t const index) const
{
  if (index < std::size(prefix_sum_bits_))
  {
    return prefix_sum_select_1_(index);
  }
  return std::size(prefix_sum_bits_);
}

std::ostream& operator<< (std::ostream& out, SparsePrefixSum const& sparse_prefix_sum)
{
  for (uint64_t i {1}, offset {}; offset != (std::size(sparse_prefix_sum.prefix_sum_bits_) - 1); ++i)
  {
    offset = sparse_prefix_sum.prefix_sum_select_1_(i);
    out << offset << ((offset != (std::size(sparse_prefix_sum.prefix_sum_bits_) - 1)) ? " " : "\n");
  }
  out << "space:\n";
  out << ProperSizeRepresentation(sdsl::size_in_bytes(sparse_prefix_sum.prefix_sum_bits_)) << "B\n";
  return out;
}

uint64_t SparsePrefixSum::Serialize
(
  std::ostream& out,
  std::shared_ptr<SpaceNode> parent,
  std::string const name
)
{
  uint64_t size_in_bytes {};
  if (!parent)
  {
    sdsl::serialize(prefix_sum_bits_, out);
    sdsl::serialize(prefix_sum_select_1_, out);
  }
  else
  {
    auto node {std::make_shared<SpaceNode>(name)};
    node->AddLeaf("prefix_sum_bits_", sdsl::serialize(prefix_sum_bits_, out));
    node->AddLeaf("prefix_sum_select_1_", sdsl::serialize(prefix_sum_select_1_, out));
    parent->AddChild(node);
    size_in_bytes = node->GetSizeInBytes();
  }
  return size_in_bytes;
}

void SparsePrefixSum::Load (std::istream& in)
{
  prefix_sum_bits_.load(in);
  prefix_sum_select_1_.load(in);
  prefix_sum_select_1_.set_vector(&prefix_sum_bits_);
  return;
}
}
