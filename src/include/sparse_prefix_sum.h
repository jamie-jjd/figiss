#pragma once

#include "utility.h"

namespace figiss
{
class SparsePrefixSum
{
public:

  SparsePrefixSum () = default;
  SparsePrefixSum (std::deque<uint64_t> const& prefix_sum);
  SparsePrefixSum (SparsePrefixSum const&) = delete;
  SparsePrefixSum (SparsePrefixSum&&);
  SparsePrefixSum& operator= (SparsePrefixSum const&) = delete;
  SparsePrefixSum& operator= (SparsePrefixSum&&);

  void Swap (SparsePrefixSum&);

  uint64_t Serialize
  (
    std::ostream& out,
    std::shared_ptr<SpaceNode> parent = nullptr,
    std::string const name = ""
  );
  void Load (std::istream& in);

  inline uint64_t operator[] (uint64_t const index) const
  {
    if (index < std::size(prefix_sum_bits_))
    {
      return prefix_sum_select_1_(index + 1);
    }
    return std::size(prefix_sum_bits_);
  }

  inline uint64_t Predecessor (uint64_t const index) const
  {
    if (index < std::size(prefix_sum_bits_))
    {
      return prefix_sum_rank_1_(index);
    }
    return prefix_sum_rank_1_(std::size(prefix_sum_bits_));
  }

  friend std::ostream& operator<< (std::ostream& out, SparsePrefixSum const& sparse_prefix_sum);

private:

  sdsl::sd_vector<> prefix_sum_bits_;
  sdsl::sd_vector<>::select_1_type prefix_sum_select_1_;
  sdsl::sd_vector<>::rank_1_type prefix_sum_rank_1_;

};

SparsePrefixSum::SparsePrefixSum (std::deque<uint64_t> const& prefix_sum)
{
  prefix_sum_bits_ = decltype(prefix_sum_bits_)(std::begin(prefix_sum), std::end(prefix_sum));
  prefix_sum_select_1_.set_vector(&prefix_sum_bits_);
  prefix_sum_rank_1_.set_vector(&prefix_sum_bits_);
}

SparsePrefixSum::SparsePrefixSum (SparsePrefixSum&& sparse_prefix_sum)
{
  if (this != &sparse_prefix_sum)
  {
    this->Swap(sparse_prefix_sum);
  }
}

SparsePrefixSum& SparsePrefixSum::operator= (SparsePrefixSum&& sparse_prefix_sum)
{
  if (this != &sparse_prefix_sum)
  {
    SparsePrefixSum temp {std::move(sparse_prefix_sum)};
    this->Swap(temp);
  }
  return *this;
}

void SparsePrefixSum::Swap (SparsePrefixSum& sparse_prefix_sum)
{
  if (this != &sparse_prefix_sum)
  {
    prefix_sum_bits_.swap(sparse_prefix_sum.prefix_sum_bits_);
    prefix_sum_select_1_.set_vector(&prefix_sum_bits_);
    sparse_prefix_sum.prefix_sum_select_1_.set_vector(&sparse_prefix_sum.prefix_sum_bits_);
    prefix_sum_rank_1_.set_vector(&prefix_sum_bits_);
    sparse_prefix_sum.prefix_sum_rank_1_.set_vector(&sparse_prefix_sum.prefix_sum_bits_);
  }
  return;
}

std::ostream& operator<< (std::ostream& out, SparsePrefixSum const& sparse_prefix_sum)
{
  auto number_of_ones {sparse_prefix_sum.prefix_sum_rank_1_(std::size(sparse_prefix_sum.prefix_sum_bits_))};
  for (uint64_t i {}; i != number_of_ones; ++i)
  {
    out << sparse_prefix_sum.prefix_sum_select_1_(i + 1) << ((i != (number_of_ones - 1)) ? " " : "\n");
  }
  {
    out << "prefix_sum_bits_:" << ProperSizeRepresentation(sdsl::size_in_bytes(sparse_prefix_sum.prefix_sum_bits_)) << "B\n";
  }
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
  }
  else
  {
    auto node {std::make_shared<SpaceNode>(name)};
    node->AddLeaf("prefix_sum_bits_", sdsl::serialize(prefix_sum_bits_, out));
    parent->AddChild(node);
    size_in_bytes = node->GetSizeInBytes();
  }
  return size_in_bytes;
}

void SparsePrefixSum::Load (std::istream& in)
{
  prefix_sum_bits_.load(in);
  prefix_sum_select_1_.set_vector(&prefix_sum_bits_);
  prefix_sum_rank_1_.set_vector(&prefix_sum_bits_);
  return;
}
}
