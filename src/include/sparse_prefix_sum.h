#pragma once

#include "utility.h"

namespace figiss
{
class SparsePrefixSum
{
public:

  SparsePrefixSum () = default;
  SparsePrefixSum (SparsePrefixSum const&) = delete;
  SparsePrefixSum (SparsePrefixSum&&);
  SparsePrefixSum (std::vector<uint64_t>& sorted_integers);
  SparsePrefixSum& operator= (SparsePrefixSum const&) = delete;
  SparsePrefixSum& operator= (SparsePrefixSum&&);
  ~SparsePrefixSum () = default;

  void Swap (SparsePrefixSum&);

  uint64_t Serialize
  (
    std::ostream &out,
    std::shared_ptr<SpaceNode> parent = nullptr,
    std::string const name = ""
  );
  void Load (std::istream &in);

  uint64_t operator[] (uint64_t const prefix_sum) const;
  uint64_t Select (uint64_t const i) const;

  std::pair<uint64_t, uint64_t> RankRange
  (
    uint64_t const smallest_prefix_sum,
    uint64_t const largest_prefix_sum
  ) const;

  friend std::ostream& operator<< (std::ostream &out, SparsePrefixSum const& sparse_prefix_sum);

private:

  uint64_t size_;
  sdsl::sd_vector<> bits_;
  sdsl::sd_vector<>::rank_1_type rank_1_;
  sdsl::sd_vector<>::select_1_type select_1_;

};

SparsePrefixSum::SparsePrefixSum (SparsePrefixSum&& sparse_prefix_sum)
{
  if (this != &sparse_prefix_sum)
  {
    size_ = std::move(sparse_prefix_sum.size_);
    bits_ = std::move(sparse_prefix_sum.bits_);
    rank_1_ = std::move(sparse_prefix_sum.rank_1_);
    rank_1_.set_vector(&bits_);
    select_1_ = std::move(sparse_prefix_sum.select_1_);
    select_1_.set_vector(&bits_);
  }
}

SparsePrefixSum::SparsePrefixSum (std::vector<uint64_t>& sorted_integers)
{
  size_ = std::size(sorted_integers);
  bits_ = decltype(bits_)(std::begin(sorted_integers), std::end(sorted_integers));
  rank_1_ = decltype(rank_1_)(&bits_);
  select_1_ = decltype(select_1_)(&bits_);
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
    std::swap(size_, sparse_prefix_sum.size_);
    bits_.swap(sparse_prefix_sum.bits_);
    rank_1_.swap(sparse_prefix_sum.rank_1_);
    rank_1_.set_vector(&bits_);
    sparse_prefix_sum.rank_1_.set_vector(&sparse_prefix_sum.bits_);
    select_1_.swap(sparse_prefix_sum.select_1_);
    select_1_.set_vector(&bits_);
    sparse_prefix_sum.select_1_.set_vector(&sparse_prefix_sum.bits_);
  }
  return;
}

uint64_t SparsePrefixSum::Serialize
(
  std::ostream &out,
  std::shared_ptr<SpaceNode> parent,
  std::string const name
)
{
  uint64_t size_in_bytes {};
  if (!parent)
  {
    sdsl::write_member(size_, out);
    sdsl::serialize(bits_, out);
    sdsl::serialize(rank_1_, out);
    sdsl::serialize(select_1_, out);
  }
  else
  {
    auto node {std::make_shared<SpaceNode>(name)};
    node->AddLeaf("size_", sdsl::write_member(size_, out));
    node->AddLeaf("bits_", sdsl::serialize(bits_, out));
    node->AddLeaf("rank_1_", sdsl::serialize(rank_1_, out));
    node->AddLeaf("select_1_", sdsl::serialize(select_1_, out));
    parent->AddChild(node);
    size_in_bytes = node->GetSizeInBytes();
  }
  return size_in_bytes;
}

void SparsePrefixSum::Load (std::istream &in)
{
  sdsl::read_member(size_, in);
  bits_.load(in);
  rank_1_.load(in);
  rank_1_.set_vector(&bits_);
  select_1_.load(in);
  select_1_.set_vector(&bits_);
  return;
}

uint64_t SparsePrefixSum::operator[] (uint64_t const prefix_sum) const
{
  if ((prefix_sum < std::size(bits_)) && bits_[prefix_sum])
  {
    return rank_1_(prefix_sum);
  }
  return 0;
}

uint64_t SparsePrefixSum::Select (uint64_t const i) const
{
  if (i <= size_)
  {
    return select_1_(i);
  }
  return (std::size(bits_) - 1);
}

std::pair<uint64_t, uint64_t> SparsePrefixSum::RankRange
(
  uint64_t const smallest_prefix_sum,
  uint64_t const largest_prefix_sum
) const
{
  std::pair<uint64_t, uint64_t> rank_range {1, 0};
  if (smallest_prefix_sum < std::size(bits_))
  {
    std::get<0>(rank_range) = rank_1_(smallest_prefix_sum);
    if (largest_prefix_sum < std::size(bits_))
    {
      std::get<1>(rank_range) =
      (
        rank_1_(largest_prefix_sum)
        - (bits_[largest_prefix_sum] == 0)
      );
    }
    else
    {
      std::get<1>(rank_range) = size_ - 1;
    }
  }
  return rank_range;
}

std::ostream& operator<< (std::ostream &out, SparsePrefixSum const &sparse_prefix_sum)
{
  out << "ranks:\n";
  for (uint64_t i {}; i != sparse_prefix_sum.size_; ++i)
  {
    out << i << ":" << sparse_prefix_sum.select_1_(i + 1) << "\n";
  }
  out << "space:\n";
  out << "size_:" << ProperSizeRepresentation(sizeof(sparse_prefix_sum.size_)) << "B\n";
  out << "bits_:" << ProperSizeRepresentation(sdsl::size_in_bytes(sparse_prefix_sum.bits_)) << "B\n";
  out << "rank_1_:" << ProperSizeRepresentation(sdsl::size_in_bytes(sparse_prefix_sum.rank_1_)) << "B\n";
  out << "select_1_:" << ProperSizeRepresentation(sdsl::size_in_bytes(sparse_prefix_sum.select_1_)) << "B\n";
  return out;
}
}
