#pragma once

#include "utility.h"

namespace gciis
{
class SymbolBucketOffsets
{
public:

  SymbolBucketOffsets () = default;
  SymbolBucketOffsets (sdsl::int_vector<> const &offsets);
  SymbolBucketOffsets& operator= (SymbolBucketOffsets &&);

  uint64_t Serialize
  (
    std::ostream &out,
    std::shared_ptr<SpaceNode> parent = nullptr,
    std::string const name = ""
  );
  void Load (std::istream &in);

  uint64_t operator[] (uint64_t const symbol) const;

  friend std::ostream& operator<< (std::ostream &out, SymbolBucketOffsets const &symbol_bucket_offsets);

private:

  sdsl::sd_vector<> offset_bits_;
  sdsl::sd_vector<>::select_1_type offset_select_1_;

};

SymbolBucketOffsets::SymbolBucketOffsets (sdsl::int_vector<> const &offsets)
{
  offset_bits_ = decltype(offset_bits_)(std::begin(offsets), std::end(offsets));
  offset_select_1_.set_vector(&offset_bits_);
}

SymbolBucketOffsets& SymbolBucketOffsets::operator= (SymbolBucketOffsets &&symbol_bucket_offsets)
{
  if (this != &symbol_bucket_offsets)
  {
    offset_bits_ = std::move(symbol_bucket_offsets.offset_bits_);
    offset_select_1_.set_vector(&offset_bits_);
  }
  return *this;
}

uint64_t SymbolBucketOffsets::operator[] (uint64_t const symbol) const
{
  if (symbol < std::size(offset_bits_))
  {
    return offset_select_1_(symbol);
  }
  return std::size(offset_bits_);
}

std::ostream& operator<< (std::ostream &out, SymbolBucketOffsets const &symbol_bucket_offsets)
{
  for (uint64_t i {1}, offset {}; offset != (std::size(symbol_bucket_offsets.offset_bits_) - 1); ++i)
  {
    offset = symbol_bucket_offsets.offset_select_1_(i);
    out << offset << ((offset != (std::size(symbol_bucket_offsets.offset_bits_) - 1)) ? " " : "\n");
  }
  out << "space:\n";
  out << ProperSizeRepresentation(sdsl::size_in_bytes(symbol_bucket_offsets.offset_bits_)) << "B\n";
  return out;
}

uint64_t SymbolBucketOffsets::Serialize
(
  std::ostream &out,
  std::shared_ptr<SpaceNode> parent,
  std::string const name
)
{
  uint64_t size_in_bytes {};
  if (!parent)
  {
    sdsl::serialize(offset_bits_, out);
    sdsl::serialize(offset_select_1_, out);
  }
  else
  {
    auto node {std::make_shared<SpaceNode>(name)};
    node->AddLeaf("offset_bits_", sdsl::serialize(offset_bits_, out));
    node->AddLeaf("offset_select_1_", sdsl::serialize(offset_select_1_, out));
    parent->AddChild(node);
    size_in_bytes = node->GetSizeInBytes();
  }
  return size_in_bytes;
}

void SymbolBucketOffsets::Load (std::istream &in)
{
  offset_bits_.load(in);
  offset_select_1_.load(in);
  offset_select_1_.set_vector(&offset_bits_);
  return;
}
}
