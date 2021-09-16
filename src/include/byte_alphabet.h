#pragma once

#include "utility.h"

namespace figiss
{
class ByteAlphabet
{
public:

  ByteAlphabet () = default;
  ByteAlphabet (sdsl::int_vector<8> const& byte_text);
  ByteAlphabet (ByteAlphabet const&) = delete;
  ByteAlphabet (ByteAlphabet&&);
  ByteAlphabet& operator= (ByteAlphabet const&) = delete;
  ByteAlphabet& operator= (ByteAlphabet&&);
  ~ByteAlphabet () = default;

  void Swap (ByteAlphabet&);

  uint64_t Serialize
  (
    std::ostream& out,
    std::shared_ptr<SpaceNode> parent = nullptr,
    std::string const name = ""
  );
  void Load (std::istream& in);

  inline uint8_t GetEffectiveAlphabetWidth () const
  {
    return effective_alphabet_width_;
  }

  inline uint64_t GetRank (uint64_t const byte) const
  {
    if ((byte < std::size(effective_alphabet_bits_)) && effective_alphabet_bits_[byte])
    {
      return effective_alphabet_rank_1_(byte);
    }
    return 0;
  }

  friend std::ostream& operator<< (std::ostream& out, ByteAlphabet const& byte_alphabet);

private:

  uint8_t effective_alphabet_width_;
  sdsl::bit_vector effective_alphabet_bits_;
  sdsl::bit_vector::rank_1_type effective_alphabet_rank_1_;

};

ByteAlphabet::ByteAlphabet (sdsl::int_vector<8> const& byte_text)
{
  std::set<uint8_t> effective_alphabet;
  for (auto const byte : byte_text)
  {
    effective_alphabet.insert(byte);
  }
  effective_alphabet_bits_.resize(*effective_alphabet.rbegin() + 1);
  sdsl::util::set_to_value(effective_alphabet_bits_, 0);
  for (auto const byte : effective_alphabet)
  {
    effective_alphabet_bits_[byte] = 1;
  }
  effective_alphabet_rank_1_ = decltype(effective_alphabet_rank_1_)(&effective_alphabet_bits_);
  effective_alphabet_width_ = sdsl::bits::hi(std::size(effective_alphabet) - 1) + 1;
}

ByteAlphabet::ByteAlphabet (ByteAlphabet&& byte_alphabet)
{
  if (this != &byte_alphabet)
  {
    effective_alphabet_width_ = std::move(byte_alphabet.effective_alphabet_width_);
    effective_alphabet_bits_ = std::move(byte_alphabet.effective_alphabet_bits_);
    effective_alphabet_rank_1_ = std::move(byte_alphabet.effective_alphabet_rank_1_);
    effective_alphabet_rank_1_.set_vector(&effective_alphabet_bits_);
  }
}

ByteAlphabet& ByteAlphabet::operator= (ByteAlphabet&& byte_alphabet)
{
  if (this != &byte_alphabet)
  {
    ByteAlphabet temp {std::move(byte_alphabet)};
    this->Swap(temp);
  }
  return *this;
}

void ByteAlphabet::Swap (ByteAlphabet& byte_alphabet)
{
  std::swap(effective_alphabet_width_, byte_alphabet.effective_alphabet_width_);
  effective_alphabet_bits_.swap(byte_alphabet.effective_alphabet_bits_);
  effective_alphabet_rank_1_.set_vector(&effective_alphabet_bits_);
  byte_alphabet.effective_alphabet_rank_1_.set_vector(&byte_alphabet.effective_alphabet_bits_);
  return;
}

uint64_t ByteAlphabet::Serialize
(
  std::ostream& out,
  std::shared_ptr<SpaceNode> parent,
  std::string const name
)
{
  uint64_t size_in_bytes {};
  if (!parent)
  {
    sdsl::write_member(effective_alphabet_width_, out);
    sdsl::serialize(effective_alphabet_bits_, out);
    sdsl::serialize(effective_alphabet_rank_1_, out);
  }
  else
  {
    auto node {std::make_shared<SpaceNode>(name)};
    node->AddLeaf("effective_alphabet_width_", sdsl::write_member(effective_alphabet_width_, out));
    node->AddLeaf("effective_alphabet_bits_", sdsl::serialize(effective_alphabet_bits_, out));
    node->AddLeaf("effective_alphabet_rank_1_", sdsl::serialize(effective_alphabet_rank_1_, out));
    parent->AddChild(node);
    size_in_bytes = node->GetSizeInBytes();
  }
  return size_in_bytes;
}

void ByteAlphabet::Load (std::istream& in)
{
  sdsl::read_member(effective_alphabet_width_, in);
  effective_alphabet_bits_.load(in);
  effective_alphabet_rank_1_.load(in);
  effective_alphabet_rank_1_.set_vector(&effective_alphabet_bits_);
  return;
}

std::ostream& operator<< (std::ostream& out, ByteAlphabet const& byte_alphabet)
{
  {
    out << "effective_alphabet_width_:";
    out << static_cast<uint64_t>(byte_alphabet.effective_alphabet_width_) << "\n";
    out << "byte_alphabet:\n";
    for (uint64_t byte {}; byte != std::size(byte_alphabet.effective_alphabet_bits_); ++byte)
    {
      if (byte_alphabet.effective_alphabet_bits_[byte])
      {
        out << byte_alphabet.GetRank(byte) << ":" << byte  << "(" << static_cast<char>(byte) << ")" << "\n";
      }
    }
  }
  {
    out << "effective_alphabet_width_:" << ProperSizeRepresentation(sizeof(byte_alphabet.effective_alphabet_width_)) << "B\n";
    out << "effective_alphabet_bits_:" << ProperSizeRepresentation(sdsl::size_in_bytes(byte_alphabet.effective_alphabet_bits_)) << "B\n";
    out << "effective_alphabet_rank_1_:" << ProperSizeRepresentation(sdsl::size_in_bytes(byte_alphabet.effective_alphabet_rank_1_)) << "B\n";
  }
  return out;
}
}
