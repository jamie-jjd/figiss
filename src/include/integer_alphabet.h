#pragma once

#include "utility.h"

namespace figiss
{
class IntegerAlphabet
{
public:

  IntegerAlphabet () = default;
  IntegerAlphabet (sdsl::int_vector<> const& text);
  IntegerAlphabet& operator= (IntegerAlphabet&&);

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

  inline uint64_t GetRank (uint64_t const integer) const
  {
    if ((integer < std::size(effective_alphabet_bits_)) && effective_alphabet_bits_[integer])
    {
      return effective_alphabet_rank_1_(integer);
    }
    return 0;
  }

  friend std::ostream& operator<< (std::ostream& out, IntegerAlphabet const& integer_alphabet);

private:

  uint8_t effective_alphabet_width_;
  sdsl::sd_vector<> effective_alphabet_bits_;
  sdsl::sd_vector<>::rank_1_type effective_alphabet_rank_1_;

};

IntegerAlphabet::IntegerAlphabet (sdsl::int_vector<> const& text)
{
  std::set<uint8_t> effective_alphabet;
  for (auto const integer : text)
  {
    effective_alphabet.insert(integer);
  }
  effective_alphabet_width_ = sdsl::bits::hi(*std::rbegin(effective_alphabet)) + 1;
  sdsl::int_vector<> sorted_integers
  (
    std::size(effective_alphabet),
    0,
    sdsl::bits::hi(*effective_alphabet.rbegin()) + 1
  );
  std::copy(std::begin(effective_alphabet), std::end(effective_alphabet), std::begin(sorted_integers));
  effective_alphabet_bits_ = decltype(effective_alphabet_bits_)(std::begin(sorted_integers), std::end(sorted_integers));
  effective_alphabet_rank_1_ = decltype(effective_alphabet_rank_1_)(&effective_alphabet_bits_);
}

IntegerAlphabet& IntegerAlphabet::operator= (IntegerAlphabet&& integer_alphabet)
{
  if (this !=& integer_alphabet)
  {
    effective_alphabet_width_ = std::move(integer_alphabet.effective_alphabet_width_);
    effective_alphabet_bits_ = std::move(integer_alphabet.effective_alphabet_bits_);
    effective_alphabet_rank_1_ = std::move(integer_alphabet.effective_alphabet_rank_1_);
    effective_alphabet_rank_1_.set_vector(&effective_alphabet_bits_);
  }
  return *this;
}

uint64_t IntegerAlphabet::Serialize
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

void IntegerAlphabet::Load (std::istream& in)
{
  sdsl::read_member(effective_alphabet_width_, in);
  effective_alphabet_bits_.load(in);
  effective_alphabet_rank_1_.load(in);
  effective_alphabet_rank_1_.set_vector(&effective_alphabet_bits_);
  return;
}

std::ostream& operator<< (std::ostream& out, IntegerAlphabet const& integer_alphabet)
{
  {
    out << "effective_alphabet_width_:";
    out << static_cast<uint64_t>(integer_alphabet.effective_alphabet_width_) << "\n";
    out << "integer_alphabet:\n";
    auto &bits {integer_alphabet.effective_alphabet_bits_};
    auto select_1 {decltype(integer_alphabet.effective_alphabet_bits_)::select_1_type(&bits)};
    auto number_of_ones {integer_alphabet.effective_alphabet_rank_1_(std::size(bits))};
    for (uint64_t i {}; i != number_of_ones; ++i)
    {
      out << i << ":" << select_1(i + 1) << "\n";
    }
  }
  {
    out << "effective_alphabet_width_:" << ProperSizeRepresentation(sizeof(integer_alphabet.effective_alphabet_width_)) << "B\n";
    out << "effective_alphabet_bits_:" << ProperSizeRepresentation(sdsl::size_in_bytes(integer_alphabet.effective_alphabet_bits_)) << "B\n";
    out << "effective_alphabet_rank_1_:" << ProperSizeRepresentation(sdsl::size_in_bytes(integer_alphabet.effective_alphabet_rank_1_)) << "B\n";
  }
  return out;
}
}
