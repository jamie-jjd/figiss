#pragma once

#include "utility.h"

namespace figiss
{
class IntegerAlphabet
{
public:

  IntegerAlphabet () = default;
  IntegerAlphabet (sdsl::int_vector<> const& text);
  IntegerAlphabet (IntegerAlphabet const&);
  IntegerAlphabet (IntegerAlphabet&&);
  IntegerAlphabet& operator= (IntegerAlphabet const&);
  IntegerAlphabet& operator= (IntegerAlphabet&&);
  ~IntegerAlphabet () = default;

  void Swap (IntegerAlphabet&);

  uint64_t Serialize
  (
    std::ostream& out,
    std::shared_ptr<SpaceNode> parent = nullptr,
    std::string const name = ""
  );
  void Load (std::istream& in);

  inline uint8_t GetAlphabetWidth () const
  {
    return sdsl::bits::hi(std::size(effective_alphabet_bits_) - 1) + 1;
  }

  inline uint8_t GetEffectiveAlphabetWidth () const
  {
    auto size {std::size(effective_alphabet_bits_)};
    if (size)
    {
      return sdsl::bits::hi(effective_alphabet_rank_1_(std::size(effective_alphabet_bits_)) - 1) + 1;
    }
    return size;
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
  auto alphabet_width {sdsl::bits::hi(*std::rbegin(effective_alphabet)) + 1};
  sdsl::int_vector<> sorted_integers(std::size(effective_alphabet), 0, alphabet_width);
  std::copy(std::begin(effective_alphabet), std::end(effective_alphabet), std::begin(sorted_integers));
  effective_alphabet_bits_ = decltype(effective_alphabet_bits_)(std::begin(sorted_integers), std::end(sorted_integers));
  effective_alphabet_rank_1_ = decltype(effective_alphabet_rank_1_)(&effective_alphabet_bits_);
}

IntegerAlphabet::IntegerAlphabet (IntegerAlphabet const& integer_alphabet)
{
  if (this != &integer_alphabet)
  {
    effective_alphabet_bits_ = integer_alphabet.effective_alphabet_bits_;
    effective_alphabet_rank_1_.set_vector(&effective_alphabet_bits_);
  }
}

IntegerAlphabet::IntegerAlphabet (IntegerAlphabet&& integer_alphabet)
{
  if (this != &integer_alphabet)
  {
    this->Swap(integer_alphabet);
  }
}

IntegerAlphabet& IntegerAlphabet::operator= (IntegerAlphabet const& integer_alphabet)
{
  if (this != &integer_alphabet)
  {
    IntegerAlphabet temp {integer_alphabet};
    this->Swap(temp);
  }
  return *this;
}

IntegerAlphabet& IntegerAlphabet::operator= (IntegerAlphabet&& integer_alphabet)
{
  if (this != &integer_alphabet)
  {
    IntegerAlphabet temp {std::move(integer_alphabet)};
    this->Swap(temp);
  }
  return *this;
}

void IntegerAlphabet::Swap (IntegerAlphabet& integer_alphabet)
{
  effective_alphabet_bits_.swap(integer_alphabet.effective_alphabet_bits_);
  effective_alphabet_rank_1_.swap(integer_alphabet.effective_alphabet_rank_1_);
  effective_alphabet_rank_1_.set_vector(&effective_alphabet_bits_);
  integer_alphabet.effective_alphabet_rank_1_.set_vector(&integer_alphabet.effective_alphabet_bits_);
  return;
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
    sdsl::serialize(effective_alphabet_bits_, out);
    sdsl::serialize(effective_alphabet_rank_1_, out);
  }
  else
  {
    auto node {std::make_shared<SpaceNode>(name)};
    node->AddLeaf("effective_alphabet_bits_", sdsl::serialize(effective_alphabet_bits_, out));
    node->AddLeaf("effective_alphabet_rank_1_", sdsl::serialize(effective_alphabet_rank_1_, out));
    parent->AddChild(node);
    size_in_bytes = node->GetSizeInBytes();
  }
  return size_in_bytes;
}

void IntegerAlphabet::Load (std::istream& in)
{
  effective_alphabet_bits_.load(in);
  effective_alphabet_rank_1_.load(in);
  effective_alphabet_rank_1_.set_vector(&effective_alphabet_bits_);
  return;
}

std::ostream& operator<< (std::ostream& out, IntegerAlphabet const& integer_alphabet)
{
  {
    out << "alphabet_width:";
    out << static_cast<uint64_t>(integer_alphabet.GetAlphabetWidth()) << "\n";
    out << "effective_alphabet_width:";
    out << static_cast<uint64_t>(integer_alphabet.GetEffectiveAlphabetWidth()) << "\n";
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
    out << "effective_alphabet_bits_:" << ProperSizeRepresentation(sdsl::size_in_bytes(integer_alphabet.effective_alphabet_bits_)) << "B\n";
    out << "effective_alphabet_rank_1_:" << ProperSizeRepresentation(sdsl::size_in_bytes(integer_alphabet.effective_alphabet_rank_1_)) << "B\n";
  }
  return out;
}
}