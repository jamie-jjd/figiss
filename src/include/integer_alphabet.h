#pragma once

#include "utility.h"

namespace figiss
{
class IntegerAlphabet
{
public:

  IntegerAlphabet () = default;
  IntegerAlphabet (std::set<uint64_t> const& alphabet);
  IntegerAlphabet (IntegerAlphabet const&) = delete;
  IntegerAlphabet (IntegerAlphabet&&);
  IntegerAlphabet& operator= (IntegerAlphabet const&) = delete;
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
  uint64_t operator[] (uint64_t const integer) const;
  uint64_t Predecessor (uint64_t const integer) const;
  uint64_t Successor (uint64_t const integer) const;
  uint64_t Select (uint64_t const index) const;

  inline uint8_t GetAlphabetWidth () const
  {
    return sdsl::bits::hi(std::size(effective_alphabet_bits_) - 1) + 1;
  }

  inline uint8_t GetEffectiveAlphabetWidth () const
  {
    return sdsl::bits::hi(effective_alphabet_rank_1_(std::size(effective_alphabet_bits_))) + 1;
  }

  friend std::ostream& operator<< (std::ostream& out, IntegerAlphabet const& integer_alphabet);

private:

  sdsl::sd_vector<> effective_alphabet_bits_;
  sdsl::sd_vector<>::rank_1_type effective_alphabet_rank_1_;
  sdsl::sd_vector<>::select_1_type effective_alphabet_select_1_;

};

IntegerAlphabet::IntegerAlphabet (std::set<uint64_t> const& alphabet)
{
  if (alphabet.find(0) != alphabet.end())
  {
    throw std::runtime_error("alphabet contains 0");
  }
  auto alphabet_width {sdsl::bits::hi(*std::rbegin(alphabet)) + 1};
  sdsl::int_vector<> sorted_integers(std::size(alphabet), 0, alphabet_width);
  std::copy(std::begin(alphabet), std::end(alphabet), std::begin(sorted_integers));
  effective_alphabet_bits_ = decltype(effective_alphabet_bits_)(std::begin(sorted_integers), std::end(sorted_integers));
  effective_alphabet_rank_1_ = decltype(effective_alphabet_rank_1_)(&effective_alphabet_bits_);
  effective_alphabet_select_1_ = decltype(effective_alphabet_select_1_)(&effective_alphabet_bits_);
}

// IntegerAlphabet::IntegerAlphabet (IntegerAlphabet const& integer_alphabet)
// {
//   if (this != &integer_alphabet)
//   {
//     effective_alphabet_bits_ = integer_alphabet.effective_alphabet_bits_;
//     effective_alphabet_rank_1_.set_vector(&effective_alphabet_bits_);
//   }
// }

IntegerAlphabet::IntegerAlphabet (IntegerAlphabet&& integer_alphabet)
{
  if (this != &integer_alphabet)
  {
    this->Swap(integer_alphabet);
  }
}

// IntegerAlphabet& IntegerAlphabet::operator= (IntegerAlphabet const& integer_alphabet)
// {
//   if (this != &integer_alphabet)
//   {
//     IntegerAlphabet temp {integer_alphabet};
//     this->Swap(temp);
//   }
//   return *this;
// }

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
  effective_alphabet_rank_1_.set_vector(&effective_alphabet_bits_);
  integer_alphabet.effective_alphabet_rank_1_.set_vector(&integer_alphabet.effective_alphabet_bits_);
  effective_alphabet_select_1_.set_vector(&effective_alphabet_bits_);
  integer_alphabet.effective_alphabet_select_1_.set_vector(&integer_alphabet.effective_alphabet_bits_);
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
  }
  else
  {
    auto node {std::make_shared<SpaceNode>(name)};
    node->AddLeaf("effective_alphabet_bits_", sdsl::serialize(effective_alphabet_bits_, out));
    parent->AddChild(node);
    size_in_bytes = node->GetSizeInBytes();
  }
  return size_in_bytes;
}

void IntegerAlphabet::Load (std::istream& in)
{
  effective_alphabet_bits_.load(in);
  effective_alphabet_rank_1_.set_vector(&effective_alphabet_bits_);
  effective_alphabet_select_1_.set_vector(&effective_alphabet_bits_);
  return;
}

uint64_t IntegerAlphabet::operator[] (uint64_t const integer) const
{
  if ((integer < std::size(effective_alphabet_bits_)) && effective_alphabet_bits_[integer])
  {
    return (effective_alphabet_rank_1_(integer) + 1);
  }
  return 0;
}

uint64_t IntegerAlphabet::Predecessor (uint64_t const integer) const
{
  if (integer < std::size(effective_alphabet_bits_))
  {
    return effective_alphabet_rank_1_(integer);
  }
  return effective_alphabet_rank_1_(std::size(effective_alphabet_bits_));
}

uint64_t IntegerAlphabet::Successor (uint64_t const integer) const
{
  if (integer < std::size(effective_alphabet_bits_) - 1)
  {
    return (effective_alphabet_rank_1_(integer + 1) + 1);
  }
  return 0;
}

uint64_t IntegerAlphabet::Select (uint64_t const index) const
{
  if ((index - 1) < effective_alphabet_rank_1_(std::size(effective_alphabet_bits_)))
  {
    return effective_alphabet_select_1_(index);
  }
  return 0;
}

std::ostream& operator<< (std::ostream& out, IntegerAlphabet const& integer_alphabet)
{
  {
    out << "alphabet_width:";
    out << static_cast<uint64_t>(integer_alphabet.GetAlphabetWidth()) << "\n";
    out << "effective_alphabet_width:";
    out << static_cast<uint64_t>(integer_alphabet.GetEffectiveAlphabetWidth()) << "\n";
    out << "integer_alphabet:\n";
    auto& bits {integer_alphabet.effective_alphabet_bits_};
    auto number_of_ones {integer_alphabet.effective_alphabet_rank_1_(std::size(bits))};
    for (uint64_t i {}; i != number_of_ones; ++i)
    {
      out << (i + 1) << ":" << integer_alphabet.Select(i + 1) << "\n";
    }
  }
  {
    out << "effective_alphabet_bits_:" << ProperSizeRepresentation(sdsl::size_in_bytes(integer_alphabet.effective_alphabet_bits_)) << "B\n";
  }
  return out;
}
}
