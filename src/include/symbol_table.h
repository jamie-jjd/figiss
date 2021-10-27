#pragma once

#include "utility.h"

namespace figiss
{
class SymbolTable
{
public:

  SymbolTable () = default;
  SymbolTable (SymbolTable const&) = delete;
  SymbolTable (SymbolTable&&);
  SymbolTable (sdsl::int_vector<8> const &byte_text);
  SymbolTable& operator= (SymbolTable const&) = delete;
  SymbolTable& operator= (SymbolTable &&);

  void Swap (SymbolTable&);

  uint64_t Serialize
  (
    std::ostream &out,
    std::shared_ptr<SpaceNode> parent = nullptr,
    std::string const name = ""
  );
  void Load (std::istream &in);

  inline uint8_t GetEffectiveAlphabetWidth () const
  {
    return effective_alphabet_width_;
  }

  uint64_t operator[] (uint64_t const byte) const;

  friend std::ostream& operator<< (std::ostream &out, SymbolTable const &symbol_table);

private:

  uint8_t effective_alphabet_width_;
  sdsl::bit_vector alphabet_bits_;
  sdsl::bit_vector::rank_1_type alphabet_rank_1_;

};

SymbolTable::SymbolTable (SymbolTable&& symbol_table)
{
  if (this != &symbol_table)
  {
    this->Swap(symbol_table);
  }
}

SymbolTable::SymbolTable (sdsl::int_vector<8> const &byte_text)
{
  alphabet_bits_.resize(*std::max_element(std::begin(byte_text), std::end(byte_text)) + 1);
  sdsl::util::set_to_value(alphabet_bits_, 0);
  for (auto const byte : byte_text)
  {
    alphabet_bits_[byte] = 1;
  }
  alphabet_rank_1_ = decltype(alphabet_rank_1_)(&alphabet_bits_);
  effective_alphabet_width_ = sdsl::bits::hi(alphabet_rank_1_(std::size(alphabet_bits_)) - 1) + 1;
}

SymbolTable& SymbolTable::operator= (SymbolTable &&symbol_table)
{
  if (this != &symbol_table)
  {
    SymbolTable temp {std::move(symbol_table)};
    this->Swap(temp);
  }
  return *this;
}

void SymbolTable::Swap (SymbolTable& symbol_table)
{
  if (this != &symbol_table)
  {
    std::swap(effective_alphabet_width_, symbol_table.effective_alphabet_width_);
    alphabet_bits_.swap(symbol_table.alphabet_bits_);
    alphabet_rank_1_.swap(symbol_table.alphabet_rank_1_);
    alphabet_rank_1_.set_vector(&alphabet_bits_);
    symbol_table.alphabet_rank_1_.set_vector(&symbol_table.alphabet_bits_);
  }
  return;
}

uint64_t SymbolTable::Serialize
(
  std::ostream &out,
  std::shared_ptr<SpaceNode> parent,
  std::string const name
)
{
  uint64_t size_in_bytes {};
  if (!parent)
  {
    sdsl::write_member(effective_alphabet_width_, out);
    sdsl::serialize(alphabet_bits_, out);
    sdsl::serialize(alphabet_rank_1_, out);
  }
  else
  {
    auto node {std::make_shared<SpaceNode>(name)};
    node->AddLeaf("effective_alphabet_width_", sdsl::write_member(effective_alphabet_width_, out));
    node->AddLeaf("alphabet_bits_", sdsl::serialize(alphabet_bits_, out));
    node->AddLeaf("alphabet_rank_1_", sdsl::serialize(alphabet_rank_1_, out));
    parent->AddChild(node);
    size_in_bytes = node->GetSizeInBytes();
  }
  return size_in_bytes;
}

void SymbolTable::Load (std::istream &in)
{
  sdsl::read_member(effective_alphabet_width_, in);
  alphabet_bits_.load(in);
  alphabet_rank_1_.load(in);
  alphabet_rank_1_.set_vector(&alphabet_bits_);
  return;
}

uint64_t SymbolTable::operator[] (uint64_t const byte) const
{
  if ((byte < std::size(alphabet_bits_)) && alphabet_bits_[byte])
  {
    return alphabet_rank_1_(byte);
  }
  return 0;
}

std::ostream& operator<< (std::ostream &out, SymbolTable const &symbol_table)
{
  {
    out << "value:\n";
    out << "effective_alphabet_width_:\n";
    out << static_cast<uint64_t>(symbol_table.effective_alphabet_width_) << "\n";
    out << "byte alphabet:\n";
    out << 0 << ":" << 0 << "\n";
    for (uint64_t i {1}; i != std::size(symbol_table.alphabet_bits_); ++i)
    {
      auto symbol {symbol_table[i]};
      if (symbol != 0)
      {
        out << symbol << ":" << i << "\n";
      }
    }
  }
  {
    out << "space:\n";
    out << "effective_alphabet_width_: " << sizeof(symbol_table.effective_alphabet_width_) << "B\n";
    out << "alphabet_bits_: " << ProperSizeRepresentation(sdsl::size_in_bytes(symbol_table.alphabet_bits_)) << "B\n";
    out << "alphabet_rank_1_: " << ProperSizeRepresentation(sdsl::size_in_bytes(symbol_table.alphabet_rank_1_)) << "B\n";
  }
  return out;
}
}
