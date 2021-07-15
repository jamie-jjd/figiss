#pragma once

#include "utility.h"

namespace project
{
class GrammarSymbolTable
{
public:

  GrammarSymbolTable () = default;
  GrammarSymbolTable (std::vector<uint64_t> &sorted_factor_integers);
  GrammarSymbolTable& operator= (GrammarSymbolTable &&);

  uint64_t Serialize
  (
    std::ostream &out,
    std::shared_ptr<SpaceNode> parent = nullptr,
    std::string const name = ""
  );
  void Load (std::istream &in);

  uint64_t operator[] (uint64_t const factor_integer) const;

  std::pair<uint64_t, uint64_t> SymbolRange
  (
    uint64_t const smallest_factor_integer,
    uint64_t const largest_factor_integer
  ) const;

  friend std::ostream& operator<< (std::ostream &out, GrammarSymbolTable const &symbol_table);

private:

  sdsl::sd_vector<> factor_integer_bits_;
  sdsl::sd_vector<>::rank_1_type factor_integer_rank_1_;

};

GrammarSymbolTable::GrammarSymbolTable (std::vector<uint64_t> &sorted_factor_integers)
{
  factor_integer_bits_ = decltype(factor_integer_bits_)(std::begin(sorted_factor_integers), std::end(sorted_factor_integers));
  factor_integer_rank_1_.set_vector(&factor_integer_bits_);
}

GrammarSymbolTable& GrammarSymbolTable::operator= (GrammarSymbolTable &&grammar_symbol_table)
{
  if (this != &grammar_symbol_table)
  {
    factor_integer_bits_ = std::move(grammar_symbol_table.factor_integer_bits_);
    factor_integer_rank_1_.set_vector(&factor_integer_bits_);
  }
  return *this;
}

uint64_t GrammarSymbolTable::Serialize
(
  std::ostream &out,
  std::shared_ptr<SpaceNode> parent,
  std::string const name
)
{
  uint64_t size_in_bytes {};
  if (!parent)
  {
    sdsl::serialize(factor_integer_bits_, out);
    sdsl::serialize(factor_integer_rank_1_, out);
  }
  else
  {
    auto node {std::make_shared<SpaceNode>(name)};
    node->AddLeaf("factor_integer_bits_", sdsl::serialize(factor_integer_bits_, out));
    node->AddLeaf("factor_integer_rank_1_", sdsl::serialize(factor_integer_rank_1_, out));
    parent->AddChild(node);
    size_in_bytes = node->GetSizeInBytes();
  }
  return size_in_bytes;
}

void GrammarSymbolTable::Load (std::istream &in)
{
  factor_integer_bits_.load(in);
  factor_integer_rank_1_.load(in);
  factor_integer_rank_1_.set_vector(&factor_integer_bits_);
  return;
}

uint64_t GrammarSymbolTable::operator[] (uint64_t const factor_integer) const
{
  if ((factor_integer < std::size(factor_integer_bits_)) && factor_integer_bits_[factor_integer])
  {
    return factor_integer_rank_1_(factor_integer);
  }
  return 0;
}

std::pair<uint64_t, uint64_t> GrammarSymbolTable::SymbolRange
(
  uint64_t const smallest_factor_integer,
  uint64_t const largest_factor_integer
) const
{
  std::pair<uint64_t, uint64_t> symbol_range {1, 0};
  if (smallest_factor_integer < std::size(factor_integer_bits_))
  {
    std::get<0>(symbol_range) = factor_integer_rank_1_(smallest_factor_integer);
    if (largest_factor_integer < std::size(factor_integer_bits_))
    {
      std::get<1>(symbol_range) =
      (
        factor_integer_rank_1_(largest_factor_integer)
        - (factor_integer_bits_[largest_factor_integer] == 0)
      );
    }
    else
    {
      std::get<1>(symbol_range) = factor_integer_rank_1_(std::size(factor_integer_bits_)) - 1;
    }
  }
  return symbol_range;
}

std::ostream& operator<< (std::ostream &out, GrammarSymbolTable const &grammar_symbol_table)
{
  out << "factor integers:\n";
  out << "0:0\n";
  for (uint64_t i {1}; i != std::size(grammar_symbol_table.factor_integer_bits_); ++i)
  {
    auto symbol {grammar_symbol_table[i]};
    if (symbol != 0)
    {
      out << symbol << ":" << i << "\n";
    }
  }
  out << "space:\n";
  out << "factor_integer_bits_: " << ProperSizeRepresentation(sdsl::size_in_bytes(grammar_symbol_table.factor_integer_bits_)) << "B\n";
  return out;
}
}
