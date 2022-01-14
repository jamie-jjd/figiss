#pragma once

#include <deque>
#include <filesystem>
#include <fstream>
#include <functional>
#include <random>

#include <sdsl/int_vector.hpp>

namespace figiss
{
class PatternCollection
{
public:

  PatternCollection () = default;
  PatternCollection (PatternCollection const&) = delete;
  PatternCollection (PatternCollection&&);
  PatternCollection
  (
    std::filesystem::path const& byte_text_path,
    std::filesystem::path const& pattern_parent_path,
    uint64_t const amount,
    uint64_t const length,
    bool const is_mutated = false
  );
  PatternCollection& operator= (PatternCollection const&) = delete;
  PatternCollection& operator= (PatternCollection&&);

  void Swap (PatternCollection&);

  void Serialize (std::filesystem::path const& pattern_path);
  void Load (std::filesystem::path const& pattern_path);

  auto begin () noexcept;
  auto end () noexcept;

  inline auto GetAmount () const noexcept
  {
    return amount_;
  }

  inline auto GetLength () const noexcept
  {
    return length_;
  }

  inline operator bool() const
  {
    return (amount_ != 0 && length_ != 0 );
  }

private:

  uint64_t amount_;
  uint64_t length_;
  sdsl::int_vector<8> concatenated_patterns_;

};

PatternCollection::PatternCollection (PatternCollection&& pattern_collection)
{
  if (this != &pattern_collection)
  {
    this->Swap(pattern_collection);
  }
}

PatternCollection::PatternCollection
(
  std::filesystem::path const& byte_text_path,
  std::filesystem::path const& pattern_parent_path,
  uint64_t const amount,
  uint64_t const length,
  bool const is_mutated
)
{
  sdsl::int_vector<8> byte_text;
  sdsl::load_vector_from_file(byte_text, byte_text_path);
  if (length > std::size(byte_text))
  {
    amount_ = length_ = 0;
    std::cerr << "\033[31mfailed at pattern generation: ";
    std::cerr
    << "length of pattern (" << length
    << ") is longer than "
    << "length of byte_text (" << std::size(byte_text)
    << ")\033[0m\n";
    return;
  }
  std::deque<uint64_t> begin_offsets;
  std::deque<std::pair<uint64_t, uint64_t>> swapped_offset_pairs;
  {
    std::cout << "generate "
    << amount << " random";
    if (is_mutated)
    {
      std::cout << " \"mutated\"";
    }
    std::cout << " patterns of length " << length
    << " from " << std::filesystem::canonical(byte_text_path).string() << "\n";
    amount_ = amount;
    length_ = length;
    concatenated_patterns_.resize(amount_ * (length_ + 1));
  }
  std::mt19937 engine {std::random_device{}()};
  std::uniform_int_distribution<uint64_t> dist_text {0, std::size(byte_text) - length_};
  std::uniform_int_distribution<uint64_t> dist_pattern {0, length_ - 1};
  auto it {std::begin(concatenated_patterns_)};
  for (uint64_t i {}; i != amount_; ++i)
  {
    begin_offsets.push_back(dist_text(engine));
    auto byte_text_first {std::next(std::begin(byte_text), begin_offsets.back())};
    std::copy(byte_text_first, std::next(byte_text_first, length_), it);
    if (is_mutated)
    {
      swapped_offset_pairs.push_back({dist_pattern(engine), dist_pattern(engine)});
      std::swap
      (
        *std::next(it, std::get<0>(swapped_offset_pairs.back())),
        *std::next(it, std::get<0>(swapped_offset_pairs.back()))
      );
    }
    it = std::next(it, length_);
    *it++ = 0;
  }
  {
    auto metadata_path
    {
      pattern_parent_path /
      byte_text_path.filename() /
      std::filesystem::path
      {
        std::to_string(amount_) + "/" +
        std::to_string(length_) + "/" +
        ((is_mutated) ? "mutated/" : "") +
        "metadata"
      }
    };
    if (!std::filesystem::exists(metadata_path.parent_path()))
    {
      std::filesystem::create_directories(metadata_path.parent_path());
    }
    std::ofstream out {metadata_path};
    std::cout << "write metadata of pattern collection to " << std::filesystem::canonical(metadata_path).string() << "\n";
    out
    << "# if pattern is mutated, mutation is done by swapping 2 characters on pattern\n"
    << "byte text path:" << std::filesystem::canonical(byte_text_path).string() << "\n"
    << "amount:" << amount_ << "\n"
    << "length:" << length_ << "\n";
    out << "begin offset of text[, left swapped offset of pattern, right swapped offset of pattern]:";
    for (uint64_t i {}; i != amount_; ++i)
    {
      out << begin_offsets[i];
      if (is_mutated)
      {
        out
        << " " << std::get<0>(swapped_offset_pairs[i])
        << " " << std::get<1>(swapped_offset_pairs[i]);
      }
      out << "\n";
    }
  }
  return;
}

PatternCollection& PatternCollection::operator= (PatternCollection&& pattern_collection)
{
  if (this != &pattern_collection)
  {
    PatternCollection temp {std::move(pattern_collection)};
    this->Swap(temp);
  }
  return *this;
}

void PatternCollection::Swap (PatternCollection& pattern_collection)
{
  if (this != &pattern_collection)
  {
    std::swap(amount_, pattern_collection.amount_);
    std::swap(length_, pattern_collection.length_);
    concatenated_patterns_.swap(pattern_collection.concatenated_patterns_);
  }
  return;
}

void PatternCollection::Serialize (std::filesystem::path const& pattern_path)
{
  if (std::size(concatenated_patterns_) == (amount_ * (length_ + 1)))
  {
    std::ofstream out {pattern_path, std::ios_base::out | std::ios_base::trunc};
    std::cout << "serialize pattern collection to " << std::filesystem::canonical(pattern_path).string() << "\n";
    sdsl::write_member(amount_, out);
    sdsl::write_member(length_, out);
    sdsl::serialize(concatenated_patterns_, out);
  }
  return;
}

void PatternCollection::Load (std::filesystem::path const& pattern_path)
{
  std::ifstream in {pattern_path};
  std::cout << "load pattern collection from " << std::filesystem::canonical(pattern_path).string() << "\n";
  sdsl::read_member(amount_, in);
  sdsl::read_member(length_, in);
  concatenated_patterns_.load(in);
  return;
}

auto PatternCollection::begin () noexcept
{
  return std::begin(concatenated_patterns_);
}

auto PatternCollection::end () noexcept
{
  return std::end(concatenated_patterns_);
}
}
