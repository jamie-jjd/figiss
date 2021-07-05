#pragma once

#include <filesystem>
#include <fstream>
#include <functional>
#include <random>

#include <sdsl/int_vector.hpp>

class Patterns
{
public:

  Patterns
  (
    std::filesystem::path const &text_path,
    uint64_t const amount,
    uint64_t const unit_size,
    bool const is_mutated = true
  );

  void Serialize (std::filesystem::path const &path);
  void Load (std::filesystem::path const &path);

private:

  uint64_t amount_;
  uint64_t unit_size_;
  sdsl::int_vector<8> labels_;

};

Patterns::Patterns
(
  std::filesystem::path const &text_path,
  uint64_t const amount,
  uint64_t const unit_size,
  bool const is_mutated
)
{
  {
    amount_ = amount;
    unit_size_ = unit_size;
  }
  sdsl::int_vector<8> text;
  sdsl::load_vector_from_file(text, text_path);
  if (unit_size_ > std::size(text))
  {
    std::cout << "warning: pattern size is replaced with text size\n";
    unit_size_ = std::size(text);
  }
  {
    std::cout
    << "generate " << amount_
    << " random patterns of size " << unit_size_
    << " from " << std::filesystem::canonical(text_path) << "\n";
  }
  std::mt19937 engine {std::random_device{}()};
  std::uniform_int_distribution<uint64_t> dist_text {0, std::size(text) - unit_size_};
  std::uniform_int_distribution<uint64_t> dist_pattern {0, unit_size_ - 1};
  labels_.resize(amount_ * unit_size_);
  auto it {std::begin(labels_)};
  for (uint64_t i {}; i != amount_; ++i)
  {
    auto text_first {std::next(std::begin(text), dist_text(engine))};
    std::copy(text_first, std::next(text_first, unit_size_), it);
    if (is_mutated)
    {
      std::swap
      (
        *std::next(it, dist_pattern(engine)),
        *std::next(it, dist_pattern(engine))
      );
    }
    it = std::next(it, unit_size_);
  }
  return;
}

void Patterns::Serialize (std::filesystem::path const &path)
{
  if (std::size(labels_) == (amount_ * unit_size_))
  {
    std::cout << "serialize patterns to " << std::filesystem::canonical(path) << "\n";
    std::ofstream fout {path, std::ios_base::out | std::ios_base::trunc};
    sdsl::write_member(amount_, fout);
    sdsl::write_member(unit_size_, fout);
    sdsl::serialize(labels_, fout);
  }
  return;
}

void Patterns::Load (std::filesystem::path const &path)
{
  std::ifstream fin {path};
  std::cout << "load patterns from " << std::filesystem::canonical(path) << "\n";
  sdsl::read_member(amount_, fin);
  sdsl::read_member(unit_size_, fin);
  labels_.load(fin);
  return;
}
