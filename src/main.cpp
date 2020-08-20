#include <algorithm>
#include <fstream>
#include <map>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <vector>

// #include "sl_types.hpp"
// #include "text_index_buckets.hpp"
#include "lib/utility.hpp"

template <typename Character_Type>
struct SL_Types
{
  constexpr static uint8_t S_type_ {1};
  constexpr static uint8_t L_type_ {0};

  std::vector<uint8_t> SL_types_;

  SL_Types (std::vector<Character_Type> const &text)
  {
    int32_t text_size {static_cast<int32_t>(text.size())};

    SL_types_.resize(text_size, S_type_);
    for (int32_t i {text_size - 2}; i != -1; --i)
    {
      if
      (
        (text[i] > text[i + 1]) ||
        (
          (text[i] == text[i + 1]) &&
          (SL_types_[i + 1] == L_type_)
        )
      )
      {
        SL_types_[i] = L_type_;
      }
    }
  }

  bool is_S_type (int32_t text_index)
  {
    return (SL_types_[text_index] == S_type_);
  }

  bool is_L_type (int32_t text_index)
  {
    return (SL_types_[text_index] == L_type_);
  }

  bool is_rightmost_L_type (int32_t text_index)
  {
    return ((SL_types_[text_index] == L_type_) && (SL_types_[text_index + 1] == S_type_));
  }

};

template <typename Character_Type>
struct Text_Index_Buckets
{
  int32_t alphabet_size_;
  std::vector<int32_t> character_counts_;
  std::vector<int32_t> bucket_begins_;
  std::vector<int32_t> bucket_ends_;

  Text_Index_Buckets (std::vector<Character_Type> const &text)
  {
    alphabet_size_ = (static_cast<int32_t>(*std::max_element(text.begin(), text.end())) + 1);

    character_counts_.resize(alphabet_size_, 0);
    for (auto const &character : text) { ++character_counts_[character]; }

    bucket_begins_.resize(alphabet_size_, 0);
    bucket_ends_.resize(alphabet_size_, 0);
  }

  void generate_bucket_begins ()
  {
    std::partial_sum(character_counts_.begin(), character_counts_.end(), bucket_begins_.begin());
    for (int32_t i {0}; i != alphabet_size_; ++i) { bucket_begins_[i] -= character_counts_[i]; }
  }

  void generate_bucket_ends ()
  {
    std::partial_sum(character_counts_.begin(), character_counts_.end(), bucket_ends_.begin());
  }

  int32_t get_end (int32_t character)
  {
    return (--bucket_ends_[character]);
  }

  int32_t get_begin (int32_t character)
  {
    return (bucket_begins_[character]++);
  }

};

template <typename Character_Type>
void text_parsing
(
  std::vector<int32_t> &alphabet_sizes,
  std::vector<int32_t> &grammar_rules_sizes,
  std::vector<Character_Type> const &text,
  std::vector<int32_t> &text_sizes
)
{

  auto SL_types {SL_Types<Character_Type>(text)};
  auto text_index_buckets {Text_Index_Buckets<Character_Type>(text)};
  auto text_indices {std::vector<int32_t>(text_sizes.back(), -1)};

  int32_t text_size {text_sizes.back()};
  auto &text_indices_size {text_size};

  text_index_buckets.generate_bucket_begins();
  for (int32_t i {text_size - 2}; i != -1; --i)
  {
    if (SL_types.is_rightmost_L_type(i))
    {
      text_indices[text_index_buckets.get_begin(text[i])] = i;
    }
  }

  for (int32_t i {1}; i != text_indices_size; ++i)
  {
    int32_t j {text_indices[i] - 1};
    if ((j >= 0) && SL_types.is_L_type(j))
    {
      text_indices[text_index_buckets.get_begin(text[j])] = j;
      text_indices[i] = -1;
    }
  }

  text_index_buckets.generate_bucket_ends();
  for (int32_t i {text_indices_size - 1}; i != 0; --i)
  {
    int32_t j {text_indices[i] - 1};
    if ((j >= 0) && SL_types.is_S_type(j))
    {
      text_indices[text_index_buckets.get_end(text[j])] = j;
      text_indices[i] = -1;
    }
  }

  alphabet_sizes.push_back(1);
  grammar_rules_sizes.push_back(0);
  text_sizes.push_back(1);

  int32_t temporary_new_text_size {(text_size / 2) + 1};
  auto temporary_new_text {std::vector<int32_t>(temporary_new_text_size, 0)};
  for (int32_t i {1}, previous_token_begin {text_size - 1}; i != text_indices_size; ++i)
  {
    int32_t current_token_begin {text_indices[i]};
    if (current_token_begin != -1)
    {
      ++(text_sizes.back());
      int32_t previous_token_index {previous_token_begin};
      int32_t current_token_index {current_token_begin};
      while
      (
        (text[previous_token_index] == text[current_token_index]) &&
        !SL_types.is_rightmost_L_type(previous_token_index) &&
        !SL_types.is_rightmost_L_type(current_token_index)
      )
      {
        ++previous_token_index;
        ++current_token_index;
      }
      if
      (
        !SL_types.is_rightmost_L_type(previous_token_index) ||
        !SL_types.is_rightmost_L_type(current_token_index)
      )
      {
        int32_t current_token_end {current_token_begin + 1};
        while (!SL_types.is_rightmost_L_type(current_token_end - 1)) { ++current_token_end; }
        ++(alphabet_sizes.back());
        grammar_rules_sizes.back() += (current_token_end - current_token_begin);
      }
      temporary_new_text[(current_token_begin + 1) / 2] = (alphabet_sizes.back() - 1);
      previous_token_begin = current_token_begin;
    }
  }

  auto &new_text {temporary_new_text};
  for (int32_t i {0}, j {0}; i != temporary_new_text_size; ++i)
  {
    if (temporary_new_text[i] != 0)
    {
      new_text[j++] = temporary_new_text[i];
    }
  }
  new_text.resize(text_sizes.back());

  if (alphabet_sizes.back() == text_sizes.back()) return;

  text_parsing<int32_t>(alphabet_sizes, grammar_rules_sizes, new_text, text_sizes);

  return;
}

int32_t get_effective_alphabet_size (std::vector<uint8_t> const &text)
{
  auto is_effective_characters {std::vector<int32_t>(*std::max_element(text.begin(), text.end()) + 1, 0)};
  for (const auto &character : text) { is_effective_characters[character] = 1; }
  return std::accumulate(is_effective_characters.begin(), is_effective_characters.end(), 0);
}

int main (int argc, char **argv)
{
  if (argc != 3) throw std::runtime_error("usage: ./main [input file] [output file]");

  auto input_file_stream {std::basic_ifstream<char>(argv[1], std::ios_base::binary)};
  input_file_stream.unsetf(std::ios_base::skipws);

  auto output_file_stream {std::basic_ofstream<char>(argv[2])};

  auto text {std::vector<uint8_t>(std::istreambuf_iterator<char>(input_file_stream), std::istreambuf_iterator<char>())};
  text.push_back(0);

  auto text_sizes {std::vector<int32_t>(1, text.size())};
  auto alphabet_sizes {std::vector<int32_t>(1, get_effective_alphabet_size(text))};
  auto grammar_rules_sizes {std::vector<int32_t>(1, 0)};

  text_parsing<uint8_t>(alphabet_sizes, grammar_rules_sizes, text, text_sizes);

  output_comma_separated_vector<int32_t>(output_file_stream, alphabet_sizes);
  output_comma_separated_vector<int32_t>(output_file_stream, grammar_rules_sizes);
  output_comma_separated_vector<int32_t>(output_file_stream, text_sizes);

  return 0;
}
