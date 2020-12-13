#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <sdsl/suffix_trees.hpp>

constexpr uint8_t S_TYPE {1};
constexpr uint8_t L_TYPE {0};

template <typename sl_type_vector>
constexpr bool is_rightmost_l_type
(
  sl_type_vector const &sl_type_vector_,
  uint64_t const &position
)
{
  return
  (
    (sl_type_vector_[position] == L_TYPE)
    &&
    (sl_type_vector_[position + 1] == S_TYPE)
  );
}

template <typename sl_type_vector>
constexpr bool is_leftmost_s_type
(
  sl_type_vector const &sl_type_vector_,
  uint64_t const &position
)
{
  return
  (
    (sl_type_vector_[position] == S_TYPE)
    &&
    (sl_type_vector_[position - 1] == L_TYPE)
  );
}

template
<
  typename text,
  typename character_bucket_range_vector
>
void calculate_cumulative_character_bucket_end
(
  text const &text_,
  character_bucket_range_vector &character_bucket_range_vector_
)
{
  std::fill
  (
    character_bucket_range_vector_.begin(),
    character_bucket_range_vector_.end(),
    0
  );
  for (auto const &character : text_)
  {
    ++character_bucket_range_vector_[character];
  }
  std::partial_sum
  (
    character_bucket_range_vector_.begin(),
    character_bucket_range_vector_.end(),
    character_bucket_range_vector_.begin()
  );
  return;
}

template
<
  typename text,
  typename character_bucket_range_vector
>
void calculate_cumulative_character_bucket_begin
(
  text const &text_,
  character_bucket_range_vector &character_bucket_range_vector_
)
{
  calculate_cumulative_character_bucket_end(text_, character_bucket_range_vector_);
  for (auto character {character_bucket_range_vector_.size() - 1}; character > 0; --character)
  {
    character_bucket_range_vector_[character] = character_bucket_range_vector_[character - 1];
  }
  character_bucket_range_vector_[0] = 0;
  return;
}

int main (int argc, char **argv)
{
  if (argc != 2)
  {
    throw std::runtime_error(std::string("usage: ") + argv[0] + " [input path]");
  }

  // read text from file
  sdsl::int_vector<8> text;
  sdsl::load_vector_from_file(text, argv[1], 1);
  sdsl::append_zero_symbol(text);

  /* initialize sl_type_vector of text */
  sdsl::bit_vector sl_type_vector(text.size(), S_TYPE);
  for (auto position {text.size() - 1}; position > 0; --position)
  {
    if
    (
      (text[position - 1] > text[position])
      ||
      (
        (text[position - 1] == text[position])
        &&
        (sl_type_vector[position] == L_TYPE)
      )
    )
    {
      sl_type_vector[position - 1] = L_TYPE;
    }
  }

  // entry of the vector is invalid if value is text.size()
  sdsl::int_vector<> text_position_vector(text.size(), text.size());

  // bucket-sort rightmost L-type characters
  auto character_upper_bound {(1ULL << text.width())};
  sdsl::int_vector<> character_bucket_range_vector(character_upper_bound, 0);
  calculate_cumulative_character_bucket_begin(text, character_bucket_range_vector);
  for (auto position {0ULL}; position < (text.size() - 1); ++position)
  {
    if (is_rightmost_l_type(sl_type_vector, position))
    {
      auto character {text[position]};
      auto current_bucket_begin {character_bucket_range_vector[character]++};
      text_position_vector[current_bucket_begin] = position;
    }
  }

  /* induce-sort L-type characters */
  for (auto rank {1ULL}; rank < text_position_vector.size(); ++rank)
  {
    auto position {text_position_vector[rank]};
    if
    (
      (position != text.size())
      &&
      (position > 0)
      &&
      (sl_type_vector[position - 1] == L_TYPE)
    )
    {
      auto character {text[position - 1]};
      auto current_bucket_begin {character_bucket_range_vector[character]++};
      text_position_vector[current_bucket_begin] = (position - 1);
      text_position_vector[rank] = text.size();
    }
  }

  /* induce-sort S-type characters */
  calculate_cumulative_character_bucket_end(text, character_bucket_range_vector);
  for (auto rank {text_position_vector.size() - 1}; rank > 0; --rank)
  {
    auto position {text_position_vector[rank]};
    if
    (
      (position != text.size())
      &&
      (position > 0)
      &&
      (sl_type_vector[position - 1] == S_TYPE)
    )
    {
      auto character {text[position - 1]};
      auto current_bucket_end {--character_bucket_range_vector[character]};
      text_position_vector[current_bucket_end] = (position - 1);
      text_position_vector[rank] = text.size();
    }
  }

  /*
    move valid text_position_vector entries to front stably
    , and reuse space of invalid entries for temporarily
    storing characters of grammar-compressed text
  */
  auto begin_spare_space
  {
    0ULL
    +
    std::distance
    (
      text_position_vector.begin(),
      std::stable_partition
      (
        text_position_vector.begin(),
        text_position_vector.end(),
        [&] (auto const &position)
        {
          return (position != text.size());
        }
      )
    )
  };

  auto non_terminal {0ULL};
  auto begin_previous_grammar_rule {text.size() - 1};
  for (auto rank {0ULL}; rank < begin_spare_space; ++rank)
  {
    auto begin_current_grammar_rule {static_cast<uint64_t>(text_position_vector[rank])};
    auto temporary_position {begin_spare_space + ((begin_current_grammar_rule + 1) / 2)};
    if
    (
      (text[begin_previous_grammar_rule] == text[begin_current_grammar_rule])
      &&
      (sl_type_vector[begin_previous_grammar_rule] == sl_type_vector[begin_current_grammar_rule])
    )
    {
      auto position_previous_grammar_rule {begin_previous_grammar_rule + 1};
      auto position_current_grammar_rule {begin_current_grammar_rule + 1};
      while // compare grammar rule
      (
        !is_leftmost_s_type(sl_type_vector, position_previous_grammar_rule)
        &&
        !is_leftmost_s_type(sl_type_vector, position_current_grammar_rule)
        &&
        (text[position_previous_grammar_rule] == text[position_current_grammar_rule])
        &&
        (sl_type_vector[position_previous_grammar_rule] == sl_type_vector[position_current_grammar_rule])
      )
      {
        ++position_previous_grammar_rule;
        ++position_current_grammar_rule;
      }
      if // same grammar rule
      (
        is_leftmost_s_type(sl_type_vector, position_previous_grammar_rule)
        &&
        is_leftmost_s_type(sl_type_vector, position_current_grammar_rule)
      )
      {
        --non_terminal;
        text_position_vector[rank] = text.size();
      }
    }
    text_position_vector[temporary_position]= ++non_terminal;
    begin_previous_grammar_rule = begin_current_grammar_rule;
  }

  /* collect characters of grammar-compressed text */
  sdsl::int_vector<> grammar_compressed_text(begin_spare_space);
  std::copy_if
  (
    std::next(text_position_vector.begin(), begin_spare_space),
    text_position_vector.end(),
    grammar_compressed_text.begin(),
    [&] (auto const &position)
    {
      return (position != text.size());
    }
  );

  std::cout << grammar_compressed_text << '\n';
  /* construct FM-index on grammar-compressed text */
  sdsl::csa_wt<sdsl::wt_int<>> grammar_compressed_text_fm_index;
  sdsl::construct_im(grammar_compressed_text_fm_index, grammar_compressed_text);

  return 0;
}
