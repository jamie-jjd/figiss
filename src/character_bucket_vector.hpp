#ifndef CHARACTER_BUCKET_VECTOR_
#define CHARACTER_BUCKET_VECTOR_

#include <algorithm>
#include <numeric>
#include <ostream>
#include <sdsl/int_vector.hpp>

template
<
  typename Text_Type = sdsl::int_vector<8>,
  typename Bucket_Vector_Type = sdsl::int_vector<>
>
class Character_Bucket_Vector
{
public:

  using Char_Type = typename Text_Type::value_type;
  using Value_Type = typename Bucket_Vector_Type::value_type;

  Character_Bucket_Vector (Text_Type const &text)
  {
    auto character_upper_bound {*std::max_element(text.begin(), text.end()) + 1};
    character_bucket_vector.resize(character_upper_bound);
  }

  void calculate_cumulative_bucket_begin (Text_Type const &text)
  {
    calculate_cumulative_bucket_end(text);
    if (character_bucket_vector.size() > 1)
    {
      for (auto i {character_bucket_vector.size() - 1}; i > 0; --i)
      {
        character_bucket_vector[i] = character_bucket_vector[i - 1];
      }
    }
    character_bucket_vector[0] = 0;
    return;
  }

  void calculate_cumulative_bucket_end (Text_Type const &text)
  {
    calculate_character_bucket_sizes(text);
    std::partial_sum
    (
      character_bucket_vector.begin(),
      character_bucket_vector.end(),
      character_bucket_vector.begin()
    );
    return;
  }

  Value_Type get_current_bucket_begin (Char_Type character)
  {
    return (character_bucket_vector[character]++);
  }

  Value_Type get_current_bucket_end (Char_Type character)
  {
    return (--character_bucket_vector[character]);
  }

  friend std::ostream& operator<< (std::ostream &out, Character_Bucket_Vector const &cbv)
  {
    for (Value_Type value : cbv.character_bucket_vector)
    {
      out << std::setw(4) << value;
    }
    return out;
  }

private:

  Bucket_Vector_Type character_bucket_vector;

  void calculate_character_bucket_sizes (Text_Type const &text)
  {
    sdsl::util::_set_zero_bits(character_bucket_vector);
    for (auto const &character : text)
    {
      ++character_bucket_vector[character];
    }
    return;
  }
};

#endif
