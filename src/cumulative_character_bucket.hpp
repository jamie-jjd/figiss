#ifndef CUMULATIVE_CHARACTER_BUCKET_
#define CUMULATIVE_CHARACTER_BUCKET_

#include <algorithm>
#include <iostream>
#include <numeric>
#include <sdsl/int_vector.hpp>

template
<
  typename Text_Type = sdsl::int_vector<8>
>
class Cumulative_Character_Bucket
{
public:

  using Char_Type = typename Text_Type::value_type;

  Cumulative_Character_Bucket (Text_Type const &text)
  {
    auto character_upper_bound {*std::max_element(text.begin(), text.end()) + 1};
    cumulative_character_bucket.resize(character_upper_bound);
  }

  void calculate_cumulative_bucket_begin (Text_Type const &text)
  {
    calculate_cumulative_bucket_end(text);
    if (cumulative_character_bucket.size() > 1)
    {
      std::copy
      (
        cumulative_character_bucket.begin(),
        std::prev(cumulative_character_bucket.end()),
        std::next(cumulative_character_bucket.begin())
      );
    }
    cumulative_character_bucket[0] = 0;
    return;
  }

  void calculate_cumulative_bucket_end (Text_Type const &text)
  {
    calculate_character_bucket(text);
    std::partial_sum
    (
      cumulative_character_bucket.begin(),
      cumulative_character_bucket.end(),
      cumulative_character_bucket.begin()
    );
    return;
  }

  uint32_t get_current_cumulative_bucket_begin (Char_Type character)
  {
    return (cumulative_character_bucket[character]++);
  }

  uint32_t get_current_cumulative_bucket_end (Char_Type character)
  {
    return (--cumulative_character_bucket[character]);
  }

  friend std::ostream& operator<< (std::ostream &out, Cumulative_Character_Bucket const &ccb)
  {
    for (auto const &value : ccb.cumulative_character_bucket)
    {
      out << value << ' ';
    }
    return out;
  }

private:

  sdsl::int_vector<32> cumulative_character_bucket;

  void calculate_character_bucket (Text_Type const &text)
  {
    sdsl::util::_set_zero_bits(cumulative_character_bucket);
    for (auto const &character : text)
    {
      ++cumulative_character_bucket[character];
    }
    return;
  }
};

#endif
