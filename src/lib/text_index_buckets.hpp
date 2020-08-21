#ifndef TEXT_INDEX_BUCKETS_
#define TEXT_INDEX_BUCKETS_

#include <algorithm>
#include <cstdint>
#include <numeric>
#include <vector>

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
    for (auto i {static_cast<int32_t>(0)}; i != alphabet_size_; ++i) { bucket_begins_[i] -= character_counts_[i]; }
  }

  void generate_bucket_ends ()
  {
    std::partial_sum(character_counts_.begin(), character_counts_.end(), bucket_ends_.begin());
  }

  int32_t get_begin (Character_Type character)
  {
    return (bucket_begins_[character]++);
  }

  int32_t get_end (Character_Type character)
  {
    return (--bucket_ends_[character]);
  }
};

#endif
