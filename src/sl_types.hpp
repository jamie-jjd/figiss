#ifndef SL_TYPE_HPP_
#define SL_TYPE_HPP_

#include <vector>

template <typename Character_Type>
struct SL_Types
{
  static constexpr uint8_t S_type_ {1};
  static constexpr uint8_t L_type_ {0};

  std::vector<uint8_t> SL_types_;
  int32_t text_size_;

  SL_Types (std::vector<Character_Type> const &text)
  {
    if (text.size() > std::numeric_limits<int32_t>::max())
    {
      throw std::runtime_error("out of indexing range");
    }

    text_size_ = static_cast<int32_t>(text.size());

    if (text_size_ > 0)
    {
      SL_types_.resize(text_size_, S_type_);
      if (text_size_ > 1)
      {
        for (auto i {text_size_ - 2}; i != -1; --i)
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
    }
  }

  bool is_S_type (int32_t text_index)
  {
    if (text_index < 0 || text_index > (text_size_ - 1))
    {
      throw std::runtime_error("out of indexing range");
    }
    return (SL_types_[text_index] == S_type_);
  }

  bool is_L_type (int32_t text_index)
  {
    if (text_index < 0 || text_index > (text_size_ - 1))
    {
      throw std::runtime_error("out of indexing range");
    }
    return (SL_types_[text_index] == L_type_);
  }

  bool is_leftmost_S_type (int32_t text_index)
  {
    if (text_index < 0 || text_index > (text_size_ - 1))
    {
      throw std::runtime_error("out of indexing range");
    }
    return
    (
      (text_index == (text_size_ - 1)) ||
      (
        (SL_types_[text_index] == S_type_) &&
        (SL_types_[text_index + 1] == L_type_)
      )
    );
  }

  bool is_rightmost_L_type (int32_t text_index)
  {
    if (text_index < 0 || text_index > (text_size_ - 1))
    {
      throw std::runtime_error("out of indexing range");
    }
    return
    ((SL_types_[text_index] == L_type_) && (SL_types_[text_index + 1] == S_type_));
  }
};

#endif
