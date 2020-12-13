#ifndef SL_TYPED_STRING_HPP_
#define SL_TYPED_STRING_HPP_

#include <ostream>
#include <sdsl/bit_vectors.hpp>

constexpr uint8_t S_TYPE {1};
constexpr uint8_t L_TYPE {0};

template
<
  typename string          = sdsl::int_vector<8>,
  typename sl_type_vector  = sdsl::bit_vector
>
class sl_typed_string
{
public:

  using size_type = typename string::size_type;

  sl_typed_string (string &&input_string_)
  {
    string_ = std::forward<string>(input_string_);
    if (!string_.empty())
    {
      sl_type_vector_.resize(string_.size());
      std::fill(sl_type_vector_.begin(), sl_type_vector_.end(), S_TYPE);
      if (string_.size() >= 2)
      {
        for (size_type position {string_.size() - 1}; position > 0; --position)
        {
          if
          (
            (string_[position - 1] > string_[position])
            ||
            (
              (string_[position - 1] == string_[position])
              &&
              (sl_type_vector_[position] == L_TYPE)
            )
          )
          {
            sl_type_vector_[position - 1] = L_TYPE;
          }
        }
      }
    }
  }

  bool is_s_type (size_type position) const
  {
    return
    (
      (position < string_.size())
      &&
      (sl_type_vector_[position] == S_TYPE)
    );
  }

  bool is_l_type (size_type position) const
  {
    return
    (
      (position < sl_type_vector_.size())
      &&
      (sl_type_vector_[position] == L_TYPE)
    );
  }

  bool is_leftmost_s_type (size_type position) const
  {
    return
    (
      (position < sl_type_vector_.size())
      &&
      (
        (
          (position > 0)
          &&
          (sl_type_vector_[position] == S_TYPE)
          &&
          (sl_type_vector_[position - 1] == L_TYPE)
        )
        ||
        (
          (position == 0)
          &&
          (sl_type_vector_[position] == S_TYPE)
        )
      )
    );
  }

  bool is_rightmost_l_type (size_type position) const
  {
    return
    (
      (position < (sl_type_vector_.size() - 1))
      &&
      (
        (sl_type_vector_[position] == L_TYPE)
        &&
        (sl_type_vector_[position + 1] == S_TYPE)
      )
    );
  }

  auto operator[] (size_type position)
  {
    return string_[position];
  }

  size_type size () const
  {
    return string_.size();
  }

  friend std::ostream& operator<< (std::ostream &out, sl_typed_string const &sl_typed_string_)
  {
    for (auto character : sl_typed_string_.string_)
    {
      out << std::setw(4) << character;
    }
    out << '\n';
    for (size_type position {0}; position < sl_typed_string_.size(); ++position)
    {
      if (sl_typed_string_.is_leftmost_s_type(position))
      {
        out << std::setw(4) << "S*";
      }
      else if (sl_typed_string_.is_rightmost_l_type(position))
      {
        out << std::setw(4) << "L*";
      }
      else if (sl_typed_string_.is_s_type(position))
      {
        out << std::setw(4) << "S";
      }
      else if (sl_typed_string_.is_l_type(position))
      {
        out << std::setw(4) << "L";
      }
      else
      {
        out << std::setw(4) << "N/A";
      }
    }
    return out;
  }

private:

  string         string_;
  sl_type_vector sl_type_vector_;

};

#endif
