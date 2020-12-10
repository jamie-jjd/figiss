#ifndef SL_TYPE_VECTOR_HPP_
#define SL_TYPE_VECTOR_HPP_

#include <ostream>
#include <sdsl/bit_vectors.hpp>

template
<
  typename Bit_Vector_Type = sdsl::bit_vector,
  typename Int_Vector_Type = sdsl::int_vector<8>
>
class SL_Type_Vector
{
public:
  using size_type = typename Bit_Vector_Type::size_type;

  static constexpr uint8_t S_type {1};
  static constexpr uint8_t L_type {0};
  Bit_Vector_Type SL_type_vector;

  SL_Type_Vector (Int_Vector_Type const &text)
  {
    if (!text.empty())
    {
      SL_type_vector.resize(text.size());
      if (text.size() > 1)
      {
        std::fill(SL_type_vector.begin(), SL_type_vector.end(), SL_Type_Vector::S_type);
        for (auto position {text.size() - 1}; position != 0; --position)
        {
          if
          (
            (text[position - 1] > text[position])
            ||
            (
              (text[position - 1] == text[position])
              &&
              (SL_type_vector[position] == SL_Type_Vector::L_type)
            )
          )
          {
            SL_type_vector[position - 1] = SL_Type_Vector::L_type;
          }
        }
      }
    }
  }

  bool is_S_type (size_type position) const
  {
    return (SL_type_vector[position] == SL_Type_Vector::S_type);
  }

  bool is_L_type (size_type position) const
  {
    return (SL_type_vector[position] == SL_Type_Vector::L_type);
  }

  bool is_leftmost_S_type (size_type position) const
  {
    return
    (
      (position == (SL_type_vector.size() - 1))
      ||
      (
        (position == 0)
        &&
        (SL_type_vector[position] == SL_Type_Vector::S_type)
      )
      ||
      (
        (SL_type_vector[position] == SL_Type_Vector::S_type)
        &&
        (SL_type_vector[position - 1] == SL_Type_Vector::L_type)
      )
    );
  }

  bool is_rightmost_L_type (size_type position) const
  {
    return
    (
      (position != (SL_type_vector.size() - 1))
      &&
      (
        (SL_type_vector[position] == SL_Type_Vector::L_type)
        &&
        (SL_type_vector[position + 1] == SL_Type_Vector::S_type)
      )
    );
  }
};

#endif
