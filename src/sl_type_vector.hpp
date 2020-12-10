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

  using Size_Type = typename Bit_Vector_Type::size_type;

  static constexpr uint8_t S_type {1};
  static constexpr uint8_t L_type {0};

  SL_Type_Vector (Int_Vector_Type const &text)
  {
    if (!text.empty())
    {
      SL_type_vector.resize(text.size());
      std::fill(SL_type_vector.begin(), SL_type_vector.end(), SL_Type_Vector::S_type);
      if (text.size() > 1)
      {
        for (Size_Type position {text.size() - 1}; position != 0; --position)
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

  bool is_S_type (Size_Type position) const
  {
    return
    (
      (position < SL_type_vector.size())
      &&
      (SL_type_vector[position] == SL_Type_Vector::S_type)
    );
  }

  bool is_L_type (Size_Type position) const
  {
    return
    (
      (position < SL_type_vector.size())
      &&
      (SL_type_vector[position] == SL_Type_Vector::L_type)
    );
  }

  bool is_leftmost_S_type (Size_Type position) const
  {
    return
    (
      (position < SL_type_vector.size())
      &&
      (
        (
          (position > 0)
          &&
          (SL_type_vector[position] == SL_Type_Vector::S_type)
          &&
          (SL_type_vector[position - 1] == SL_Type_Vector::L_type)
        )
        ||
        (
          (position == 0)
          &&
          (SL_type_vector[position] == SL_Type_Vector::S_type)
        )
      )
    );
  }

  bool is_rightmost_L_type (Size_Type position) const
  {
    return
    (
      (position < (SL_type_vector.size() - 1))
      &&
      (
        (SL_type_vector[position] == SL_Type_Vector::L_type)
        &&
        (SL_type_vector[position + 1] == SL_Type_Vector::S_type)
      )
    );
  }

  friend std::ostream& operator<< (std::ostream &out, SL_Type_Vector const &sltv)
  {
    for (Size_Type position {0}; position != sltv.SL_type_vector.size(); ++position)
    {
      out << std::setw(4) << position;
    }
    std::cout << '\n';
    for (Size_Type position {0}; position != sltv.SL_type_vector.size(); ++position)
    {
      if (sltv.is_leftmost_S_type(position))
      {
        out << std::setw(4) << "S*";
      }
      else if (sltv.is_rightmost_L_type(position))
      {
        out << std::setw(4) << "L*";
      }
      else if (sltv.is_S_type(position))
      {
        out << std::setw(4) << "S";
      }
      else if (sltv.is_L_type(position))
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

  Bit_Vector_Type SL_type_vector;

};

#endif
