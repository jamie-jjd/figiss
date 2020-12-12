#ifndef SL_TYPED_STRING_HPP_
#define SL_TYPED_STRING_HPP_

#include <ostream>
#include <sdsl/bit_vectors.hpp>

constexpr uint8_t TYPE_S {1};
constexpr uint8_t TYPE_L {0};

template
<
  typename String_Type          = sdsl::int_vector<8>,
  typename SL_Type_Vector_Type  = sdsl::bit_vector
>
class SL_Typed_String
{
public:

  using Size_Type = typename String_Type::size_type;

  SL_Typed_String (String_Type &input_string)
  {
    string = std::move(input_string);
    if (!string.empty())
    {
      SL_type_vector.resize(string.size());
      std::fill(SL_type_vector.begin(), SL_type_vector.end(), TYPE_S);
      if (string.size() >= 2)
      {
        for (Size_Type position {string.size() - 1}; position > 0; --position)
        {
          if
          (
            (string[position - 1] > string[position])
            ||
            (
              (string[position - 1] == string[position])
              &&
              (SL_type_vector[position] == TYPE_L)
            )
          )
          {
            SL_type_vector[position - 1] = TYPE_L;
          }
        }
      }
    }
  }

  bool is_S_type (Size_Type position) const
  {
    return
    (
      (position < string.size())
      &&
      (SL_type_vector[position] == TYPE_S)
    );
  }

  bool is_L_type (Size_Type position) const
  {
    return
    (
      (position < SL_type_vector.size())
      &&
      (SL_type_vector[position] == TYPE_L)
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
          (SL_type_vector[position] == TYPE_S)
          &&
          (SL_type_vector[position - 1] == TYPE_L)
        )
        ||
        (
          (position == 0)
          &&
          (SL_type_vector[position] == TYPE_S)
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
        (SL_type_vector[position] == TYPE_L)
        &&
        (SL_type_vector[position + 1] == TYPE_S)
      )
    );
  }

  auto operator[] (Size_Type position)
  {
    return string[position];
  }

  Size_Type size () const
  {
    return string.size();
  }

  friend std::ostream& operator<< (std::ostream &out, SL_Typed_String const &SL_typed_string)
  {
    for (auto character : SL_typed_string.string)
    {
      out << std::setw(4) << character;
    }
    out << '\n';
    for (Size_Type position {0}; position < SL_typed_string.size(); ++position)
    {
      if (SL_typed_string.is_leftmost_S_type(position))
      {
        out << std::setw(4) << "S*";
      }
      else if (SL_typed_string.is_rightmost_L_type(position))
      {
        out << std::setw(4) << "L*";
      }
      else if (SL_typed_string.is_S_type(position))
      {
        out << std::setw(4) << "S";
      }
      else if (SL_typed_string.is_L_type(position))
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

  String_Type         string;
  SL_Type_Vector_Type SL_type_vector;

};

#endif
