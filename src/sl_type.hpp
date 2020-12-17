#ifndef SL_TYPE_HPP_
#define SL_TYPE_HPP_

constexpr uint8_t S_TYPE {1};
constexpr uint8_t L_TYPE {0};

template <typename sl_type_vector_iterator_type>
constexpr bool is_rightmost_l_type (sl_type_vector_iterator_type it)
{
  return
  (
    (*it == L_TYPE)
    &&
    (*std::next(it) == S_TYPE)
  );
}

template <typename sl_type_vector_iterator_type>
constexpr bool is_leftmost_s_type (sl_type_vector_iterator_type it)
{
  return
  (
    (*it == S_TYPE)
    &&
    (*std::prev(it) == L_TYPE)
  );
}

template
<
  typename string_type,
  typename sl_type_vector_type
>
void calculate_sl_type_vector
(
  string_type const &string,
  sl_type_vector_type &sl_type_vector
)
{
  sl_type_vector.resize(string.size());
  std::fill
  (
    std::begin(sl_type_vector),
    std::end(sl_type_vector),
    S_TYPE
  );
  auto string_rit {std::prev(std::end(string))};
  auto string_rend {std::begin(string)};
  auto rit {std::prev(std::end(sl_type_vector))};
  while (string_rit != string_rend)
  {
    if
    (
      (*std::prev(string_rit) > *string_rit)
      ||
      (
        (*std::prev(string_rit) == *string_rit)
        &&
        (*rit == L_TYPE)
      )
    )
    {
      *std::prev(rit) = L_TYPE;
    }
    --rit;
    --string_rit;
  }
  return;
}

#endif
