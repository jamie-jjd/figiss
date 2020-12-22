#ifndef SL_TYPE_HPP_
#define SL_TYPE_HPP_

constexpr uint8_t S {1};
constexpr uint8_t L {0};

template <typename sl_types_iterator_type>
constexpr bool is_rightmost_l (sl_types_iterator_type it)
{
  return
  (
    (*it == L)
    &&
    (*std::next(it) == S)
  );
}

template <typename sl_types_iterator_type>
constexpr bool is_leftmost_s (sl_types_iterator_type it)
{
  return
  (
    (*it == S)
    &&
    (*std::prev(it) == L)
  );
}

template
<
  typename string_type,
  typename sl_types_type
>
void calculate_sl_types
(
  string_type const &string,
  sl_types_type &sl_types
)
{
  sl_types.resize(string.size());
  std::fill
  (
    std::begin(sl_types),
    std::end(sl_types),
    S
  );
  auto string_rit {std::prev(std::end(string))};
  auto string_rend {std::begin(string)};
  auto rit {std::prev(std::end(sl_types))};
  while (string_rit != string_rend)
  {
    if
    (
      (*std::prev(string_rit) > *string_rit)
      ||
      (
        (*std::prev(string_rit) == *string_rit)
        &&
        (*rit == L)
      )
    )
    {
      *std::prev(rit) = L;
    }
    --rit;
    --string_rit;
  }
  return;
}

#endif
