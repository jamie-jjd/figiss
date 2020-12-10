#ifndef GRAMMAR_COMPRESSED_INDEX_HPP_
#define GRAMMAR_COMPRESSED_INDEX_HPP_

#include <sdsl/wavelet_trees.hpp>

template
<
  typename Text_Type,
  typename Pattern_Type,
  typename Grammar_Rule_Set_Type,
  typename Grammar_Compressed_Text_Type
>
class Grammar_Compressed_Index
{
public:

  using Size_Type = typename Text_Type::size_type;

  Grammar_Compressed_Index (char const *file)
  {
    Text_Type text;
    sdsl::load_vector_from_file(text, file, 1);
    sdsl::append_zero_symbol(text);
    parse_text(text, grammar_rule_set, grammar_compressed_text);
  }

  Size_Type count (Pattern_Type const &pattern)
  {
    return 0;
  }

private:

  Grammar_Rule_Set_Type         grammar_rule_set;
  Grammar_Compressed_Text_Type  grammar_compressed_text;

};

#endif
