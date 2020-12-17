#include <algorithm>
#include <memory>
#include <numeric>
#include <stdexcept>

#include <sdsl/suffix_trees.hpp>

#include "calculate_grammar.hpp"
#include "grammar_rule_trie.hpp"
#include "sl_type.hpp"

template
<
  typename grammar_compressed_text_type,
  typename grammar_compressed_text_index_type,
  typename temporary_grammar_compressed_text_iterator_type
>
void calculate_grammar_compressed_text_index
(
  grammar_compressed_text_type &grammar_compressed_text,
  grammar_compressed_text_index_type &grammar_compressed_text_index,
  temporary_grammar_compressed_text_iterator_type temporary_grammar_compressed_text_begin,
  temporary_grammar_compressed_text_iterator_type temporary_grammar_compressed_text_end
)
{
  grammar_compressed_text.resize
  (
    std::distance
    (
      temporary_grammar_compressed_text_begin,
      temporary_grammar_compressed_text_end
    )
  );
  std::copy
  (
    temporary_grammar_compressed_text_begin,
    temporary_grammar_compressed_text_end,
    grammar_compressed_text.begin()
  );
  sdsl::construct_im
  (
    grammar_compressed_text_index,
    grammar_compressed_text
  );
  return;
}

int main (int argc, char **argv)
{
  if (argc != 2)
  {
    throw std::runtime_error(std::string("usage: ") + argv[0] + " [input path]");
  }

  sdsl::int_vector<8> text;
  sdsl::bit_vector sl_type_vector;
  sdsl::int_vector<> text_distance_vector;
  sdsl::int_vector<> grammar_compressed_text;
  sdsl::csa_wt<sdsl::wt_int<>> grammar_compressed_text_index;

  sdsl::load_vector_from_file(text, argv[1], 1);
  sdsl::append_zero_symbol(text);

  auto
  [
    grammar_rule_begin_distance_vector_begin,
    grammar_rule_begin_distance_vector_end,
    temporary_grammar_compressed_text_begin,
    temporary_grammar_compressed_text_end
  ]
  {
    calculate_grammar
    (
      text,
      sl_type_vector,
      text_distance_vector
    )
  };

  calculate_grammar_compressed_text_index
  (
    grammar_compressed_text,
    grammar_compressed_text_index,
    temporary_grammar_compressed_text_begin,
    temporary_grammar_compressed_text_end
  );

  grammar_rule_trie<> grammar_rule_trie
  (
    text,
    sl_type_vector,
    grammar_rule_begin_distance_vector_begin,
    grammar_rule_begin_distance_vector_end
  );

  return 0;
}
