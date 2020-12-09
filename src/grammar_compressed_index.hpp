#ifndef GRAMMAR_COMPRESSED_INDEX_HPP_
#define GRAMMAR_COMPRESSED_INDEX_HPP_

#include <vector>
#include <sdsl/int_vector.hpp>
#include <sdsl/wavelet_trees.hpp>

class Grammar_Compressed_Index
{
public:

  Grammar_Compressed_Index (char const *file)
  {
    sdsl::int_vector<8> text;
    sdsl::load_vector_from_file(text, file, 1);

    std::vector<sdsl::int_vector<8>> terminal_sequence_vector;
    sdsl::int_vector<32> non_terminal_vector;
    parse_text(text, terminal_sequence_vector, nontermial_text);

    construct(grammar_rules, terminal_sequence_vector);
    construct(fm_index_nontermial_text, nontermial_text)
  }

  size_t count (sdsl::int_vector<8> const &pattern)
  {
    return 0;
  }

  csa_wt<> grammar_rules;
  csa_wt<> fm_index_nontermial_text;
};

#endif
