#ifndef GRAMMAR_COMPRESSED_INDEX_HPP_
#define GRAMMAR_COMPRESSED_INDEX_HPP_

#include <vector>
#include <sdsl/int_vector.hpp>
#include <sdsl/wavelet_trees.hpp>


void construct_non_terminal_text_index
(
  csa_wt<wt_int<>> &non_terminal_text_index,
  sdsl::int_vector<32> &non_terminal_text
)
{
  sdsl::constrcut_im(non_terminal_text_index, non_terminal_text, 4);
  return;
}

void construct_grammar_rules
(
  csa_wt<> &grammar_rules,
  std::vector<sdsl::int_vector<8>> const &terminal_substring_vector
)
{
  // TODO: constrcut grammar_rules by terminal_substring_vector
  return;
}

class Grammar_Compressed_Index
{
public:

  Grammar_Compressed_Index (char const *file)
  {
    sdsl::int_vector<8> text;
    sdsl::load_vector_from_file(text, file, 1);

    std::vector<sdsl::int_vector<8>> terminal_substring_vector;
    sdsl::int_vector<32> nonterminal_vector;
    parse_text(text, terminal_substring_vector, non_terminal_text);

    construct_grammar_rules(grammar_rules, terminal_substring_vector);
    construct_non_terminal_text_index(non_terminal_text_index, non_terminal_text)
  }

  size_t count (sdsl::int_vector<8> const &pattern)
  {
    // TODO: details of couting
    return 0;
  }

  csa_wt<> grammar_rules;
  csa_wt<wt_int<>> non_terminal_text_index;
};

#endif
