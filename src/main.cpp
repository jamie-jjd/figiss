#include <algorithm>
#include <memory>
#include <numeric>
#include <stdexcept>

#include <sdsl/suffix_trees.hpp>

#include "grammar_rule_trie.hpp"

constexpr uint8_t S_TYPE {1};
constexpr uint8_t L_TYPE {0};

int main (int argc, char **argv)
{
  if (argc != 2)
  {
    throw std::runtime_error(std::string("usage: ") + argv[0] + " [input path]");
  }

  sdsl::int_vector<8> text;
  sdsl::int_vector<> compressed_text;
  sdsl::sd_vector<> sorted_lex_compressed_text_range_vector;
  sdsl::wt_int<> compressed_bwt;
  sdsl::wt_int<> colex_compressed_bwt;

  grammar_rule_trie<> grammar_rule_set;

  sdsl::load_vector_from_file(text, argv[1], 1);
  sdsl::append_zero_symbol(text);

  size_t size_compressed_text {0};
  uint8_t prev_sl_type {S_TYPE};
  auto text_rfirst {std::prev(std::end(text), 2)};
  auto text_rit {text_rfirst};
  auto text_rlast {std::prev(std::begin(text))};
  while (text_rit != text_rlast)
  {
    ++size_compressed_text;
    prev_sl_type = L_TYPE;
    --text_rit;
    while
    (
      (text_rit != text_rlast)
      &&
      ! (
          (prev_sl_type == S_TYPE)
          &&
          (*text_rit > *std::next(text_rit))
        )
    )
    {
      if
      (
        (prev_sl_type == L_TYPE)
        &&
        (*text_rit < *std::next(text_rit))
      )
      {
        prev_sl_type = S_TYPE;
      }
      --text_rit;
    }
    grammar_rule_set.insert_grammar_rule
    (
      std::begin(text),
      std::next(text_rit),
      std::next(text_rfirst),
      1
    );
    grammar_rule_set.insert_grammar_rule
    (
      std::begin(text),
      text_rfirst,
      text_rit,
      -1
    );
  }
  grammar_rule_set.print(std::begin(text), grammar_rule_set.root, 1);
  grammar_rule_set.print(std::begin(text), grammar_rule_set.colex_root, -1);
  return 0;
}
