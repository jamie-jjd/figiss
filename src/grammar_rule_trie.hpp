#ifndef GRAMMAR_RULE_TRIE_HPP_
#define GRAMMAR_RULE_TRIE_HPP_

#include <map>
#include <memory>
#include <sdsl/int_vector.hpp>

template <grammar_rule_vector_type>
struct grammar_rule_trie_node
{
};

template
<
  grammar_rule_vector_type = sdsl::int_vector<>
>
struct grammar_rule_trie
{
  using difference_type = typename grammar_rule_vector_type::difference_type;
  using character_type = typename grammar_rule_vector_type::value_type;
  using size_type = typename grammar_rule_vector_type::size_type;

  struct trie_node
  {
    using edge_range_type = std::pair<difference_type, difference_type>;
    using non_terminal_range_type = std::pair<size_type, size_type>;
    using trie_node_pointer_type = std::shared_ptr<trie_node>;
    using branch_type = std::map<character_type, trie_node_pointer_type>;

    edge_range_type edge_range;
    branch_type branch;
    non_terminal_range_type non_terminal_range;

    grammar_rule_trie_node
    (
      edge_range_type const &edge_range_
    ) : edge_range {edge_range_}
    {
    }

    grammar_rule_trie_node
    (
      edge_range_type const &edge_range_,
      size_type const &non_terminal
    )
    : edge_range {edge_range_},
      non_terminal_range {{non_terminal, non_terminal}}
    {
    }
  };

  using trie_node_pointer_type = std::shared_ptr<trie_node>;

  grammar_rule_vector_type grammar_rule_vector;
  trie_node_pointer_type root;
  trie_node_pointer_type reverse_root;

  grammar_rule_trie (grammar_rule_vector_type &&grammar_rule_vector_)
  : grammar_rule_vector {std::move(grammar_rule_vector_)}
  {
    insert_grammar_rule_set();
    insert_reverse_grammar_rule_set();
  }

  void insert_grammar_rule_vector ()
  {
    size_type non_terminal {1};
    for (auto it {std::next(grammar_rule_vector.begin())}; it != grammar_rule_vector.end(); std::advance(it))
    {
      auto node {root};
      while (true)
      {
        auto child {node->branch[*it]};
        if (child == node->branch.end())
        {
          node->branch[character] = std::make_shared<grammar_rule_trie_node_type>({distance, end}, non_terminal);
          return;
        }
        auto distance_edge {std::get<0>(child->edge_range)};
        auto end_edge {std::get<1>(child->edge_range)};
        while
        (
          (distance_edge < end_edge)
          &&
          (distance < end)
          &&
          (grammar_rule_vector[distance_edge] == grammar_rule_vector[distance])
        )
        {
          ++distance_child;
          ++distance;
        }
        if (distance_edge < end_edge)
        {
        }
      }
    }
    return;
  }

  void insert_reverse_grammar_rule_vector ()
  {

  }
};

#endif
