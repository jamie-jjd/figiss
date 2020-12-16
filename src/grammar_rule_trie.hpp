#ifndef GRAMMAR_RULE_TRIE_HPP_
#define GRAMMAR_RULE_TRIE_HPP_

#include <map>
#include <memory>
#include <sdsl/int_vector.hpp>

template
<
  grammar_rule_size_vector_type = sdsl::int_vector<>,
  grammar_rule_vector_type = sdsl::int_vector<>
>
struct grammar_rule_trie
{
  using difference_type = typename grammar_rule_vector_type::difference_type;
  using character_type = typename grammar_rule_vector_type::value_type;

  struct trie_node
  {
    using range_type = std::pair<difference_type, difference_type>;
    using trie_node_pointer_type = std::shared_ptr<trie_node>;
    using branch_type = std::map<character_type, trie_node_pointer_type>;

    range_type edge_range;
    branch_type branch;
    range_type non_terminal_range;

    grammar_rule_trie_node
    (
      range_type edge_range_
    )
    : edge_range {edge_range_},
    {
    }

    grammar_rule_trie_node
    (
      range_type edge_range_,
      difference_type non_terminal
    )
    : edge_range {edge_range_},
      non_terminal_range {{non_terminal, non_terminal}}
    {
    }
  };

  using trie_node_pointer_type = std::shared_ptr<trie_node>;
  using non_terminal_vector_type = sdsl::int_vector<>;

  grammar_rule_size_vector_type grammar_rule_size_vector;
  grammar_rule_vector_type grammar_rule_vector;
  trie_node_pointer_type grammar_rule_trie_root;
  trie_node_pointer_type reverse_grammar_rule_trie_root;
  non_terminal_permutation_vector_type non_terminal_permutation_vector;

  grammar_rule_trie
  (
    grammar_rule_size_vector_type &&grammar_rule_size_vector_,
    grammar_rule_vector_type &&grammar_rule_vector_
  )
  : grammar_rule_size_vector {std::move(grammar_rule_size_vector_)},
    grammar_rule_vector {std::move(grammar_rule_vector_)}
  {
    insert_grammar_rule_vector();
    calculate_contiguous_non_terminal_range();
    insert_reverse_grammar_rule_set();
    calculate_non_terminal_permutation_vector();
  }

  void insert_grammar_rule_vector ()
  {
    difference_type non_terminal {1};
    auto distance_it {grammar_rule_size_vector.begin()};
    auto grammar_rule_begin {grammar_rule_vector.begin()};
    while (grammar_rule_begin != grammar_rule_vector.end())
    {
      auto grammar_rule_it {grammar_rule_begin};
      auto grammar_rule_end {std::next(grammar_rule_begin, *distance_it)};
      auto node {grammar_rule_trie_root};
      while (true)
      {
        auto child {node->branch[*grammar_rule_it]};
        if (child == node->branch.end())
        {
          node->branch[*grammar_rule_it] = std::make_shared<trie_node>
          (
            {
              std::distance(grammar_rule_vector.begin(), grammar_rule_it);
              std::distance(grammar_rule_vector.begin(), grammar_rule_end);
            },
            non_terminal
          );
          break;
        }
        else
        {
          auto original_grammar_rule_it {grammar_rule_it};
          auto edge_it {std::next(grammar_rule_vector.begin(), std::get<0>(child->edge_range))};
          auto edge_end {std::next(grammar_rule_vector.begin(), std::get<1>(child->edge_range))};
          while
          (
            (*grammar_rule_it == *edge_it)
            &&
            (*edge_it != edge_end)
          )
          {
            std::advence(grammar_rule_it);
            std::advence(edge_it);
          }
          if (edge_it == edge_end)
          {
            node = child;
          }
          else
          {
            auto edge_it_distance {std::distance(grammar_rule_vector.begin(), edge_it)};
            auto internal_node {std::make_shared<trie_node>({std::get<0>(child->edge_range), edge_it_distance})};
            std::get<0>(child->edge_range) = edge_it_distance;
            node->branch[*original_grammar_rule_it] = internal_node;
            internal_node->branch[*edge_it] = child;
            internal_node->branch[*grammar_rule_it] = std::make_shared<trie_node>
            (
              {
                std::distance(grammar_rule_vector.begin(), grammar_rule_it);
                std::distance(grammar_rule_vector.begin(), grammar_rule_end);
              },
              non_terminal
            );
            break;
          }
        }
      }
      ++non_terminal;
      std::advance(distance_it);
      grammar_rule_begin = grammar_rule_end;
    }
    return;
  }

  void calculate_contiguous_non_terminal_range ()
  {
    calculate_non_terminal_range(grammar_trie_root);
  }

  auto calculate_non_terminal_range (trie_node_pointer_type node)
  {
    if (std::get<0>(node->non_terminal_range) != 0)
    {
      std::get<0>(node->non_terminal_range) = std::get<0>(get_non_terminal_range(node->branch.begin()));
      std::get<1>(node->non_terminal_range) = std::get<1>(get_non_terminal_range(std::prev(node->branch.end())));
      for (auto child {std::next(node->branch.begin())}; child != std::prev(node->branch.end()); ++child)
      {
        get_non_terminal_range(child);
      }
    }
    return node->edge_range;
  }

  void insert_reverse_grammar_rule_vector ()
  {
    difference_type non_terminal {1};
    auto distance_it {grammar_rule_size_vector.begin()};
    auto grammar_rule_begin {grammar_rule_vector.begin()};
    while (grammar_rule_begin != grammar_rule_vector.end())
    {
      auto grammar_rule_it {grammar_rule_begin};
      auto grammar_rule_end {std::next(grammar_rule_begin, *distance_it)};
      auto node {grammar_rule_trie_root};
      while (true)
      {
        auto child {node->branch[*grammar_rule_it]};
        if (child == node->branch.end())
        {
          node->branch[*grammar_rule_it] = std::make_shared<trie_node>
          (
            {
              std::distance(grammar_rule_vector.begin(), grammar_rule_it);
              std::distance(grammar_rule_vector.begin(), grammar_rule_end);
            },
            non_terminal
          );
          break;
        }
        else
        {
          auto original_grammar_rule_it {grammar_rule_it};
          auto edge_it {std::next(grammar_rule_vector.begin(), std::get<0>(child->edge_range))};
          auto edge_end {std::next(grammar_rule_vector.begin(), std::get<1>(child->edge_range))};
          while
          (
            (*grammar_rule_it == *edge_it)
            &&
            (*edge_it != edge_end)
          )
          {
            std::advence(grammar_rule_it);
            std::advence(edge_it);
          }
          if (edge_it == edge_end)
          {
            node = child;
          }
          else
          {
            auto edge_it_distance {std::distance(grammar_rule_vector.begin(), edge_it)};
            auto internal_node {std::make_shared<trie_node>({std::get<0>(child->edge_range), edge_it_distance})};
            std::get<0>(child->edge_range) = edge_it_distance;
            node->branch[*original_grammar_rule_it] = internal_node;
            internal_node->branch[*edge_it] = child;
            internal_node->branch[*grammar_rule_it] = std::make_shared<trie_node>
            (
              {
                std::distance(grammar_rule_vector.begin(), grammar_rule_it);
                std::distance(grammar_rule_vector.begin(), grammar_rule_end);
              },
              non_terminal
            );
            break;
          }
        }
      }
      ++non_terminal;
      std::advance(distance_it);
      grammar_rule_begin = grammar_rule_end;
    }
  }

  void calculate_non_terminal_permutation_vector ()
  {
    calculate_non_terminal_permutation_range(reverse_grammar_trie_root);
    return;
  }

  void calculate_non_terminal_permutation_range (trie_node_pointer_type)
  {
    return;
  }

};

#endif
