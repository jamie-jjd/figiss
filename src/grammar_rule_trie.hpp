#ifndef GRAMMAR_RULE_TRIE_HPP_
#define GRAMMAR_RULE_TRIE_HPP_

#include <map>
#include <memory>
#include <sdsl/int_vector.hpp>

template
<
  typename grammar_rule_size_vector_type = sdsl::int_vector<>,
  typename grammar_rule_vector_type = sdsl::int_vector<8>
>
struct grammar_rule_trie
{
  using difference_type = typename grammar_rule_vector_type::difference_type;
  using character_type = typename grammar_rule_vector_type::value_type;
  using size_type = typename grammar_rule_vector_type::size_type;

  struct trie_node
  {
    using edge_range_type = std::pair<difference_type, difference_type>;
    using non_terminal_related_range_type = std::pair<size_type, size_type>;
    using trie_node_pointer_type = std::shared_ptr<trie_node>;
    using branch_type = std::map<character_type, trie_node_pointer_type>;

    edge_range_type edge_range;
    branch_type branch;
    non_terminal_related_range_type non_terminal_related_range;

    trie_node (edge_range_type edge_range_)
    : edge_range {edge_range_}
    {
    }

    trie_node
    (
      edge_range_type edge_range_,
      size_type non_terminal
    )
    : edge_range {edge_range_},
      non_terminal_related_range {{non_terminal, non_terminal}}
    {
    }
  };

  using trie_node_pointer_type = std::shared_ptr<trie_node>;
  using non_terminal_permutation_vector_type = grammar_rule_size_vector_type;

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
    insert_reverse_grammar_rule_vector();
    calculate_non_terminal_permutation_vector();
  }

  void insert_grammar_rule_vector ()
  {
    size_type non_terminal {1};
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
              std::distance(grammar_rule_vector.begin(), grammar_rule_it),
              std::distance(grammar_rule_vector.begin(), grammar_rule_end)
            },
            non_terminal
          );
          break;
        }
        else
        {
          auto &new_child {node->branch[*grammar_rule_it]};
          auto edge_it {std::next(grammar_rule_vector.begin(), std::get<0>(child->edge_range))};
          auto edge_end {std::next(grammar_rule_vector.begin(), std::get<1>(child->edge_range))};
          while
          (
            (*edge_it != edge_end)
            &&
            (*grammar_rule_it == *edge_it)
          )
          {
            std::advance(grammar_rule_it);
            std::advance(edge_it);
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
            new_child = internal_node;
            internal_node->branch[*edge_it] = child;
            internal_node->branch[*grammar_rule_it] = std::make_shared<trie_node>
            (
              {
                std::distance(grammar_rule_vector.begin(), grammar_rule_it),
                std::distance(grammar_rule_vector.begin(), grammar_rule_end)
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
    calculate_non_terminal_related_range(grammar_rule_trie_root);
    return;
  }

  auto calculate_non_terminal_range (trie_node_pointer_type node)
  {
    if (std::get<0>(node->non_terminal_related_range) == 0)
    {
      std::get<0>(node->non_terminal_related_range) = std::get<0>
      (
        calculate_non_terminal_range(node->branch.begin())
      );
      for
      (
        auto child {std::next(node->branch.begin())};
        child != std::prev(node->branch.begin());
        std::advance(child)
      )
      {
        calculate_non_terminal_range(child);
      }
      std::get<1>(node->non_terminal_related_range) = std::get<1>
      (
        calculate_non_terminal_range(std::prev(node->branch.end()))
      );
    }
    return node->non_terminal_related_range;
  }

  void insert_reverse_grammar_rule_vector ()
  {
    size_type non_terminal {grammar_rule_size_vector.size()};
    auto distance_rit {std::prev(grammar_rule_size_vector.end())};
    auto grammar_rule_rbegin {std::prev(grammar_rule_vector.end())};
    while (grammar_rule_rbegin != std::prev(grammar_rule_vector.begin()))
    {
      auto grammar_rule_rit {grammar_rule_rbegin};
      auto grammar_rule_rend {std::prev(grammar_rule_rbegin, *distance_rit)};
      auto node {reverse_grammar_rule_trie_root};
      while (true)
      {
        auto child {node->branch[*grammar_rule_rit]};
        if (child == node->branch.end())
        {
          node->branch[*grammar_rule_rit] = std::make_shared<trie_node>
          (
            {
              std::distance(grammar_rule_vector.begin(), grammar_rule_rit),
              std::distance(grammar_rule_vector.begin(), grammar_rule_rend)
            },
            non_terminal
          );
          break;
        }
        else
        {
          auto &new_child {node->branch[*grammar_rule_rit]};
          auto edge_rit {std::next(grammar_rule_vector.begin(), std::get<0>(child->edge_range))};
          auto edge_rend {std::next(grammar_rule_vector.begin(), std::get<1>(child->edge_range))};
          while
          (
            (grammar_rule_rit != grammar_rule_rend)
            &&
            (edge_rit != edge_rend)
            &&
            (*grammar_rule_rit == *edge_rit)
          )
          {
            std::advance(grammar_rule_rit, -1);
            std::advance(edge_rit, -1);
          }
          if (edge_rit == edge_rend)
          {
            node = child;
          }
          else
          {
            auto edge_rit_distance {std::distance(grammar_rule_vector.begin(), edge_rit)};
            auto internal_node {std::make_shared<trie_node>({std::get<0>(child->edge_range), edge_rit_distance})};
            std::get<0>(child->edge_range) = edge_rit_distance;
            new_child = internal_node;
            internal_node->branch[*edge_rit] = child;
            if (grammar_rule_rit != grammar_rule_rend)
            {
              internal_node->branch[*grammar_rule_rit] = std::make_shared<trie_node>
              (
                {
                  std::distance(grammar_rule_vector.begin(), grammar_rule_rit),
                  std::distance(grammar_rule_vector.begin(), grammar_rule_rend)
                },
                non_terminal
              );
            }
            else
            {
              internal_node->non_terminal_related_range = {non_terminal, non_terminal};
            }
            break;
          }
        }
      }
      --non_terminal;
      std::advance(distance_rit, -1);
      grammar_rule_rbegin = grammar_rule_rend;
    }
  }

  void calculate_non_terminal_permutation_vector ()
  {
    calculate_non_terminal_permutation_vector_range(reverse_grammar_rule_trie_root);
    return;
  }

  auto calculate_non_terminal_permutation_vector_range (trie_node_pointer_type node)
  {
    static difference_type distance {0};
    if (std::get<0>(node->non_terminal_related_range) != 0)
    {
      non_terminal_permutation_vector[distance] = std::get<0>(node->non_terminal_related_range);
      node->non_terminal_related_range = {distance, distance};
      ++distance;
    }
    else
    {
      std::get<0>(node->non_terminal_related_range) = std::get<0>
      (
        calculate_non_terminal_permutation_vector_range(node->branch.begin())
      );
      for
      (
        auto child {std::next(node->branch_begin)};
        child != std::prev(node->branch_end);
        std::advance(child)
      )
      {
        calculate_non_terminal_range(child);
      }
      std::get<1>(node->non_terminal_related_range) = std::get<1>
      (
        calculate_non_terminal_permutation_vector_range(std::prev(node->branch.end()))
      );
    }
    return node->non_terminal_related_range;
  }
};

#endif
