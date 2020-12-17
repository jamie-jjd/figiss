#ifndef GRAMMAR_RULE_TRIE_HPP_
#define GRAMMAR_RULE_TRIE_HPP_

#include <map>
#include <memory>

#include <sdsl/int_vector.hpp>

#include "sl_type.hpp"

template <typename T>
void f (T)
{
  std::cout << __PRETTY_FUNCTION__ << '\n';
}

template
<
  typename grammar_rule_size_vector_type,
  typename sl_type_vector_type,
  typename grammar_rule_begin_distance_vector_iterator_type
>
void calculate_grammar_rule_size_vector
(
  grammar_rule_size_vector_type &grammar_rule_size_vector,
  sl_type_vector_type const &sl_type_vector,
  grammar_rule_begin_distance_vector_iterator_type grammar_rule_begin_distance_vector_begin,
  grammar_rule_begin_distance_vector_iterator_type grammar_rule_begin_distance_vector_end
)
{
  grammar_rule_size_vector.resize
  (
    std::distance
    (
      grammar_rule_begin_distance_vector_begin,
      grammar_rule_begin_distance_vector_end
    )
  );
  auto begin_distance_it {grammar_rule_begin_distance_vector_begin};
  for
  (
    auto grammar_rule_size_it {std::begin(grammar_rule_size_vector)};
    grammar_rule_size_it != std::end(grammar_rule_size_vector);
    ++grammar_rule_size_it
  )
  {
    auto sl_type_vector_current_begin {std::next(std::begin(sl_type_vector), *begin_distance_it)};
    auto sl_type_vector_it {sl_type_vector_current_begin};
    do
    {
      ++sl_type_vector_it;
    }
    while (!is_leftmost_s_type(sl_type_vector_it));
    *grammar_rule_size_it = std::distance(sl_type_vector_current_begin, sl_type_vector_it);
    ++begin_distance_it;
  }
  return;
}

template
<
  typename text_type,
  typename grammar_rule_size_vector_type,
  typename grammar_rule_vector_type,
  typename grammar_rule_begin_distance_vector_iterator_type
>
void calculate_grammar_rule_vector
(
  text_type const &text,
  grammar_rule_size_vector_type const &grammar_rule_size_vector,
  grammar_rule_vector_type &grammar_rule_vector,
  grammar_rule_begin_distance_vector_iterator_type begin_distance_begin
)
{
  grammar_rule_vector.resize
  (
    std::accumulate
    (
      std::begin(grammar_rule_size_vector),
      std::end(grammar_rule_size_vector),
      0
    )
  );
  auto grammar_rule_vector_it {std::begin(grammar_rule_vector)};
  auto begin_distance_it {begin_distance_begin};
  for (auto const &grammar_rule_size : grammar_rule_size_vector)
  {
    auto grammar_rule_it {std::next(std::begin(text), *begin_distance_it)};
    auto grammar_rule_end {std::next(std::begin(text), *begin_distance_it + grammar_rule_size)};
    while (grammar_rule_it != grammar_rule_end)
    {
      *grammar_rule_vector_it = *grammar_rule_it;
      ++grammar_rule_vector_it;
      ++grammar_rule_it;
    }
    ++begin_distance_it;
  }
  return;
}

template
<
  typename grammar_rule_size_vector_type = sdsl::int_vector<>,
  typename grammar_rule_vector_type = sdsl::int_vector<8>
>
struct grammar_rule_trie
{
  using difference_type = typename grammar_rule_vector_type::difference_type;
  using size_type = typename grammar_rule_vector_type::size_type;
  using character_type = typename grammar_rule_vector_type::value_type;

  struct trie_node
  {
    using trie_node_pointer_type = std::shared_ptr<trie_node>;
    using branch_type = std::map<character_type, trie_node_pointer_type>;

    difference_type edge_begin_distance;
    difference_type edge_end_distance;
    branch_type branch;
    size_type non_terminal_begin_distance;
    size_type non_terminal_end_distance;

    trie_node () = default;

    trie_node
    (
      difference_type edge_begin_distance_,
      difference_type edge_end_distance_,
      size_type non_terminal_value
    )
    : edge_begin_distance {edge_begin_distance_},
      edge_end_distance {edge_end_distance_},
      non_terminal_begin_distance {non_terminal_value},
      non_terminal_end_distance {non_terminal_value}
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

  template
  <
    typename text_type,
    typename sl_type_vector_type,
    typename grammar_rule_begin_distance_vector_iterator_type
  >
  grammar_rule_trie
  (
    text_type const &text,
    sl_type_vector_type const &sl_type_vector,
    grammar_rule_begin_distance_vector_iterator_type grammar_rule_begin_distance_vector_begin,
    grammar_rule_begin_distance_vector_iterator_type grammar_rule_begin_distance_vector_end
  )
  {
    calculate_grammar_rule_size_vector
    (
      grammar_rule_size_vector,
      sl_type_vector,
      grammar_rule_begin_distance_vector_begin,
      grammar_rule_begin_distance_vector_end
    );
    calculate_grammar_rule_vector
    (
      text,
      grammar_rule_size_vector,
      grammar_rule_vector,
      grammar_rule_begin_distance_vector_begin
    );
    insert_grammar_rule_vector();
    calculate_non_terminal_range();
    insert_reverse_grammar_rule_vector();
    calculate_non_terminal_distance_range_and_permutation_vector();
  }

  void insert_grammar_rule_vector ()
  {
    grammar_rule_trie_root = std::make_shared<trie_node>();
    size_type non_terminal {1};
    auto grammar_rule_size_it {std::begin(grammar_rule_size_vector)};
    auto grammar_rule_begin {std::begin(grammar_rule_vector)};
    while (grammar_rule_begin != std::end(grammar_rule_vector))
    {
      auto grammar_rule_it {grammar_rule_begin};
      auto grammar_rule_end {std::next(grammar_rule_begin, *grammar_rule_size_it)};
      auto node {grammar_rule_trie_root};
      while (true)
      {
        auto character {*grammar_rule_it};
        if (node->branch.find(character) == std::end(node->branch))
        {
          node->branch[character] = std::make_shared<trie_node>
          (
            std::distance(std::begin(grammar_rule_vector), grammar_rule_it),
            std::distance(std::begin(grammar_rule_vector), grammar_rule_end),
            non_terminal
          );
          break;
        }
        else
        {
          auto child {node->branch[character]};
          auto edge_it {std::next(std::begin(grammar_rule_vector), child->edge_begin_distance)};
          auto edge_end {std::next(std::begin(grammar_rule_vector), child->edge_end_distance)};
          while
          (
            (edge_it != edge_end)
            &&
            (*grammar_rule_it == *edge_it)
          )
          {
            ++grammar_rule_it;
            ++edge_it;
          }
          if (edge_it == edge_end)
          {
            node = child;
          }
          else
          {
            auto edge_it_distance {std::distance(std::begin(grammar_rule_vector), edge_it)};
            auto internal_node {std::make_shared<trie_node>(child->edge_begin_distance, edge_it_distance, 0)};
            child->edge_begin_distance = edge_it_distance;
            node->branch[character] = internal_node;
            internal_node->branch[*edge_it] = child;
            internal_node->branch[*grammar_rule_it] = std::make_shared<trie_node>
            (
              std::distance(std::begin(grammar_rule_vector), grammar_rule_it),
              std::distance(std::begin(grammar_rule_vector), grammar_rule_end),
              non_terminal
            );
            break;
          }
        }
      }
      ++non_terminal;
      ++grammar_rule_size_it;
      grammar_rule_begin = grammar_rule_end;
    }
    return;
  }

  void calculate_non_terminal_range ()
  {
    calculate_non_terminal_range(grammar_rule_trie_root);
    return;
  }

  void calculate_non_terminal_range (trie_node_pointer_type node)
  {
    if (!node->branch.empty())
    {
      auto branch_it {std::begin(node->branch)};
      auto branch_end {std::end(node->branch)};
      while (branch_it != branch_end)
      {
        calculate_non_terminal_range(std::get<1>(*branch_it));
        ++branch_it;
      }
      auto first_child {std::get<1>(*std::begin(node->branch))};
      auto last_child {std::get<1>(*std::prev(std::end(node->branch)))};
      if (node->non_terminal_begin_distance == 0)
      {
        node->non_terminal_begin_distance = first_child->non_terminal_begin_distance;
      }
      node->non_terminal_end_distance = last_child->non_terminal_end_distance;
    }
    return;
  }

  void insert_reverse_grammar_rule_vector ()
  {
    reverse_grammar_rule_trie_root = std::make_shared<trie_node>();
    size_type non_terminal {std::size(grammar_rule_size_vector)};
    auto grammar_rule_size_rit {std::prev(std::end(grammar_rule_size_vector))};
    auto grammar_rule_rbegin {std::prev(std::end(grammar_rule_vector))};
    while (grammar_rule_rbegin != std::prev(std::begin(grammar_rule_vector)))
    {
      auto grammar_rule_rit {grammar_rule_rbegin};
      auto grammar_rule_rend {std::prev(grammar_rule_rbegin, *grammar_rule_size_rit)};
      auto node {reverse_grammar_rule_trie_root};
      while (true)
      {
        auto character {*grammar_rule_rit};
        if (node->branch.find(character) == std::end(node->branch))
        {
          node->branch[character] = std::make_shared<trie_node>
          (
            std::distance(std::begin(grammar_rule_vector), grammar_rule_rit),
            std::distance(std::begin(grammar_rule_vector), grammar_rule_rend),
            non_terminal
          );
          break;
        }
        else
        {
          auto child {node->branch[character]};
          auto edge_rit {std::next(std::begin(grammar_rule_vector), child->edge_begin_distance)};
          auto edge_rend {std::next(std::begin(grammar_rule_vector), child->edge_end_distance)};
          while
          (
            (grammar_rule_rit != grammar_rule_rend)
            &&
            (edge_rit != edge_rend)
            &&
            (*grammar_rule_rit == *edge_rit)
          )
          {
            --grammar_rule_rit;
            --edge_rit;
          }
          if (edge_rit == edge_rend)
          {
            node = child;
          }
          else
          {
            auto edge_rit_distance {std::distance(std::begin(grammar_rule_vector), edge_rit)};
            auto internal_node {std::make_shared<trie_node>(child->edge_begin_distance, edge_rit_distance, 0)};
            child->edge_begin_distance = edge_rit_distance;
            node->branch[character] = internal_node;
            internal_node->branch[*edge_rit] = child;
            if (grammar_rule_rit != grammar_rule_rend)
            {
              internal_node->branch[*grammar_rule_rit] = std::make_shared<trie_node>
              (
                std::distance(std::begin(grammar_rule_vector), grammar_rule_rit),
                std::distance(std::begin(grammar_rule_vector), grammar_rule_rend),
                non_terminal
              );
            }
            else
            {
              internal_node->non_terminal_begin_distance = non_terminal;
              internal_node->non_terminal_end_distance = non_terminal;
            }
            break;
          }
        }
      }
      --non_terminal;
      --grammar_rule_size_rit;
      grammar_rule_rbegin = grammar_rule_rend;
    }
    return;
  }

  void calculate_non_terminal_distance_range_and_permutation_vector ()
  {
    non_terminal_permutation_vector.resize(std::size(grammar_rule_size_vector));
    calculate_non_terminal_distance_range_and_permutation_vector(reverse_grammar_rule_trie_root);
    return;
  }

  void calculate_non_terminal_distance_range_and_permutation_vector (trie_node_pointer_type node)
  {
    static size_type distance {std::size(non_terminal_permutation_vector)};
    if (node->non_terminal_end_distance != 0)
    {
      non_terminal_permutation_vector[--distance] = node->non_terminal_end_distance;
      node->non_terminal_begin_distance = node->non_terminal_end_distance = distance;
    }
    if (!node->branch.empty())
    {
      auto branch_rit = std::prev(std::end(node->branch));
      auto branch_rend = std::prev(std::begin(node->branch));
      while (branch_rit != branch_rend)
      {
        calculate_non_terminal_distance_range_and_permutation_vector(std::get<1>(*branch_rit));
        --branch_rit;
      }
      auto last_child {std::get<1>(*std::prev(std::end(node->branch)))};
      auto first_child {std::get<1>(*std::begin(node->branch))};
      if (node->non_terminal_end_distance == 0)
      {
        node->non_terminal_end_distance = last_child->non_terminal_end_distance;
      }
      node->non_terminal_begin_distance = first_child->non_terminal_begin_distance;
    }
    return;
  }

  void print_trie (trie_node_pointer_type node, bool is_reversed)
  {
    static size_t depth {0};
    if (node != nullptr)
    {
      std::cout << depth << ": ";
      if (node->edge_begin_distance != node->edge_end_distance)
      {
        auto it {std::next(std::begin(grammar_rule_vector), node->edge_begin_distance)};
        auto end {std::next(std::begin(grammar_rule_vector), node->edge_end_distance)};
        while (it != end)
        {
          std::cout << *it;
          if (is_reversed)
          {
            --it;
          }
          else
          {
            ++it;
          }
        }
        std::cout << ' ';
      }
      std::cout << '(' << node->non_terminal_begin_distance << ", " << node->non_terminal_end_distance << ")\n";
      for (auto character_child_pair : node->branch)
      {
        ++depth;
        print_trie(std::get<1>(character_child_pair), is_reversed);
        --depth;
      }
    }
    return;
  }
};

#endif
