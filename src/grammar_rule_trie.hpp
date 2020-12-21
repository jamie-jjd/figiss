#ifndef GRAMMAR_RULE_TRIE_HPP_
#define GRAMMAR_RULE_TRIE_HPP_

#include <map>
#include <memory>

#include <sdsl/int_vector.hpp>

template <typename T>
void f (T)
{
  std::cout << __PRETTY_FUNCTION__ << '\n';
}

template
<
  typename rule_size_vector_type,
  typename sl_type_vector_type,
  typename grammar_rule_begin_distance_vector_iterator_type
>
void calculate_rule_size_vector
(
  rule_size_vector_type &rule_size_vector,
  sl_type_vector_type const &sl_type_vector,
  grammar_rule_begin_distance_vector_iterator_type grammar_rule_begin_distance_vector_begin,
  grammar_rule_begin_distance_vector_iterator_type grammar_rule_begin_distance_vector_end
)
{
  rule_size_vector.resize
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
    auto grammar_rule_size_it {std::begin(rule_size_vector)};
    grammar_rule_size_it != std::end(rule_size_vector);
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
  typename rule_size_vector_type,
  typename rule_vector_type,
  typename grammar_rule_begin_distance_vector_iterator_type
>
void calculate_rule_vector
(
  text_type const &text,
  rule_size_vector_type const &rule_size_vector,
  rule_vector_type &rule_vector,
  grammar_rule_begin_distance_vector_iterator_type begin_distance_begin
)
{
  rule_vector.resize
  (
    std::accumulate
    (
      std::begin(rule_size_vector),
      std::end(rule_size_vector),
      0
    )
  );
  auto rule_vector_it {std::begin(rule_vector)};
  auto begin_distance_it {begin_distance_begin};
  for (auto const &grammar_rule_size : rule_size_vector)
  {
    auto grammar_rule_it {std::next(std::begin(text), *begin_distance_it)};
    auto grammar_rule_end {std::next(std::begin(text), *begin_distance_it + grammar_rule_size)};
    while (grammar_rule_it != grammar_rule_end)
    {
      *rule_vector_it = *grammar_rule_it;
      ++rule_vector_it;
      ++grammar_rule_it;
    }
    ++begin_distance_it;
  }
  return;
}

template
<
  typename rule_size_vector_type = sdsl::int_vector<>,
  typename rule_vector_type = sdsl::int_vector<8>
>
struct grammar_rule_trie
{
  using difference_type = typename rule_vector_type::difference_type;
  using size_type = typename rule_vector_type::size_type;
  using character_type = typename rule_vector_type::value_type;

  struct trie_node;
  using trie_node_pointer_type = std::shared_ptr<trie_node>;
  using branch_type = std::map<character_type, trie_node_pointer_type>;

  struct trie_node
  {
    difference_type edge_begin_distance;
    difference_type edge_end_distance;
    branch_type branch;
    size_type non_terminal_lower_bound;
    size_type non_terminal_upper_bound;

    trie_node () = default;

    trie_node
    (
      difference_type edge_begin_distance_,
      difference_type edge_end_distance_,
      size_type non_terminal
    )
    : edge_begin_distance {edge_begin_distance_},
      edge_end_distance {edge_end_distance_},
      non_terminal_lower_bound {non_terminal},
      non_terminal_upper_bound {non_terminal}
    {
    }
  };

  rule_size_vector_type rule_size_vector;
  rule_vector_type rule_vector;
  trie_node_pointer_type root;
  trie_node_pointer_type colex_root;

  grammar_rule_trie ()
  : root {std::make_shared<trie_node>()},
    colex_root {std::make_shared<trie_node>()}
  {

  }

  template
  <
    typename text_iterator_type
  >
  void insert_grammar_rule
  (
    text_iterator_type text_begin,
    text_iterator_type text_it,
    text_iterator_type text_last,
    int difference
  )
  {
    auto node {root};
    if (difference < 0)
    {
      node = colex_root;
    }
    while (text_it != text_last)
    {
      auto character {*text_it};
      if (node->branch.find(character) == std::end(node->branch))
      {
        node->branch[character] = std::make_shared<trie_node>
        (
          std::distance(text_begin, text_it),
          std::distance(text_begin, text_last),
          1
        );
        return;
      }
      auto child {node->branch[character]};
      auto edge_it {std::next(text_begin, child->edge_begin_distance)};
      auto edge_end {std::next(text_begin, child->edge_end_distance)};
      while
      (
        (text_it != text_last)
        &&
        (edge_it != edge_end)
        &&
        (*text_it == *edge_it)
      )
      {
        text_it += difference;
        edge_it += difference;
      }
      if (edge_it == edge_end)
      {
        node = child;
      }
      else
      {
        auto edge_it_distance {std::distance(text_begin, edge_it)};
        auto internal_node {std::make_shared<trie_node>(child->edge_begin_distance, edge_it_distance, 0)};
        child->edge_begin_distance = edge_it_distance;
        internal_node->branch[*edge_it] = child;
        if (text_it != text_last)
        {
          internal_node->branch[*text_it] = std::make_shared<trie_node>
          (
            std::distance(text_begin, text_it),
            std::distance(text_begin, text_last),
            1
          );
        }
        else
        {
          internal_node->non_terminal_lower_bound = 1;
          internal_node->non_terminal_upper_bound = 1;
        }
        node->branch[character] = internal_node;
        return;
      }
    }
    return;
  }

  void calculate_non_terminal_range ()
  {
    calculate_non_terminal_range(root);
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
      if (node->non_terminal_lower_bound == 0)
      {
        node->non_terminal_lower_bound = first_child->non_terminal_lower_bound;
      }
      node->non_terminal_upper_bound = last_child->non_terminal_upper_bound;
    }
    return;
  }

  template
  <
    typename text_iterator_type
  >
  void print
  (
    text_iterator_type text_begin,
    trie_node_pointer_type node,
    int difference
  )
  {
    static size_t depth {0};
    if (node != nullptr)
    {
      std::cout << depth << ": ";
      if (node->edge_begin_distance != node->edge_end_distance)
      {
        auto it {std::next(text_begin, node->edge_begin_distance)};
        auto end {std::next(text_begin, node->edge_end_distance)};
        while (it != end)
        {
          std::cout << *it;
          it += difference;
        }
        std::cout << ' ';
      }
      std::cout << '(' << node->non_terminal_lower_bound << ", " << node->non_terminal_upper_bound << ")\n";
      for (auto character_child_pair : node->branch)
      {
        ++depth;
        print(text_begin, std::get<1>(character_child_pair), difference);
        --depth;
      }
    }
    return;
  }
};

#endif
