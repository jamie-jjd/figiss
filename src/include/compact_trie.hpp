#ifndef TRIE_HPP_
#define TRIE_HPP_

#include <cstdint>
#include <fstream>
#include <map>
#include <memory>
#include <vector>

template <typename Character_Type>
struct Compact_Trie_Node
{
  using Compact_Trie_Node_Type = Compact_Trie_Node<Character_Type>;
  using Compact_Trie_Node_Pointer_Type = std::shared_ptr<Compact_Trie_Node>;
  using Compact_Trie_Node_Map_Type = std::map<Character_Type, Compact_Trie_Node_Pointer_Type>;

  int32_t edge_begin_;
  int32_t edge_end_;
  Compact_Trie_Node_Map_Type children_;

  Compact_Trie_Node (int32_t const edge_begin, int32_t const edge_end)
  {
    edge_begin_ = edge_begin;
    edge_end_ = edge_end;
  }

  std::pair<int32_t, int32_t> get_edge_range () const
  {
    return {edge_begin_, edge_end_};
  }

  void set_edge_range (int32_t const edge_begin, int32_t const edge_end)
  {
    edge_begin_ = edge_begin;
    edge_end_ = edge_end;
  }

  void insert
  (
    std::vector<Character_Type> const &text,
    int32_t const substring_begin,
    int32_t const substring_end
  )
  {
    if (substring_begin == substring_end)
    {
      children_[0] = std::make_shared<Compact_Trie_Node_Type>(substring_begin, substring_end);
      return;
    }

    if (children_.find(text[substring_begin]) == children_.cend())
    {
      children_[text[substring_begin]] = std::make_shared<Compact_Trie_Node_Type>(substring_begin, substring_end);
      children_[text[substring_begin]]->insert(text, substring_end, substring_end);
      return;
    }

    auto child_node_pointer {children_[text[substring_begin]]};
    auto const [edge_begin, edge_end] {child_node_pointer->get_edge_range()};
    auto edge_index {edge_begin};
    auto substring_index {substring_begin};
    while
    (
      (edge_index != edge_end) &&
      (substring_index != substring_end) &&
      (text[edge_index] == text[substring_index])
    )
    {
      ++edge_index;
      ++substring_index;
    }
    if (edge_index != edge_end)
    {
      auto new_child_node_pointer {std::make_shared<Compact_Trie_Node_Type>(edge_begin, edge_index)};
      child_node_pointer->set_edge_range(edge_index, edge_end);
      new_child_node_pointer->children_[text[edge_index]] = child_node_pointer;
      new_child_node_pointer->insert(text, substring_index, substring_end);
      children_[text[edge_begin]] = new_child_node_pointer;
    }
    else
    {
      child_node_pointer->insert(text, substring_index, substring_end);
    }
    return;
  }

  int32_t size () const
  {
    int32_t size_ {edge_end_ - edge_begin_};
    for (auto const &[character, child_node_pointer] : children_)
    {
      if (character != 0)
      {
        size_ += child_node_pointer->size();
      }
      else
      {
        ++size_;
      }
    }
    return size_;
  }

  friend std::ofstream& operator<< (std::ofstream &output_file_stream, Compact_Trie_Node const &compact_trie_node)
  {
    auto const &children {compact_trie_node.children_};
    if (!children.empty())
    {
      output_file_stream << '(' << compact_trie_node.edge_begin_ << ',' << compact_trie_node.edge_end_ << ")->(";
      auto const [first_character, first_child_node_pointer] {*(children.cbegin())};
      if (first_character != 0)
      {
        output_file_stream << first_character << ':';
        output_file_stream << *first_child_node_pointer;
      }
      else
      {
        output_file_stream << 0;
      }
      for (auto children_iterator {++(children.cbegin())}; children_iterator != children.cend(); ++children_iterator)
      {
        auto const [character, child_node_pointer] {*children_iterator};
        output_file_stream << ',' << character << ':';
        output_file_stream << *child_node_pointer;
      }
      output_file_stream << ')';
    }
    return output_file_stream;
  }
};

template <typename Character_Type>
struct Compact_Trie
{
  using Compact_Trie_Node_Type = Compact_Trie_Node<Character_Type>;
  using Compact_Trie_Node_Pointer_Type = std::shared_ptr<Compact_Trie_Node_Type>;

  Compact_Trie_Node_Pointer_Type root_;

  Compact_Trie ()
  {
    root_ = std::make_shared<Compact_Trie_Node_Type>(-1, -1);
  }

  void clear ()
  {
    root_ = std::make_shared<Compact_Trie_Node_Type>(-1, -1);
  }

  void insert
  (
    std::vector<Character_Type> const &text,
    int32_t const substring_begin,
    int32_t const substring_end
  )
  {
    root_->insert(text, substring_begin, substring_end);
  }

  int32_t size () const
  {
    return root_->size();
  }

  friend std::ofstream& operator<< (std::ofstream &output_file_stream, Compact_Trie const &compact_trie)
  {
    output_file_stream << *(compact_trie.root_) << '\n';
    return output_file_stream;
  }
};

#endif
