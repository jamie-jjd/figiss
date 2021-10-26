#pragma once

#include <sys/wait.h>

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <map>
#include <memory>

namespace figiss
{
template <typename File, typename Container>
void Print
(
  Container const& container,
  File& file,
  int8_t const step = 1,
  std::string const& separator = " ",
  std::string const& endmarker = "\n"
)
{
  if (std::size(container) != 0)
  {
    auto it {std::begin(container)};
    auto prev_end {std::prev(std::end(container))};
    if (step < 0)
    {
      it = std::prev(std::end(container));
      prev_end = std::begin(container);
    }
    while (it != prev_end)
    {
      file << static_cast<uint64_t>(*it) << separator;
      std::advance(it, step);
    }
    file << static_cast<uint64_t>(*prev_end) << endmarker;
  }
  return;
}

template <typename File, typename Iterator>
void Print
(
  Iterator first,
  Iterator last,
  File& file,
  int8_t const step = 1,
  std::string const& separator = " ",
  std::string const& endmarker = "\n"
)
{
  if (std::distance(first, last) * step > 0)
  {
    auto it {first};
    auto prev_last {std::prev(last, step)};
    while (it != prev_last)
    {
      file << static_cast<uint64_t>(*it) << separator;
      std::advance(it, step);
    }
    file << static_cast<uint64_t>(*prev_last) << endmarker;
  }
  return;
}

template <typename Size>
std::string ProperSizeRepresentation (Size const size)
{
  std::stringstream stringstream;
  stringstream << std::fixed << std::setprecision(2);
  if (size >= (1ULL << 30))
  {
    stringstream << (static_cast<double>(size) / (1ULL << 30));
    return (stringstream.str() + "Gi");
  }
  else if (size >= (1ULL << 20))
  {
    stringstream << (static_cast<double>(size) / (1ULL << 20));
    return (stringstream.str() + "Mi");
  }
  else if (size >= (1ULL << 10))
  {
    stringstream << (static_cast<double>(size) / (1ULL << 10));
    return (stringstream.str() + "Ki");
  }
  else
  {
    stringstream << static_cast<double>(size);
    return stringstream.str();
  }
  return "0.00";
}

class SpaceNode
{
public:

  SpaceNode () = default;
  SpaceNode (std::string const& name, uint64_t const size_in_bytes = 0);

  void AccumalateSizeInBytes (uint64_t const size_in_bytes);
  void AddChild (std::shared_ptr<SpaceNode> child);
  void AddLeaf (std::string const& name, uint64_t const size_in_bytes);

  inline uint64_t GetSizeInBytes () const
  {
    return size_in_bytes_;
  }

  friend std::ostream& operator<<
  (
    std::ostream& out,
    std::pair<std::shared_ptr<SpaceNode>, bool> pair
  );

private:

  std::string name_;
  uint64_t size_in_bytes_;
  std::deque<std::shared_ptr<SpaceNode>> children_;

};

SpaceNode::SpaceNode (std::string const& name, uint64_t const size_in_bytes)
: name_ {name},
  size_in_bytes_ {size_in_bytes},
  children_ {}
{
}

void SpaceNode::AccumalateSizeInBytes (uint64_t const size_in_bytes)
{
  size_in_bytes_ += size_in_bytes;
  return;
}

void SpaceNode::AddChild (std::shared_ptr<SpaceNode> child)
{
  children_.emplace_back(child);
  size_in_bytes_ += child->GetSizeInBytes();
  return;
}

void SpaceNode::AddLeaf (std::string const& name, uint64_t const size_in_bytes)
{
  auto node {std::make_shared<SpaceNode>(name, size_in_bytes)};
  children_.emplace_back(node);
  size_in_bytes_ += size_in_bytes;
  return;
}

std::ostream& operator<<
(
  std::ostream& out,
  std::pair<std::shared_ptr<SpaceNode>, bool> pair
)
{
  std::deque<std::pair<std::shared_ptr<SpaceNode>, uint64_t>> nodes;
  auto root {std::get<0>(pair)};
  auto is_proper {std::get<1>(pair)};
  nodes.emplace_back(root, 0);
  while (!nodes.empty())
  {
    auto const node {std::get<0>(nodes.back())};
    auto const depth {std::get<1>(nodes.back())};
    nodes.pop_back();
    std::string whitespaces(depth * 2, ' ');
    out << whitespaces << node->name_ << ":";
    if (is_proper)
    {
      out << ProperSizeRepresentation(node->size_in_bytes_);
    }
    else
    {
      out << node->size_in_bytes_;
    }
    out << "B\n";
    for (auto it {std::rbegin(node->children_)}; it != std::rend(node->children_); ++it)
    {
      nodes.emplace_back(*it, depth + 1);
    }
  }
  return out;
}

void DecompressCompressedCorpus
(
  std::filesystem::path const& compressed_corpus_path,
  std::filesystem::path const& corpus_path
)
{
  if (!std::filesystem::exists(corpus_path))
  {
    std::filesystem::create_directories(corpus_path);
  }
  for (auto const& entry : std::filesystem::directory_iterator(compressed_corpus_path))
  {
    if (entry.is_regular_file())
    {
      auto compressed_text_path {entry.path()};
      auto decompressed_text_to_path {corpus_path / entry.path().stem()};
      if (!std::filesystem::exists(decompressed_text_to_path))
      {
        if (compressed_text_path.extension() == ".7z")
        {
          {
            auto code {system(("p7zip -k -d " + compressed_text_path.string()).c_str())};
            if (WIFSIGNALED(code) && (WTERMSIG(code) == SIGINT || WTERMSIG(code) == SIGQUIT))
            {
              throw std::runtime_error
              (
                "\033[31mfailed at decompressing " +
                std::filesystem::canonical(compressed_text_path).string() +
                "\033[0m\n"
              );
            }
          }
          {
            auto decompressed_text_from_path {decompressed_text_to_path.filename()};
            auto code {system(("mv " + decompressed_text_from_path.string() + " " + decompressed_text_to_path.string()).c_str())};
            if (WIFSIGNALED(code) && (WTERMSIG(code) == SIGINT || WTERMSIG(code) == SIGQUIT))
            {
              throw std::runtime_error("\033[31mfailed at moving decompressed text");
            }
          }
        }
        else if (compressed_text_path.extension() == ".xz")
        {
          {
            auto code {system(("xz -kdv -T0 " + compressed_text_path.string()).c_str())};
            if (WIFSIGNALED(code) && (WTERMSIG(code) == SIGINT || WTERMSIG(code) == SIGQUIT))
            {
              throw std::runtime_error
              (
                "\033[31mfailed at decompressing " +
                std::filesystem::canonical(compressed_text_path).string() +
                "\033[0m\n"
              );
            }
            {
              auto decompressed_text_from_path {compressed_corpus_path / compressed_text_path.stem()};
              auto code {system(("mv " + decompressed_text_from_path.string() + " " + decompressed_text_to_path.string()).c_str())};
              if (WIFSIGNALED(code) && (WTERMSIG(code) == SIGINT || WTERMSIG(code) == SIGQUIT))
              {
                throw std::runtime_error("\033[31mfailed at moving decompressed text");
              }
            }
          }
        }
      }
    }
  }
  return;
}
}
