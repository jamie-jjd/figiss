#pragma once

#include <sys/wait.h>

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <map>
#include <memory>
#include <random>

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

void GeneratePrefix
(
  std::filesystem::path const& byte_text_path,
  uint64_t const size_in_megabytes
)
{
  sdsl::int_vector<8> byte_text;
  sdsl::load_vector_from_file(byte_text, byte_text_path);
  auto prefix_path
  {
    std::filesystem::path
    {
      byte_text_path.filename().string()
      + "." + std::to_string(size_in_megabytes) + "mb"
    }
  };
  std::fstream out {prefix_path, std::ios_base::out | std::ios_base::trunc};
  for (uint64_t i {}; i != (size_in_megabytes * (1ULL << 20)); ++i)
  {
    out << byte_text[i];
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
  std::cout
  << "\n\tdecompress " << std::filesystem::canonical(compressed_corpus_path).string() << "/* to "
  << std::filesystem::canonical(corpus_path).string() << "/*\n"
  << "\t(xz and p7zip are required for decompression)\n";
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

void GenerateGenerationalDnaSequence
(
  uint64_t const size,
  uint64_t const amount_copies,
  uint64_t const mutative_rate,
  std::filesystem::path const& path
)
{
  if (!path.parent_path().empty() && !std::filesystem::exists(path.parent_path()))
  {
    std::filesystem::create_directories(path.parent_path());
  }
  auto out {std::ofstream{path}};
  auto engine {std::mt19937_64(std::random_device{}())};
  auto const bases {std::array<char, 4>({'A', 'C', 'G', 'T'})};
  auto dist_base {std::uniform_int_distribution<uint64_t>(0, 3)};
  auto dist_other_base {std::uniform_int_distribution<uint64_t>(1, 3)};
  auto dist_mutation {std::uniform_int_distribution<uint64_t>(0, 99)};
  auto base_sequence {std::string()};
  for (uint64_t i {}; i != size; ++i)
  {
    base_sequence.push_back(dist_base(engine));
  }
  for (uint64_t i {}; i != amount_copies; ++i)
  {
    for (uint64_t j {}; j != size; ++j)
    {
      auto dice {dist_mutation(engine)};
      if (dice < mutative_rate)
      {
        if (dice & 1)
        {
          base_sequence.push_back((base_sequence[j] + dist_other_base(engine)) % 4);
        }
      }
      else
      {
        base_sequence.push_back(base_sequence[j]);
      }
    }
  }
  for (auto& base : base_sequence)
  {
    base = bases[base];
  }
  out << base_sequence;
  return;
}

void PrintTextParameters (std::filesystem::path const& byte_text_path)
{
  auto parameter_path {std::filesystem::path{"../data/parameter/" + byte_text_path.filename().string() + "/parameter"}};
  if (!std::filesystem::exists(parameter_path.parent_path()))
  {
    std::filesystem::create_directories(parameter_path.parent_path());
  }
  auto out {std::ofstream(parameter_path)};
  std::cout << "write text parameters to " + std::filesystem::canonical(parameter_path).string() << "\n";
  {
    sdsl::int_vector<8> byte_text;
    sdsl::load_vector_from_file(byte_text, byte_text_path);
    sdsl::append_zero_symbol(byte_text);
    out << "n:" << std::size(byte_text) << "\n";
    {
      sdsl::int_vector<> buffer;
      sdsl::qsufsort::construct_sa(buffer, byte_text);
      uint64_t r {};
      uint64_t prev_character {std::size(byte_text)};
      for (auto it {std::begin(buffer)}; it != std::end(buffer); ++it)
      {
        if (*it != 0)
        {
          *it = byte_text[*it - 1];
        }
        if (*it != prev_character)
        {
          prev_character = *it;
          ++r;
        }
      }
      out << "r:" << r << "\n";
    }
    {
      std::set<uint8_t> byte_alphabet;
      for (auto const& byte : byte_text)
      {
        byte_alphabet.insert(byte);
      }
      out << "\u03C3:" << std::size(byte_alphabet) << "\n";
    }
  }
}
}
