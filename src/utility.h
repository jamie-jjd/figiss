#ifndef UTILITY_HPP_
#define UTILITY_HPP_

#include <chrono>
#include <filesystem>
#include <fstream>
#include <functional>
#include <map>
#include <memory>
#include <random>
#include <string>
#include <unordered_map>
#include <unordered_set>

// #include <cmath>

namespace fs = std::filesystem;
using namespace std::chrono;

namespace grammar_compressed_index
{
inline void generate_pattern
(
  fs::path const text_path,
  uint64_t const pattern_amount,
  uint64_t const pattern_size
)
{
  std::ifstream text_file {text_path};
  sdsl::int_vector<> text;
  text.load(text_file);
  auto path {fs::path{"../data/pattern"} / text_path.parent_path().filename()};
  if (!fs::exists(path))
  {
    fs::create_directories(path);
  }
  if (pattern_size < std::size(text))
  {
    path /= fs::path
    {
      text_path.filename().string()
      + "." + std::to_string(pattern_amount)
      + "." + std::to_string(pattern_size)
    };
    std::ofstream pattern_file {path};
    std::mt19937 engine {std::random_device{}()};
    std::uniform_int_distribution<uint64_t> distribution(0, std::size(text) - pattern_size);
    auto random_begin_offset {std::bind(distribution, engine)};
    sdsl::int_vector<> pattern(pattern_size, 0, text.width());
    pattern_file << pattern_amount << "\n";
    for (uint64_t i {0}; i != pattern_amount; ++i)
    {
      auto text_it {std::next(std::begin(text), random_begin_offset())};
      auto pattern_it {std::begin(pattern)};
      auto pattern_end {std::end(pattern)};
      while (pattern_it != pattern_end)
      {
        *pattern_it = *text_it;
        ++text_it;
        ++pattern_it;
      }
      pattern.serialize(pattern_file);
    }
  }
  else
  {
    path /= fs::path{text_path.filename().string() + ".1." + std::to_string(std::size(text))};
    std::ofstream pattern_file {path};
    std::cout
    << "pattern size >= text size\n"
    << "generated pattern is the text\n";
    pattern_file << "1\n";
    text.serialize(pattern_file);
  }
  return;
}

enum class record_key
{
#if defined(LOG_GCI_COUNT) || defined(LOG_VERBOSE_GCI_COUNT)
  GCI_COUNT,
#endif
#if defined(LOG_VERBOSE_GCI_COUNT)
  CALCULATE_SL_FACTOR,
  LOOKUP_GRAMMAR_RULE,
  WAVELET_TREE_OPERATIONS_L,
  WAVELET_TREE_OPERATIONS_S,
#endif
#ifdef LOG_FMI_COUNT
  FMI_COUNT,
#endif
  DUMMY
};

struct timer
{
  using clock_ = steady_clock;
  using record = std::tuple<std::string, clock_::duration, clock_::time_point>;

  static std::unordered_map<record_key, record> records;

  static void register_record
  (
    record_key const key,
    std::string const &category
  )
  {
    records[key] = {category, {}, {}};
    return;
  }

  static void resume (record_key const key)
  {
    std::get<2>(records[key]) = clock_::now();
    return;
  }

  static void pause (record_key const key)
  {
    auto &record {records[key]};
    std::get<1>(record) += (clock_::now() - std::get<2>(record));
    return;
  }

  static void print (fs::path const &records_path)
  {
    std::ofstream records_file {records_path};
    fout << "category,seconds\n";
    for (auto const &record_pair : records)
    {
      auto const &record {std::get<1>(record_pair)};
      records_file
      << std::get<0>(record)
      << ","
      << duration_cast<duration<double>>(std::get<1>(record)).count()
      << "\n";
    }
    records.clear();
    return;
  }
};

std::unordered_map<record_key, timer::record> timer::records {};

void generate_canonical_text (fs::path text_path)
{
  sdsl::int_vector<8> text;
  sdsl::load_vector_from_file(text, text_path.string());
  fs::path path {fs::path("../input/sdsl") / text_path.parent_path().filename()};
  if (!fs::exists(path))
  {
    fs::create_directories(path);
  }
  std::ofstream fout {path / (text_path.filename().string() + ".sdsl")};
  sdsl::int_vector<8> codebook(256, 0);
  for (auto const &character : text)
  {
    codebook[character] = 1;
  }
  std::partial_sum(std::begin(codebook), std::end(codebook), std::begin(codebook));
  for (auto &character : text)
  {
    character = codebook[character];
  }
  sdsl::append_zero_symbol(text);
  sdsl::util::bit_compress(text);
  text.serialize(fout);
  return;
}

template <typename text_type>
void print_number_runs
(
  std::ofstream &fout,
  text_type &text
)
{
  uint64_t prev_character {0};
  uint64_t number_runs {0};
  for (auto const &character : text)
  {
    if (character != prev_character)
    {
      ++number_runs;
      prev_character = character;
    }
  }
  fout << "number_runs(r)," << number_runs << "\n";
  fout << "compression_ratio(r/n)" << (static_cast<double>(number_runs) / std::size(text));
  std::cout << "number_runs(r)," << number_runs << "\n";
  std::cout << "compression_ratio(r/n)" << (static_cast<double>(number_runs) / std::size(text));
  return;
}

template <typename text_type>
void print_alphabet_size
(
  std::ofstream &fout,
  text_type &text
)
{
  std::unordered_set<uint64_t> alphabet;
  for (auto const &character : text)
  {
    alphabet.insert(character);
  }
  fout << "alphabet_size(sigma)," << std::size(alphabet) << "\n";
  std::cout << "alphabet_size(sigma)," << std::size(alphabet) << "\n";
  return;
}

template <typename text_type>
void print_0th_entropy
(
  std::ofstream &fout,
  text_type &text
)
{
  std::map<uint64_t, uint64_t> alphabet_count;
  for (auto const &character : text)
  {
    auto it {alphabet_count.find(character)};
    if (it != alphabet_count.end())
    {
      ++(*it);
    }
    else
    {
      alphabet_count[character] = 1;
    }
  }
  double entropy {};
  auto log_text_size {std::log(std::size(text))};
  for (auto const &character_count : alphabet_count)
  {
    auto count {static_cast<double>(std::get<1>(character_count))};
    entropy += (count * (log_text_size - std::log(count)));
  }
  entropy /= std::size(text);
  fout << "0th_entropy(H_0)," << entropy << "\n";
  fout << "compression_ratio(sigma/H_0)," << (entropy / std::size(alphabet_count)) << "\n";
  std::cout << "0th_entropy(H_0)," << entropy << "\n";
  std::cout << "compression_ratio(sigma/H_0)," << (entropy / std::size(alphabet_count)) << "\n";
  return;
}

template <typename text_type>
struct k_mer_trie
{
  struct node
  {
    uint64_t edge_begin_offset;
    uint64_t edge_end_offset;
    std::map<uint64_t, std::shared_ptr<node>> branches;
    uint64_t context_size;
    std::map<uint64_t, uint64_t> alphabet_count;

    node () = default;
    node
    (
      uint64_t begin_offset,
      uint64_t end_offset
    )
    : edge_begin_offset {begin_offset},
      edge_end_offset {end_offset},
      branches {},
      context_size {},
      alphabet_count {}
    {
    }
  };

  std::shared_ptr<node> root;
  std::unordered_set<uint64_t> alphabet;
  double entropy;

  k_mer_trie (text_type const &text, uint64_t k)
  : root {std::make_shared<node>()},
    alphabet {},
    entropy {}
  {
    auto text_begin {std::begin(text)};
    auto text_first {text_begin};
    auto text_last {std::prev(std::end(text), k)};
    for (auto text_it {text_first}; text_it != text_last; ++text_it)
    {
      auto k_mer_begin {text_it};
      auto k_mer_end {std::next(k_mer_begin, k)};
      insert(text_begin, k_mer_begin, k_mer_end);
      alphabet.insert(*text_it);
    }
    for (auto text_it {text_last}; text_it != std::end(text); ++text_it)
    {
      alphabet.insert(*text_it);
    }
    entropy = calculate_cumulative_0th_entropy(root) / std::size(text);
  }

  template <typename text_iteraotr_type>
  void insert
  (
    text_iteraotr_type text_begin,
    text_iteraotr_type k_mer_it,
    text_iteraotr_type k_mer_end
  )
  {
    auto current_node {root};
    while (k_mer_it != k_mer_end)
    {
      auto character {*k_mer_it};
      auto branches_it {current_node->branches.find(character)};
      if (branches_it == std::end(current_node->branches))
      {
        current_node = current_node->branches[character] = std::make_shared<node>
        (
          std::distance(text_begin, k_mer_it),
          std::distance(text_begin, k_mer_end)
        );
        k_mer_it = k_mer_end;
      }
      else
      {
        auto child {std::get<1>(*branches_it)};
        auto edge_it {std::next(text_begin, child->edge_begin_offset)};
        auto edge_end {std::next(text_begin, child->edge_end_offset)};
        while
        (
          (edge_it != edge_end)
          &&
          (k_mer_it != k_mer_end)
          &&
          (*k_mer_it == *edge_it)
        )
        {
          ++k_mer_it;
          ++edge_it;
        }
        if
        (
          (edge_it != edge_end)
          &&
          (k_mer_it != k_mer_end)
        )
        {
          auto edge_it_offset {std::distance(text_begin, edge_it)};
          auto internal_node {std::make_shared<node>(child->edge_begin_offset, edge_it_offset)};
          child->edge_begin_offset = edge_it_offset;
          internal_node->branches[*edge_it] = child;
          current_node = internal_node->branches[*k_mer_it] = std::make_shared<node>
          (
            std::distance(text_begin, k_mer_it),
            std::distance(text_begin, k_mer_end)
          );
          std::get<1>(*branches_it) = internal_node;
        }
        else
        {
          current_node = child;
        }
      }
    }
    ++(current_node->context_size);
    auto next_character {*k_mer_end};
    auto alphabet_count_it {current_node->alphabet_count.find(next_character)};
    if (alphabet_count_it == std::end(current_node->alphabet_count))
    {
      current_node->alphabet_count[next_character] = 1;
    }
    else
    {
      ++std::get<1>(*alphabet_count_it);
    }
    return;
  }

  double calculate_cumulative_0th_entropy (std::shared_ptr<node> current_node)
  {
    double cumulative_0th_entropy {};
    if (current_node->context_size == 0)
    {
      auto branches_it {std::begin(current_node->branches)};
      auto branches_end {std::end(current_node->branches)};
      while (branches_it != branches_end)
      {
        cumulative_0th_entropy += calculate_cumulative_0th_entropy(std::get<1>(*branches_it));
        ++branches_it;
      }
    }
    else
    {
      auto log_context_size {std::log(static_cast<double>(current_node->context_size))};
      for (auto const &character_count : current_node->alphabet_count)
      {
        auto count {static_cast<double>(std::get<1>(character_count))};
        cumulative_0th_entropy += count * (log_context_size - std::log(count));
      }
    }
    return cumulative_0th_entropy;
  }
};

template <typename text_type>
void print_kth_entropy
(
  std::ofstream &fout,
  text_type &text,
  uint64_t k
)
{
  k_mer_trie trie(text, k);
  fout << k << "th_entropy(H_" << k << ")," << trie.entropy << "\n";
  fout << "compression_ratio," << (trie.entropy / std::size(trie.alphabet)) << "\n";
  std::cout << k << "th_entropy(H_" << k << ")," << trie.entropy << "\n";
  std::cout << "compression_ratio," << (trie.entropy / std::size(trie.alphabet)) << "\n";
  return;
}

template <typename text_type>
void print_p7zip_size
(
  std::ofstream &fout,
  text_type &text
)
{
  return;
}

template <typename text_type>
void print_bzip2_size
(
  std::ofstream &fout,
  text_type &text
)
{
  return;
}

template <typename text_type>
void print_repair_size
(
  std::ofstream &fout,
  text_type &text
)
{
  return;
}

void print_text_attributes (fs::path path)
{
  std::ifstream input_file {path};
  std::ofstream output_file {fs::path{"../output/attribute"} / (path.filename().string() + ".csv")};
  sdsl::int_vector<> text;
  text.load(input_file);
  output_file << "text_size(n)," << std::size(text) << "\n";
  std::cout << "text_size(n)," << std::size(text) << "\n";
  print_number_runs(output_file, text);
  print_alphabet_size(output_file, text);
  print_0th_entropy(output_file, text);
  for (uint64_t power {0}; power != 4; ++power)
  {
    print_kth_entropy(output_file, text, (1 << power));
  }
  print_p7zip_size(output_file, text);
  print_bzip2_size(output_file, text);
  print_repair_size(output_file, text);
  return;
}
}

#endif
