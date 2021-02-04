#ifndef UTILITY_HPP_
#define UTILITY_HPP_

#include <chrono>
#include <fstream>
#include <functional>
#include <random>
#include <string>
#include <map>
#include <set>

#include <cmath>

namespace gci
{
namespace util
{
void generate_pattern
(
  std::string const text_path,
  std::string const pattern_path,
  uint64_t const pattern_number,
  uint64_t const pattern_size
)
{
  sdsl::int_vector<8> text;
  sdsl::load_vector_from_file(text, text_path);
  if (pattern_size < std::size(text))
  {
    std::ofstream pattern_output {pattern_path};
#ifndef LOG_VERBOSE_GCI_COUNT
    std::mt19937 engine {std::random_device{}()};
    std::uniform_int_distribution<uint64_t> distribution(0, std::size(text) - pattern_size);
    auto random_begin_dist {std::bind(distribution, engine)};
#endif
    sdsl::int_vector<8> pattern(pattern_size);
    pattern_output.write((char*)(&pattern_number), sizeof(pattern_number));
    for (uint64_t i {0}; i != pattern_number; ++i)
    {
#ifndef LOG_VERBOSE_GCI_COUNT
      auto text_it {std::next(std::begin(text), random_begin_dist())};
#else
      auto text_it {std::next(std::begin(text), (i * pattern_size) % (std::size(text) - pattern_size))};
#endif
      auto pattern_it {std::begin(pattern)};
      auto pattern_end {std::end(pattern)};
      while (pattern_it != pattern_end)
      {
        *pattern_it = *text_it;
        ++text_it;
        ++pattern_it;
      }
      pattern.serialize(pattern_output);
    }
  }
  else
  {
    std::cout << "pattern size >= text size\n";
    std::cout << "generated pattern is one copy of the text\n";
    std::ofstream pattern_output {pattern_path};
    uint64_t pattern_number_ {1};
    pattern_output.write((char*)(&pattern_number_), sizeof(pattern_number_));
    text.serialize(pattern_output);
  }
  return;
}

std::string basename (std::string const path)
{
  auto basename_begin_dist {path.find_last_of("/")};
  if (basename_begin_dist != std::string::npos)
  {
    ++basename_begin_dist;
  }
  else
  {
    basename_begin_dist = 0;
  }
  return path.substr(basename_begin_dist);
}

enum class category_id
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
  using inner_clock = std::chrono::high_resolution_clock;

  static std::map
  <
    category_id,
    std::tuple
    <
      std::string,
      std::chrono::nanoseconds,
      inner_clock::time_point
    >
  >
  records;

  static void register_category
  (
    category_id id,
    std::string const category
  )
  {
    records[id] =
    {
      category,
      std::chrono::nanoseconds::zero(),
      inner_clock::time_point::min()
    };
    return;
  }

  static void resume (category_id id)
  {
    std::get<2>(records[id]) = inner_clock::now();
    return;
  }

  static void pause (category_id id)
  {
    std::get<1>(records[id]) +=
    std::chrono::duration_cast<std::chrono::nanoseconds>
    (
      inner_clock::now()
      -
      std::get<2>(records[id])
    );
    return;
  }

  static void print
  (
    std::string const path,
    std::string const unit = "seconds"
  )
  {
    std::ofstream fout {path};
    fout << "category," << unit << "\n";
    for (auto const &record : records)
    {
      auto category {std::get<0>(std::get<1>(record))};
      auto duration {std::get<1>(std::get<1>(record))};
      fout << category << ",";
      if (unit == "minutes")
      {
        fout << std::chrono::duration_cast<std::chrono::minutes>(duration).count() << "\n";
      }
      else if (unit == "seconds")
      {
        fout << std::chrono::duration_cast<std::chrono::seconds>(duration).count() << "\n";
      }
      else if (unit == "milliseconds")
      {
        fout << std::chrono::duration_cast<std::chrono::milliseconds>(duration).count() << "\n";
      }
      else if (unit == "microseconds")
      {
        fout << std::chrono::duration_cast<std::chrono::microseconds>(duration).count() << "\n";
      }
      else
      {
        fout << duration.count() << "\n";
      }
    }
    records.clear();
    return;
  }
};

std::map
<
  category_id,
  std::tuple
  <
    std::string,
    std::chrono::nanoseconds,
    timer::inner_clock::time_point
  >
>
timer::records
{};

void generate_canonical_text (std::string const path)
{
  sdsl::int_vector<8> text;
  sdsl::load_vector_from_file(text, path);
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
  std::ofstream output_file {"../input/text/sdsl_vector/" + basename(path) + ".sdsl"};
  sdsl::append_zero_symbol(text);
  sdsl::util::bit_compress(text);
  text.serialize(output_file);
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
  std::cout << "number_runs(r)," << number_runs << "\n";
  return;
}

template <typename text_type>
void print_alphabet_size
(
  std::ofstream &fout,
  text_type &text
)
{
  std::set<uint64_t> alphabet;
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
    if (alphabet_count.find(character) != alphabet_count.end())
    {
      ++alphabet_count[character];
    }
    else
    {
      alphabet_count[character] = 1;
    }
  }
  double entropy {0.0};
  auto log_text_size {std::log(std::size(text))};
  for (auto const &character_count : alphabet_count)
  {
    auto count {static_cast<double>(std::get<1>(character_count))};
    entropy += (count * (log_text_size - std::log(count)));
  }
  entropy /= std::size(text);
  fout << "0th_entropy(H_0)," << entropy << "\n";
  fout << "compression_ratio," << (entropy / std::size(alphabet_count)) << "\n";
  std::cout << "0th_entropy(H_0)," << entropy << "\n";
  std::cout << "compression_ratio," << (entropy / std::size(alphabet_count)) << "\n";
  return;
}

template <typename text_type>
void print_kth_entropy
(
  std::ofstream &fout,
  text_type &text,
  uint64_t k
)
{
  std::set<uint64_t> alphabet;
  std::map<text_type, std::map<uint64_t, uint64_t>> k_mers_alphabet_count;
  text_type k_mer(k);
  for (uint64_t i {0}; i != k; ++i)
  {
    k_mer[i] = text[i];
    alphabet.insert(text[i]);
  }
  auto character {text[k]};
  k_mers_alphabet_count[k_mer][character] = 1;
  alphabet.insert(character);
  for (uint64_t i {k + 1}; i != std::size(text); ++i)
  {
    for (uint64_t j {1}; j != k; ++j)
    {
      k_mer[j - 1] = k_mer[j];
    }
    k_mer[k - 1] = character;
    character = text[i];
    alphabet.insert(character);
    if (k_mers_alphabet_count.find(k_mer) != k_mers_alphabet_count.end())
    {
      auto &alphabet_count {k_mers_alphabet_count[k_mer]};
      if (alphabet_count.find(character) != alphabet_count.end())
      {
        ++alphabet_count[character];
      }
      else
      {
        alphabet_count[character] = 1;
      }
    }
    else
    {
      k_mers_alphabet_count[k_mer][character] = 1;
    }
  }
  double kth_entropy {0.0};
  for (auto const &k_mer_alphabet_count : k_mers_alphabet_count)
  {
    auto &alphabet_count {std::get<1>(k_mer_alphabet_count)};
    uint64_t k_mer_context_size {0};
    for (auto const &character_count: alphabet_count)
    {
      k_mer_context_size += std::get<1>(character_count);
    }
    auto log_k_mer_context_size {std::log(static_cast<double>(k_mer_context_size))};
    for (auto const &character_count: alphabet_count)
    {
      auto count {static_cast<double>(std::get<1>(character_count))};
      kth_entropy += count * (log_k_mer_context_size - std::log(count));
    }
  }
  kth_entropy /= std::size(text);
  fout << k << "th_entropy(H_" << k << ")," << kth_entropy << "\n";
  fout << "compression_ratio," << (kth_entropy / std::size(alphabet)) << "\n";
  std::cout << k << "th_entropy(H_" << k << ")," << kth_entropy << "\n";
  std::cout << "compression_ratio," << (kth_entropy / std::size(alphabet)) << "\n";
  return;
}

void print_text_attributes (std::string const path)
{
  std::ifstream input_file {path};
  std::ofstream output_file {"../output/attribute/" + basename(path) + ".csv"};
  sdsl::int_vector<> text;
  text.load(input_file);
  input_file.close();
  output_file << "text_size(n)," << std::size(text) << "\n";
  std::cout << "text_size(n)," << std::size(text) << "\n";
  print_number_runs(output_file, text);
  print_alphabet_size(output_file, text);
  print_0th_entropy(output_file, text);
  for (uint64_t power {0}; power != 4; ++power)
  {
    print_kth_entropy(output_file, text, (1 << power));
  }
  return;
}

}
}

#endif
