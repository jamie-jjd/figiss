#ifndef BENCHMARK_HPP_
#define BENCHMARK_HPP_

#include <tudocomp_stat/StatPhase.hpp>

#include "grammar_compressed_index.hpp"
#include "pattern.hpp"
#include "utility.hpp"

namespace gci
{
void load_gc_and_fm_index
(
  gc_index &index,
  sdsl::csa_wt<> &fm_index,
  std::string const text_path
)
{
  std::string gc_index_path {"../input/index/" + util::basename(text_path) + ".gci"};
  std::ifstream gc_index_input {gc_index_path};
  if (!gc_index_input.is_open())
  {
    construct(index, text_path);
    serialize(index, gc_index_path);
  }
  else
  {
    load(index, gc_index_path);
  }
  std::string fm_index_path {"../input/index/" + util::basename(text_path) + ".fmi"};
  std::ifstream fm_index_input {fm_index_path};
  if (!fm_index_input.is_open())
  {
    sdsl::int_vector<8> text;
    sdsl::load_vector_from_file(text, text_path);
    sdsl::construct_im(fm_index, text);
    std::ofstream fm_index_output {fm_index_path};
    sdsl::serialize(fm_index, fm_index_output);
  }
  else
  {
    sdsl::load(fm_index, fm_index_input);
  }
  return;
}

void benchmark_gc_index_count
(
  gc_index &index,
  std::string const pattern_path
)
{
  std::ifstream pattern_input {pattern_path};
  std::ofstream json_output {"../output/count/" + util::basename(pattern_path) + ".gci.json"};
  sdsl::int_vector<8> pattern;
  uint64_t pattern_number {0};
  pattern_input.read((char*)(&pattern_number), sizeof(pattern_number));
  tdc::StatPhase phases {"gc_index_count"};
  for (uint64_t i {0}; i != pattern_number; ++i)
  {
    tdc::StatPhase::pause_tracking();
    pattern.load(pattern_input);
    tdc::StatPhase::resume_tracking();
    tdc::StatPhase::wrap
    (
      "",
      [&] ()
      {
        gci::count(index, std::begin(pattern), std::end(pattern));
      }
    );
  }
  phases.to_json().str(json_output);
  return;
}

void benchmark_fm_index_count
(
  sdsl::csa_wt<> &fm_index,
  std::string const pattern_path
)
{
  std::ifstream pattern_input {pattern_path};
  std::ofstream json_output {"../output/count/" + util::basename(pattern_path) + ".fmi.json"};
  sdsl::int_vector<8> pattern;
  uint64_t pattern_number {0};
  pattern_input.read((char*)(&pattern_number), sizeof(pattern_number));
  tdc::StatPhase phases {"fm_index_count"};
  for (uint64_t i {0}; i != pattern_number; ++i)
  {
    tdc::StatPhase::pause_tracking();
    pattern.load(pattern_input);
    tdc::StatPhase::resume_tracking();
    tdc::StatPhase::wrap
    (
      "",
      [&] ()
      {
        sdsl::count(fm_index, std::begin(pattern), std::end(pattern));
      }
    );
  }
  phases.to_json().str(json_output);
  return;
}

void benchmark_count
(
  std::string const text_path,
  uint64_t const pattern_number,
  uint64_t const pattern_size
)
{
  auto min_pattern_size {calculate_max_sl_factor_size(text_path)};
  if (pattern_size < min_pattern_size)
  {
    throw std::runtime_error
    (
      std::string{"pattern size should be at least "}
      + std::to_string(min_pattern_size)
      + " (characters) for this text"
    );
  }
  std::string pattern_path
  {
    "../input/pattern/"
    + util::basename(text_path) + "."
    + std::to_string(pattern_number) + "."
    + std::to_string(pattern_size)
  };
  generate_pattern
  (
    text_path,
    pattern_path,
    pattern_number,
    pattern_size
  );
  gc_index index;
  sdsl::csa_wt<> fm_index;
  load_gc_and_fm_index(index, fm_index, text_path);
  benchmark_gc_index_count(index, pattern_path);
  benchmark_fm_index_count(fm_index, pattern_path);
  return;
}

void test_count (std::string const text_path)
{
  gc_index index;
  sdsl::csa_wt<> fm_index;
  load_gc_and_fm_index(index, fm_index, text_path);
  sdsl::int_vector<8> text;
  sdsl::load_vector_from_file(text, text_path);
  uint64_t max_sl_factor_size {calculate_max_sl_factor_size(text)};
  std::string const pattern_path {"pattern.sample"};
  for (uint64_t multiple {2}; multiple != 7; ++multiple)
  {
    uint64_t pattern_number {100};
    generate_pattern
    (
      text_path,
      pattern_path,
      pattern_number,
      static_cast<uint64_t>(max_sl_factor_size * multiple * 0.5)
    );
    std::ifstream pattern_input {pattern_path};
    sdsl::int_vector<8> pattern;
    pattern_input.read((char*)(&pattern_number), sizeof(pattern_number));
    for (uint64_t pattern_id {0}; pattern_id != pattern_number; ++pattern_id)
    {
      pattern.load(pattern_input);
      auto fm_count {sdsl::count(fm_index, std::begin(pattern), std::end(pattern))};
      auto gc_count {gci::count(index, std::begin(pattern), std::end(pattern))};
      if (fm_count != gc_count)
      {
        throw std::runtime_error
        (
          "failed at pattern size: " +
          std::to_string(static_cast<uint64_t>(max_sl_factor_size * multiple * 0.5))
        );
      }
    }
    std::cout << "pass pattern size: " << std::to_string(static_cast<uint64_t>(max_sl_factor_size * multiple * 0.5)) << '\n';
  }
  for (uint64_t divisor {5}; divisor != 1; --divisor)
  {
    if ((std::size(text) / divisor) > max_sl_factor_size)
    {
      uint64_t pattern_number {1};
      generate_pattern
      (
        text_path,
        pattern_path,
        pattern_number,
        (std::size(text) / divisor)
      );
      std::ifstream pattern_input {pattern_path};
      sdsl::int_vector<8> pattern;
      pattern_input.read((char*)(&pattern_number), sizeof(pattern_number));
      for (uint64_t pattern_id {0}; pattern_id != pattern_number; ++pattern_id)
      {
        pattern.load(pattern_input);
        auto fm_count {sdsl::count(fm_index, std::begin(pattern), std::end(pattern))};
        auto gc_count {gci::count(index, std::begin(pattern), std::end(pattern))};
        if (fm_count != gc_count)
        {
          throw std::runtime_error
          (
            "failed at pattern size: " +
            std::to_string(std::size(text) / divisor)
          );
        }
      }
      std::cout << "pass pattern size: " + std::to_string(std::size(text) / divisor) << '\n';
    }
  }
  auto fm_count {sdsl::count(fm_index, std::begin(text), std::end(text))};
  auto gc_count {gci::count(index, std::begin(text), std::end(text))};
  if (fm_count != gc_count)
  {
    throw std::runtime_error("failed at full text");
  }
  std::cout << "pass full text\n";
  sdsl::remove(pattern_path);
  return;
}
}

#endif
