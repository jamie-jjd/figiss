#ifndef BENCHMARK_HPP_
#define BENCHMARK_HPP_

#include <tudocomp_stat/StatPhase.hpp>

#include "grammar_compressed_index.hpp"
#include "pattern.hpp"
#include "utility.hpp"

namespace gci
{
void benchmark_gc_index_count
(
  gc_index &index,
  std::string const pattern_path
)
{
  std::ifstream patterns_input {pattern_path};
  std::ofstream json_output {"../output/count/" + util::basename(pattern_path) + "_gc_index.json"};
  sdsl::int_vector<8> pattern;
  uint32_t pattern_number {0};
  patterns_input >> pattern_number;
  tdc::StatPhase phases {"gc_index count"};
  for (uint32_t i {0}; i != pattern_number; ++i)
  {
    tdc::StatPhase::pause_tracking();
    pattern.load(patterns_input);
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
  std::ifstream patterns_input {pattern_path};
  std::ofstream json_output {"../output/count/" + util::basename(pattern_path) + "_fm_index.json"};
  sdsl::int_vector<8> pattern;
  uint32_t pattern_number {0};
  patterns_input >> pattern_number;
  tdc::StatPhase phases {"fm_index count"};
  for (uint32_t i {0}; i != pattern_number; ++i)
  {
    tdc::StatPhase::pause_tracking();
    pattern.load(patterns_input);
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
  gc_index &index,
  sdsl::csa_wt<> &fm_index,
  std::string const text_path,
  uint32_t const pattern_number,
  uint32_t const pattern_size
)
{
  auto min_pattern_size {gci::calculate_max_sl_factor_size(text_path)};
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
    + util::basename(text_path) + "_"
    + std::to_string(pattern_number) + "_"
    + std::to_string(pattern_size)
  };
  generate_patterns
  (
    text_path,
    pattern_path,
    pattern_number,
    pattern_size
  );
  benchmark_gc_index_count(index, pattern_path);
  benchmark_fm_index_count(fm_index, pattern_path);
  sdsl::remove(pattern_path);
  return;
}

}

#endif
