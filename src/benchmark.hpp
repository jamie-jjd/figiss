#ifndef BENCHMARK_HPP_
#define BENCHMARK_HPP_

#ifdef LOG
#include <tudocomp_stat/StatPhase.hpp>
#endif
#include "grammar_compressed_index.hpp"
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
    {
#ifdef LOG_FMI_CONSTRUCT
      std::ofstream json_output {"../output/construct/" + util::basename(text_path) + ".fmi.json"};
      tdc::StatPhase phases {"construct_fm_index"};
#endif
      sdsl::int_vector<8> text;
#ifdef LOG_FMI_CONSTRUCT
      tdc::StatPhase::wrap
      (
        "",
        [&] ()
#endif
        {
          sdsl::load_vector_from_file(text, text_path);
          sdsl::construct_im(fm_index, text);
        }
#ifdef LOG_FMI_CONSTRUCT
      );
      phases.to_json().str(json_output);
#endif
    }
    {
      std::ofstream fm_index_output {fm_index_path};
#ifdef LOG_FMI_SERIALIZE
      std::ofstream json_output {"../output/serialize/" + util::basename(text_path) + ".fmi.json"};
      tdc::StatPhase phases {"serialize_fm_index"};
      tdc::StatPhase::wrap
      (
        "",
        [&] ()
#endif
        {
          sdsl::serialize(fm_index, fm_index_output);
        }
#ifdef LOG_FMI_SERIALIZE
      );
      phases.to_json().str(json_output);
#endif
    }
  }
  else
  {
#ifdef LOG_FMI_LOAD
    std::ofstream json_output {"../output/load/" + util::basename(text_path) + ".fmi.json"};
    tdc::StatPhase phases {"load_fm_index"};
    tdc::StatPhase::wrap
    (
      "",
      [&] ()
#endif
      {
        sdsl::load(fm_index, fm_index_input);
      }
#ifdef LOG_FMI_LOAD
    );
    phases.to_json().str(json_output);
#endif
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
  sdsl::int_vector<8> pattern;
  uint64_t pattern_number {0};
  pattern_input.read((char*)(&pattern_number), sizeof(pattern_number));
#if defined(LOG_GCI_COUNT) || defined(LOG_VERBOSE_GCI_COUNT)
  gci::util::timer::register_category(gci::util::category_id::GCI_COUNT, "gc_index_count");
#endif
#if defined(LOG_VERBOSE_GCI_COUNT)
  gci::util::timer::register_category(gci::util::category_id::CALCULATE_SL_FACTOR, "calculate_sl_factor");
  gci::util::timer::register_category(gci::util::category_id::WAVELET_TREE_OPERATIONS_L, "wavelet_tree_operations_L");
  gci::util::timer::register_category(gci::util::category_id::LOOKUP_GRAMMAR_RULE, "lookup_grammar_rule");
  gci::util::timer::register_category(gci::util::category_id::WAVELET_TREE_OPERATIONS_S, "wavelet_tree_operations_S");
#endif
  for (uint64_t i {0}; i != pattern_number; ++i)
  {
    pattern.load(pattern_input);
    gci::count(index, std::begin(pattern), std::end(pattern));
  }
#if defined(LOG_GCI_COUNT)
  gci::util::timer::print
  (
    "../output/count/" + util::basename(pattern_path) + ".gci.csv",
    "milliseconds"
  );
#elif defined(LOG_VERBOSE_GCI_COUNT)
  gci::util::timer::print
  (
    "../output/count/" + util::basename(pattern_path) + ".gci.verbose.csv",
    "milliseconds"
  );
#endif
  return;
}

void benchmark_fm_index_count
(
  sdsl::csa_wt<> &fm_index,
  std::string const pattern_path
)
{
  std::ifstream pattern_input {pattern_path};
  sdsl::int_vector<8> pattern;
  uint64_t pattern_number {0};
  pattern_input.read((char*)(&pattern_number), sizeof(pattern_number));
#ifdef LOG_FMI_COUNT
  gci::util::timer::register_category(gci::util::category_id::FMI_COUNT, "fm_index_count");
#endif
  for (uint64_t i {0}; i != pattern_number; ++i)
  {
    pattern.load(pattern_input);
#ifdef LOG_FMI_COUNT
    gci::util::timer::resume(gci::util::category_id::FMI_COUNT);
#endif
    sdsl::count(fm_index, std::begin(pattern), std::end(pattern));
#ifdef LOG_FMI_COUNT
  gci::util::timer::pause(gci::util::category_id::FMI_COUNT);
#endif
  }
#ifdef LOG_FMI_COUNT
  gci::util::timer::print
  (
    "../output/count/" + util::basename(pattern_path) + ".fmi.csv",
    "milliseconds"
  );
#endif
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
  gci::util::generate_pattern
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
  std::string const pattern_path {"../input/pattern/pattern.sample"};
  for (uint64_t multiple {2}; multiple != 11; ++multiple)
  {
    uint64_t pattern_number {1000 / multiple};
    gci::util::generate_pattern
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
    std::cout << "succeed at pattern size: " << std::to_string(static_cast<uint64_t>(max_sl_factor_size * multiple * 0.5)) << '\n';
  }
  for (uint64_t divisor {5}; divisor != 1; --divisor)
  {
    if ((std::size(text) / divisor) >= max_sl_factor_size)
    {
      // uint64_t pattern_number {10 * divisor};
      uint64_t pattern_number {1};
      gci::util::generate_pattern
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
      std::cout << "succeed at pattern size: " + std::to_string(std::size(text) / divisor) << '\n';
    }
  }
  auto fm_count {sdsl::count(fm_index, std::begin(text), std::end(text))};
  auto gc_count {gci::count(index, std::begin(text), std::end(text))};
  if (fm_count != gc_count)
  {
    throw std::runtime_error("failed at text size");
  }
  std::cout << "succeed at text size\n";
  sdsl::remove(pattern_path);
  return;
}
}

#endif
