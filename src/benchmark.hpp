#ifndef BENCHMARK_HPP_
#define BENCHMARK_HPP_

#include <tudocomp_stat/StatPhase.hpp>

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
      std::ofstream json_output {"../output/construct/" + util::basename(text_path) + ".fmi.json"};
      tdc::StatPhase phases {"construct_fm_index"};
      sdsl::int_vector<8> text;
      tdc::StatPhase::wrap
      (
        "",
        [&] ()
        {
          sdsl::load_vector_from_file(text, text_path);
          sdsl::construct_im(fm_index, text);
        }
      );
      phases.to_json().str(json_output);
    }
    {
      std::ofstream fm_index_output {fm_index_path};
      std::ofstream json_output {"../output/serialize/" + util::basename(text_path) + ".fmi.json"};
      tdc::StatPhase phases {"serialize_fm_index"};
      tdc::StatPhase::wrap
      (
        "",
        [&] ()
        {
          sdsl::serialize(fm_index, fm_index_output);
        }
      );
      phases.to_json().str(json_output);
    }
  }
  else
  {
    std::ofstream json_output {"../output/load/" + util::basename(text_path) + ".fmi.json"};
    tdc::StatPhase phases {"load_fm_index"};
    tdc::StatPhase::wrap
    (
      "",
      [&] ()
      {
        sdsl::load(fm_index, fm_index_input);
      }
    );
    phases.to_json().str(json_output);
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
  std::string const time_path {"../output/count/" + util::basename(pattern_path) + ".wm.sub.gci.csv"};
  gci::util::timer timer;
  sdsl::int_vector<8> pattern;
  uint64_t pattern_number {0};
  pattern_input.read((char*)(&pattern_number), sizeof(pattern_number));
  timer.reset("gc_index_count");
  timer.resume("gc_index_count");
  for (uint64_t i {0}; i != pattern_number; ++i)
  {
    timer.pause("gc_index_count");
    pattern.load(pattern_input);
    timer.resume("gc_index_count");
    gci::count(index, std::begin(pattern), std::end(pattern), timer);
  }
  timer.pause("gc_index_count");
  timer.print(time_path, "microseconds");
  return;
}

void benchmark_fm_index_count
(
  sdsl::csa_wt<> &fm_index,
  std::string const pattern_path
)
{
  std::ifstream pattern_input {pattern_path};
  std::string const time_path {"../output/count/" + util::basename(pattern_path) + ".fmi.csv"};
  gci::util::timer timer;
  sdsl::int_vector<8> pattern;
  uint64_t pattern_number {0};
  pattern_input.read((char*)(&pattern_number), sizeof(pattern_number));
  timer.reset("fm_index_count");
  timer.resume("fm_index_count");
  for (uint64_t i {0}; i != pattern_number; ++i)
  {
    timer.pause("fm_index_count");
    pattern.load(pattern_input);
    timer.resume("fm_index_count");
    sdsl::count(fm_index, std::begin(pattern), std::end(pattern));
  }
  timer.pause("fm_index_count");
  timer.print(time_path, "microseconds");
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
  gci::util::timer timer;
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
      auto gc_count {gci::count(index, std::begin(pattern), std::end(pattern), timer)};
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
      uint64_t pattern_number {10 * divisor};
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
        auto gc_count {gci::count(index, std::begin(pattern), std::end(pattern), timer)};
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
  auto gc_count {gci::count(index, std::begin(text), std::end(text), timer)};
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
