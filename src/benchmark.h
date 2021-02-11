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
  fs::path text_path
)
{
  auto index_path {fs::path{"../input/index"} / text_path.parent_path().filename()};
  if (!fs::exists(index_path))
  {
    fs::create_directories(index_path);
  }
  index_path /= text_path.filename();
  auto gc_index_path {index_path / (text_path.filename().string() + ".gci")};
  if (!fs::exists(gc_index_path))
  {
    construct(index, text_path);
    serialize(index, gc_index_path);
  }
  load(index, gc_index_path);
  auto fm_index_path {gc_index_path.replace_extension(".fmi")};
  if (!fs::exists(fm_index_path))
  {
    {
#ifdef LOG_FMI_CONSTRUCT
      auto json_path {fs::path("../data/fmi/construct") / text_path.parent_path()};
      if (!fs::exists(json_path))
      {
        fs::create_directories(json_path);
      }
      std::ofstream json_file {json_path / (text_path.filename().string() + ".json")};
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
          sdsl::load_vector_from_file(text, text_path.string());
          sdsl::construct_im(fm_index, text);
        }
#ifdef LOG_FMI_CONSTRUCT
      );
      phases.to_json().str(json_file);
#endif
    }
    {
      std::ofstream fm_index_file {fm_index_path};
#ifdef LOG_FMI_SERIALIZE
      auto json_path {fs::path("../data/fmi/serialize") / text_path.parent_path()};
      if (!fs::exists(json_path))
      {
        fs::create_directories(json_path);
      }
      std::ofstream json_file {json_path / (text_path.filename().string() + ".json")};
      tdc::StatPhase phases {"serialize_fm_index"};
      tdc::StatPhase::wrap
      (
        "",
        [&] ()
#endif
        {
          sdsl::serialize(fm_index, fm_index_file);
        }
#ifdef LOG_FMI_SERIALIZE
      );
      phases.to_json().str(json_file);
#endif
    }
  }
  {
    std::ifstream fm_index_file {fm_index_path};
#ifdef LOG_FMI_LOAD
    auto json_path {fs::path("../data/fmi/load") / text_path.parent_path()};
    if (!fs::exists(json_path))
    {
      fs::create_directories(json_path);
    }
    std::ofstream json_file {json_path / (text_path.filename().string() + ".json")};
    tdc::StatPhase phases {"load_fm_index"};
    tdc::StatPhase::wrap
    (
      "",
      [&] ()
#endif
      {
        sdsl::load(fm_index, fm_index_file);
      }
#ifdef LOG_FMI_LOAD
    );
    phases.to_json().str(json_file);
#endif
  }
  return;
}

void benchmark_gc_index_count
(
  gc_index &index,
  fs::path pattern_path
)
{
  std::ifstream pattern_file {pattern_path};
  sdsl::int_vector<> pattern;
  uint64_t pattern_amount {0};
  pattern_file >> pattern_amount;
#if defined(LOG_GCI_COUNT) || defined(LOG_VERBOSE_GCI_COUNT)
  gci::timer::register_record(gci::record_key::GCI_COUNT, "gc_index_count");
#endif
#if defined(LOG_VERBOSE_GCI_COUNT)
  gci::timer::register_record(gci::record_key::CALCULATE_SL_FACTOR, "calculate_sl_factor");
  gci::timer::register_record(gci::record_key::WAVELET_TREE_OPERATIONS_L, "wavelet_tree_operations_L");
  gci::timer::register_record(gci::record_key::LOOKUP_GRAMMAR_RULE, "lookup_grammar_rule");
  gci::timer::register_record(gci::record_key::WAVELET_TREE_OPERATIONS_S, "wavelet_tree_operations_S");
#endif
  for (uint64_t i {0}; i != pattern_amount; ++i)
  {
    pattern.load(pattern_file);
    gci::count(index, std::begin(pattern), std::end(pattern));
  }
#if defined(LOG_GCI_COUNT)
  auto csv_path {fs::path("../data/gci/count")};
#elif defined(LOG_VERBOSE_GCI_COUNT)
  auto csv_path {fs::path("../data/gci/count/verbose")};
#endif
#if defined(LOG_GCI_COUNT) || defined(LOG_VERBOSE_GCI_COUNT)
  if (!fs::exists(csv_path))
  {
    fs::create_directories(csv_path);
  }
  gci::timer::print(csv_path / (pattern_path.filename().string() + ".csv"));
#endif
  return;
}

void benchmark_fm_index_count
(
  sdsl::csa_wt<> &fm_index,
  fs::path pattern_path
)
{
  std::ifstream pattern_file {pattern_path};
  sdsl::int_vector<> pattern;
  uint64_t pattern_amount {0};
  pattern_file >> pattern_amount;
#ifdef LOG_FMI_COUNT
  gci::timer::register_record(gci::record_key::FMI_COUNT, "fm_index_count");
#endif
  for (uint64_t i {0}; i != pattern_amount; ++i)
  {
    pattern.load(pattern_file);
#ifdef LOG_FMI_COUNT
    gci::timer::resume(gci::record_key::FMI_COUNT);
#endif
    sdsl::count(fm_index, std::begin(pattern), std::end(pattern));
#ifdef LOG_FMI_COUNT
  gci::timer::pause(gci::record_key::FMI_COUNT);
#endif
  }
#ifdef LOG_FMI_COUNT
  gci::timer::print
  (
    "../output/count/" + gci::basename(pattern_path) + ".fmi.csv",
    "milliseconds"
  );
#endif
  return;
}

void benchmark_count
(
  fs::path text_path,
  uint64_t const pattern_amount,
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
    + gci::basename(text_path) + "."
    + std::to_string(pattern_amount) + "."
    + std::to_string(pattern_size)
  };
  gci::generate_pattern
  (
    text_path,
    pattern_path,
    pattern_amount,
    pattern_size
  );
  gc_index index;
  sdsl::csa_wt<> fm_index;
  load_gc_and_fm_index(index, fm_index, text_path);
  benchmark_gc_index_count(index, pattern_path);
  benchmark_fm_index_count(fm_index, pattern_path);
  sdsl::remove(pattern_path);
  return;
}

void test_count (fs::path text_path)
{
  gc_index index;
  sdsl::csa_wt<> fm_index;
  load_gc_and_fm_index(index, fm_index, text_path);
  sdsl::int_vector<> text;
  sdsl::load_vector_from_file(text, text_path);
  uint64_t max_sl_factor_size {calculate_max_sl_factor_size(text)};
  fs::path pattern_path {"../input/pattern/pattern.sample"};
  for (uint64_t multiple {2}; multiple != 11; ++multiple)
  {
    uint64_t pattern_amount {1000 / multiple};
    gci::generate_pattern
    (
      text_path,
      pattern_path,
      pattern_amount,
      static_cast<uint64_t>(max_sl_factor_size * multiple * 0.5)
    );
    std::ifstream pattern_file {pattern_path};
    sdsl::int_vector<8> pattern;
    pattern_file.read((char*)(&pattern_amount), sizeof(pattern_amount));
    for (uint64_t pattern_id {0}; pattern_id != pattern_amount; ++pattern_id)
    {
      pattern.load(pattern_file);
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
      // uint64_t pattern_amount {10 * divisor};
      uint64_t pattern_amount {1};
      gci::generate_pattern
      (
        text_path,
        pattern_path,
        pattern_amount,
        (std::size(text) / divisor)
      );
      std::ifstream pattern_file {pattern_path};
      sdsl::int_vector<8> pattern;
      pattern_file.read((char*)(&pattern_amount), sizeof(pattern_amount));
      for (uint64_t pattern_id {0}; pattern_id != pattern_amount; ++pattern_id)
      {
        pattern.load(pattern_file);
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
