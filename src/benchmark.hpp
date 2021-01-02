#ifndef BENCHMARK_HPP_
#define BENCHMARK_HPP_

#include <fstream>
#include <tudocomp_stat/StatPhase.hpp>

#include "grammar_compressed_index.hpp"
#include "pattern.hpp"
#include "utility.hpp"

namespace gci
{
void benchmark_gc_index_count
(
  gc_index &index,
  std::string const text_path,
  uint32_t const number_patterns,
  uint32_t const pattern_size
)
{
  auto min_pattern_size {gci::calculate_max_sl_factor_size(text_path)};
  if (pattern_size < min_pattern_size)
  {
    throw std::runtime_error
    (
      std::string{"pattern size must be at least "}
      + std::to_string(min_pattern_size)
      + " (characters) for this text"
    );
  }
  gci::generate_patterns
  (
    text_path,
    "../input/patterns",
    number_patterns,
    pattern_size
  );
  std::fstream fin ("../input/patterns", std::ios_base::in | std::ios_base::binary);
  sdsl::int_vector<8> pattern;
  std::ofstream fout (std::string{"../output/count/"} + util::basename(text_path) + ".json");
  tdc::StatPhase phases {""};
  tdc::StatPhase::wrap
  (
    "gc_index count",
    [&] ()
    {
      for (uint32_t i {0}; i != number_patterns; ++i)
      {
        tdc::StatPhase::pause_tracking();
        pattern.load(fin);
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
    }
  );
  phases.to_json().str(fout);
  return;
}

}

#endif
