#pragma once

#include <chrono>

#include <faster-minuter/wt_fbb.hpp>
#include <rlcsa/rlcsa.h>
#include <sdsl/csa_wt.hpp>
#include <sdsl/suffix_array_algorithm.hpp>
#include <tudocomp_stat/StatPhase.hpp>

#include "index.h"
#include "pattern_collection.h"
#include "utility.h"

namespace figiss
{
std::string ProperTimeRepresentation (double const nanoseconds)
{
  std::stringstream stringstream;
  stringstream << std::fixed << std::setprecision(2);
  if (nanoseconds < 1'000.0)
  {
    stringstream << nanoseconds;
    return (stringstream.str() + "ns");
  }
  else if (nanoseconds < 1'000'000.0)
  {
    stringstream << (nanoseconds / 1'000.0);
    return (stringstream.str() + "us");
  }
  else if (nanoseconds < 1'000'000'000.0)
  {
    stringstream << (nanoseconds / 1'000'000.0);
    return (stringstream.str() + "ms");
  }
  else
  {
    stringstream << (nanoseconds / 1'000'000'000.0);
    return (stringstream.str() + "s");
  }
  return "0.00";
}

void PrintIndexSpace
(
  std::filesystem::path const& index_path,
  std::filesystem::path const& index_space_path,
  bool const is_proper_representation = false
)
{
  if (!std::filesystem::exists(index_space_path.parent_path()))
  {
    std::filesystem::create_directories(index_space_path.parent_path());
  }
  std::ofstream out {index_space_path};
  std::cout << "write index space to " << std::filesystem::canonical(index_space_path).string() << "\n";
  {
    if (is_proper_representation)
    {
      out << "index space:" << ProperSizeRepresentation(std::filesystem::file_size(index_path)) << "B\n";
    }
    else
    {
      out << "index space(B):" << std::filesystem::file_size(index_path);
    }
  }
  return;
}

void PrintDetailedIndexSpace
(
  std::filesystem::path const& index_path,
  std::filesystem::path const& detailed_index_space_path,
  bool const is_proper_representation = false
)
{
  if (!std::filesystem::exists(detailed_index_space_path.parent_path()))
  {
    std::filesystem::create_directories(detailed_index_space_path.parent_path());
  }
  std::ofstream out {detailed_index_space_path};
  {
    Index index;
    index.Load(index_path);
    auto root {std::make_shared<SpaceNode>("index")};
    {
      auto dummy_index_path {std::filesystem::path{"_.figiss"}};
      index.Serialize(dummy_index_path, root);
      std::filesystem::remove(dummy_index_path);
    }
    std::cout << "write detailed index space to " << std::filesystem::canonical(detailed_index_space_path).string() << "\n";
    out << std::make_pair(root, is_proper_representation);
  }
  return;
}

void PrintConstructionTimeAndPeakMemory
(
  std::filesystem::path const& construction_path,
  std::string const& json,
  bool is_proper_representation = false
)
{
  if (!std::filesystem::exists(construction_path.parent_path()))
  {
    std::filesystem::create_directories(construction_path.parent_path());
  }
  auto begin_time {std::stoull(json.substr(json.find_last_of(":", json.find_last_of(":") - 1) + 1))};
  auto end_time {std::stoull(json.substr(json.find_last_of(":", json.find_last_of(":", json.find_last_of(":") - 1) - 1) + 1))};
  auto duration {end_time - begin_time + 1};
  auto peak_memory {std::stoull(json.substr(json.find(":", json.find(":", json.find(":") + 1) + 1) + 1))};
  std::cout
  << "time:" << ProperTimeRepresentation(duration * 1'000'000) << "\n"
  << "peak memory:" << ProperSizeRepresentation(peak_memory) << "B\n";
  std::ofstream out {construction_path};
  std::cout << "write construction time and peak memory to " << std::filesystem::canonical(construction_path).string() << "\n";
  if (is_proper_representation)
  {
    out
    << "time:" << ProperTimeRepresentation(duration * 1'000'000) << "\n"
    << "peak memory:" << ProperSizeRepresentation(peak_memory) << "B\n";
  }
  else
  {
    out
    << "time(ns):" << (duration * 1'000'000)<< "\n"
    << "peak memory(B):" << peak_memory << "\n";
  }
}

template <typename Index>
void ConstructAndSerialize
(
  std::filesystem::path const& byte_text_path,
  uint8_t const max_factor_size,
  std::string const& index_name
)
{
  Index index;
  std::cout << "--- construct & serialize " << index_name << " (\u03BB = " << std::to_string(max_factor_size) << ") ---\n";
  auto root {tdc::StatPhase("construction")};
  tdc::StatPhase::wrap
  (
    index_name.c_str(),
    [&] ()
    {
      index = Index {byte_text_path, max_factor_size};
    }
  );
  auto middle_path {std::filesystem::path{byte_text_path.filename().string() + "/" + index_name + "/" + std::to_string(max_factor_size)}};
  {
    auto construction_path {std::filesystem::path{"../data/construction"} / middle_path / std::filesystem::path{"construction"}};
    PrintConstructionTimeAndPeakMemory(construction_path, root.to_json().str());
  }
  {
    auto parameter_path {std::filesystem::path{"../data/parameter"} / middle_path / std::filesystem::path{"parameter"}};
    index.PrintParameters(parameter_path);
  }
  {
    auto index_path {std::filesystem::path{"../data/index"} / middle_path / std::filesystem::path{"index"}};
    index.Serialize(index_path);
    {
      auto index_space_path {std::filesystem::path{"../data/index_space"} / middle_path / std::filesystem::path{"index_space"}};
      PrintIndexSpace(index_path, index_space_path);
    }
    {
      auto detailed_index_space_path {std::filesystem::path{"../data/index_space"} / middle_path / std::filesystem::path{"detailed_index_space"}};
      PrintDetailedIndexSpace(index_path, detailed_index_space_path);
    }
  }
  return;
}

template <typename Index>
void ConstructAndSerialize
(
  std::filesystem::path const& byte_text_path,
  std::string const& index_name
);

template <>
void ConstructAndSerialize
<sdsl::csa_wt<wt_fbb<>, 0xFFFF'FFFF, 0xFFFF'FFFF>>
(
  std::filesystem::path const& byte_text_path,
  std::string const& index_name
)
{
  sdsl::csa_wt<wt_fbb<>, 0xFFFF'FFFF, 0xFFFF'FFFF> index;
  std::cout << "--- construct & serialize " << index_name << " ---\n";
  std::cout << "construct " << index_name << " of " << std::filesystem::canonical(byte_text_path).string() << "\n";
  auto root {tdc::StatPhase("construction")};
  tdc::StatPhase::wrap
  (
    index_name.c_str(),
    [&] ()
    {
      sdsl::int_vector<8> byte_text;
      sdsl::load_vector_from_file(byte_text, byte_text_path);
      sdsl::construct_im(index, byte_text);
    }
  );
  auto middle_path {byte_text_path.filename() / std::filesystem::path{index_name}};
  {
    auto construction_path {std::filesystem::path{"../data/construction"} / middle_path / std::filesystem::path{"construction"}};
    PrintConstructionTimeAndPeakMemory(construction_path, root.to_json().str());
  }
  {
    auto index_path {std::filesystem::path{"../data/index"} / middle_path / std::filesystem::path{"index"}};
    if (!std::filesystem::exists(index_path.parent_path()))
    {
      std::filesystem::create_directories(index_path.parent_path());
    }
    {
      std::ofstream out {index_path};
      std::cout << "serialize " + index_name + " to " << std::filesystem::canonical(index_path).string() << "\n";
      index.serialize(out);
    }
    {
      auto index_space_path {std::filesystem::path{"../data/index_space"} / middle_path / std::filesystem::path{"index_space"}};
      PrintIndexSpace(index_path, index_space_path);
    }
  }
  return;
}

template <>
void ConstructAndSerialize
<sdsl::csa_wt<sdsl::wt_rlmn<>, 0xFFFF'FFFF, 0xFFFF'FFFF>>
(
  std::filesystem::path const& byte_text_path,
  std::string const& index_name
)
{
  sdsl::csa_wt<sdsl::wt_rlmn<>, 0xFFFF'FFFF, 0xFFFF'FFFF> index;
  std::cout << "--- construct & serialize " << index_name << " ---\n";
  std::cout << "construct " << index_name << " of " << std::filesystem::canonical(byte_text_path).string() << "\n";
  auto root {tdc::StatPhase("construction")};
  tdc::StatPhase::wrap
  (
    index_name.c_str(),
    [&] ()
    {
      sdsl::int_vector<8> byte_text;
      sdsl::load_vector_from_file(byte_text, byte_text_path);
      sdsl::construct_im(index, byte_text);
    }
  );
  auto middle_path {byte_text_path.filename() / std::filesystem::path{index_name}};
  {
    auto construction_path {std::filesystem::path{"../data/construction"} / middle_path / std::filesystem::path{"construction"}};
    PrintConstructionTimeAndPeakMemory(construction_path, root.to_json().str());
  }
  {
    auto index_path {std::filesystem::path{"../data/index"} / middle_path / std::filesystem::path{"index"}};
    if (!std::filesystem::exists(index_path.parent_path()))
    {
      std::filesystem::create_directories(index_path.parent_path());
    }
    {
      std::ofstream out {index_path};
      std::cout << "serialize " + index_name + " to " << std::filesystem::canonical(index_path).string() << "\n";
      index.serialize(out);
    }
    {
      auto index_space_path {std::filesystem::path{"../data/index_space"} / middle_path / std::filesystem::path{"index_space"}};
      PrintIndexSpace(index_path, index_space_path);
    }
  }
  return;
}

template <>
void ConstructAndSerialize
<CSA::RLCSA>
(
  std::filesystem::path const& byte_text_path,
  std::string const& index_name
)
{
  std::cout << "--- construct & serialize " << index_name << " ---\n";
  std::cout << "construct " << index_name << " of " << std::filesystem::canonical(byte_text_path).string() << "\n";
  auto root {tdc::StatPhase("construction")};
  tdc::StatPhase::wrap
  (
    index_name.c_str(),
    [&] ()
    {
      sdsl::int_vector<8> byte_text;
      sdsl::load_vector_from_file(byte_text, byte_text_path);
      sdsl::append_zero_symbol(byte_text);
      CSA::RLCSA index((CSA::uchar*)byte_text.data(), std::size(byte_text), 32, 0, 1, false);
    }
  );
  auto middle_path {byte_text_path.filename() / std::filesystem::path{index_name}};
  {
    auto construction_path {std::filesystem::path{"../data/construction"} / middle_path / std::filesystem::path{"construction"}};
    PrintConstructionTimeAndPeakMemory(construction_path, root.to_json().str());
  }
  {
    auto index_base_path {std::filesystem::path{"../data/index"} / middle_path / byte_text_path.filename()};
    auto index_path {std::filesystem::path{index_base_path.string() + ".rlcsa.array"}};
    if (!std::filesystem::exists(index_path.parent_path()))
    {
      std::filesystem::create_directories(index_path.parent_path());
    }
    {
      sdsl::int_vector<8> byte_text;
      sdsl::load_vector_from_file(byte_text, byte_text_path);
      sdsl::append_zero_symbol(byte_text);
      CSA::RLCSA index((CSA::uchar*)byte_text.data(), std::size(byte_text), 32, 0, 1, false);
      index.writeTo(index_base_path.string());
    }
    {
      auto index_space_path {std::filesystem::path{"../data/index_space"} / middle_path / std::filesystem::path{"index_space"}};
      PrintIndexSpace(index_path, index_space_path);
    }
  }
  return;
}

template <typename Index, typename Iterator>
uint64_t Count (Index& index, Iterator begin, Iterator end)
{
  return index.Count(begin, end);
}

template <typename Iterator>
uint64_t Count (sdsl::csa_wt<sdsl::wt_rlmn<>, 0xFFFF'FFFF, 0xFFFF'FFFF>& index, Iterator begin, Iterator end)
{
  return sdsl::count(index, begin, end);
}

template <typename Iterator>
uint64_t Count (sdsl::csa_wt<wt_fbb<>, 0xFFFF'FFFF, 0xFFFF'FFFF>& index, Iterator begin, Iterator end)
{
  return sdsl::count(index, begin, end);
}

template <typename Iterator>
uint64_t Count (CSA::RLCSA& index, Iterator begin, Iterator end)
{
  auto it {std::prev(end)};
  auto rend {std::prev(begin)};
  auto pattern_range {index.getCharRange(*it)};
  while (!CSA::isEmpty(pattern_range) && --it != rend)
  {
    pattern_range = index.LF(pattern_range, *it);
  }
  if (!CSA::isEmpty(pattern_range))
  {
    index.convertToSARange(pattern_range);
    return CSA::length(pattern_range);
  }
  return 0;
}

template <typename Index>
void MeasureCountingTime
(
  std::filesystem::path const& pattern_path,
  Index& index,
  std::filesystem::path const& time_path,
  bool const is_proper_representation = false
)
{
  if (!std::filesystem::exists(pattern_path))
  {
    std::cerr << "\033[31mfailed at finding pattern path\033[0m\n";
    std::cerr << std::filesystem::canonical(pattern_path);
    return;
  }
  PatternCollection pattern_collection {};
  {
    pattern_collection.Load(pattern_path);
    std::cout << "warm up cache by counting on pattern collection once\n";
    auto begin {std::begin(pattern_collection)};
    while (begin != std::end(pattern_collection))
    {
      auto end {std::next(begin, pattern_collection.GetLength())};
      Count(index, begin, end);
      begin = end;
    }
  }
  {
    auto begin {std::begin(pattern_collection)};
    std::cout << "count on pattern collection once again and measure average counting time\n";
    auto begin_time {std::chrono::steady_clock::now()};
    while (begin != std::end(pattern_collection))
    {
      auto end {std::next(begin, pattern_collection.GetLength())};
      Count(index, begin, end);
      begin = end;
    }
    auto duration {static_cast<double>((std::chrono::steady_clock::now() - begin_time).count())};
    duration /= (pattern_collection.GetAmount() * pattern_collection.GetLength());
    {
      if (!std::filesystem::exists(time_path.parent_path()))
      {
        std::filesystem::create_directories(time_path.parent_path());
      }
      std::cout << "average counting time per character:";
      if (is_proper_representation)
      {
        std::cout << ProperTimeRepresentation(duration) << "\n";
      }
      else
      {
        std::cout << std::fixed << std::setprecision(2);
        std::cout << duration << "ns\n";
      }
      std::ofstream out {time_path};
      std::cout << "write counting time to " << std::filesystem::canonical(time_path) << "\n";
      if (is_proper_representation)
      {
        out << "average counting time per character:" << ProperTimeRepresentation(duration) << "\n";
      }
      else
      {
        out << "average counting time per character(ns):" << duration << "\n";
      }
    }
  }
  return;
}
}
