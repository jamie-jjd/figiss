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
  if (nanoseconds >= 1'000'000'000.0)
  {
    stringstream << (nanoseconds / 1'000'000'000.0);
    return (stringstream.str() + "s");
  }
  else if (nanoseconds >= 1'000'000.0)
  {
    stringstream << (nanoseconds / 1'000'000.0);
    return (stringstream.str() + "ms");
  }
  else if (nanoseconds >= 1'000.0)
  {
    stringstream << (nanoseconds / 1'000.0);
    return (stringstream.str() + "us");
  }
  else
  {
    stringstream << nanoseconds;
    return (stringstream.str() + "ns");
  }
  return "0.00";
}

void PrintSpace
(
  std::filesystem::path const& index_path,
  bool const is_proper_representation = false
)
{
  auto space_path {std::filesystem::path{"../data/space/" + index_path.filename().string() + ".space"}};
  if (!std::filesystem::exists(space_path.parent_path()))
  {
    std::filesystem::create_directories(space_path.parent_path());
  }
  std::ofstream out {space_path};
  std::cout << "write index space to " << std::filesystem::canonical(space_path).string() << "\n";
  {
    if (is_proper_representation)
    {
      out << "index space," << ProperSizeRepresentation(std::filesystem::file_size(index_path)) << "B\n";
    }
    else
    {
      out << "index space (B)," << std::filesystem::file_size(index_path);
    }
  }
  return;
}

void PrintDetailedSpace
(
  std::filesystem::path const& index_path,
  Index<>& index,
  bool const is_proper_representation = false
)
{
  auto detailed_space_path {std::filesystem::path{"../data/space/" + index_path.filename().string() + ".detailed.space"}};
  if (!std::filesystem::exists(detailed_space_path.parent_path()))
  {
    std::filesystem::create_directories(detailed_space_path.parent_path());
  }
  std::ofstream out {detailed_space_path};
  {
    auto root {std::make_shared<SpaceNode>("index")};
    index.Load(index_path);
    {
      auto dummy_index_path {std::filesystem::path{"_.figiss"}};
      index.Serialize(dummy_index_path, root);
      std::filesystem::remove(dummy_index_path);
    }
    std::cout << "write detailed space to " << std::filesystem::canonical(detailed_space_path).string() << "\n";
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
  auto begin_time {std::stoull(json.substr(json.find_last_of(":", json.find_last_of(":") - 1) + 1))};
  auto end_time {std::stoull(json.substr(json.find_last_of(":", json.find_last_of(":", json.find_last_of(":") - 1) - 1) + 1))};
  auto duration {end_time - begin_time + 1};
  auto peak_memory {std::stoull(json.substr(json.find(":", json.find(":", json.find(":") + 1) + 1) + 1))};
  std::cout
  << "construction time:" << ProperTimeRepresentation(duration) << "\n"
  << "peak memory:" << ProperSizeRepresentation(peak_memory) << "B\n";
  std::ofstream out {construction_path};
  std::cout << "write construction time and peak memory to " << std::filesystem::canonical(construction_path).string() << "\n";
  if (is_proper_representation)
  {
    out
    << "construction time," << ProperTimeRepresentation(duration) << "\n"
    << "peak memory," << ProperSizeRepresentation(peak_memory) << "B\n";
  }
  else
  {
    out
    << "construction time (ns)," << duration << "\n"
    << "peak memory (B)," << peak_memory << "\n";
  }
}

template <typename Index>
void ConstructAndSerialize
(
  std::filesystem::path const& byte_text_path,
  std::filesystem::path const& construction_parent_path,
  std::filesystem::path const& index_parent_path
)
{
  Index index;
  std::cout << "--- construct & serialize figiss ---\n";
  auto root {tdc::StatPhase("construction")};
  tdc::StatPhase::wrap
  (
    "figiss",
    [&] ()
    {
      index = Index{byte_text_path};
    }
  );
  {
    auto construction_path {construction_parent_path / std::filesystem::path{byte_text_path.filename().string() + ".figiss.construction"}};
    PrintConstructionTimeAndPeakMemory(construction_path, root.to_json().str());
  }
  {
    auto index_path {index_parent_path / std::filesystem::path{byte_text_path.filename().string() + ".figiss"}};
    index.Serialize(index_path);
    PrintSpace(index_path);
    PrintDetailedSpace(index_path, index);
  }
  return;
}

template <>
void ConstructAndSerialize
<sdsl::csa_wt<sdsl::wt_rlmn<>, 0xFFFF'FFFF, 0xFFFF'FFFF>>
(
  std::filesystem::path const& byte_text_path,
  std::filesystem::path const& construction_parent_path,
  std::filesystem::path const& index_parent_path
)
{
  sdsl::csa_wt<sdsl::wt_rlmn<>, 0xFFFF'FFFF, 0xFFFF'FFFF> index;
  std::cout << "--- construct & serialize rlfm ---\n";
  std::cout << "construct rlfm of " << std::filesystem::canonical(byte_text_path).string() << "\n";
  auto root {tdc::StatPhase("construction")};
  tdc::StatPhase::wrap
  (
    "rlfm",
    [&] ()
    {
      sdsl::int_vector<8> byte_text;
      sdsl::load_vector_from_file(byte_text, byte_text_path);
      sdsl::construct_im(index, byte_text);
    }
  );
  {
    auto construction_path {construction_parent_path / std::filesystem::path{byte_text_path.filename().string() + ".rlfm.construction"}};
    PrintConstructionTimeAndPeakMemory(construction_path, root.to_json().str());
  }
  {
    auto index_path {index_parent_path / std::filesystem::path{byte_text_path.filename().string() + ".rlfm"}};
    std::ofstream out {index_path};
    std::cout << "serialize rlfm to " << std::filesystem::canonical(index_path).string() << "\n";
    index.serialize(out);
    PrintSpace(index_path);
  }
  return;
}

template <>
void ConstructAndSerialize
<sdsl::csa_wt<wt_fbb<>, 0xFFFF'FFFF, 0xFFFF'FFFF>>
(
  std::filesystem::path const& byte_text_path,
  std::filesystem::path const& construction_parent_path,
  std::filesystem::path const& index_parent_path
)
{
  sdsl::csa_wt<wt_fbb<>, 0xFFFF'FFFF, 0xFFFF'FFFF> index;
  std::cout << "--- construct & serialize faster-minuter ---\n";
  std::cout << "construct faster-minuter of " << std::filesystem::canonical(byte_text_path).string() << "\n";
  auto root {tdc::StatPhase("construction")};
  tdc::StatPhase::wrap
  (
    "faster-minuter",
    [&] ()
    {
      sdsl::int_vector<8> byte_text;
      sdsl::load_vector_from_file(byte_text, byte_text_path);
      sdsl::construct_im(index, byte_text);
    }
  );
  {
    auto construction_path {construction_parent_path / std::filesystem::path{byte_text_path.filename().string() + ".faster-minuter.construction"}};
    PrintConstructionTimeAndPeakMemory(construction_path, root.to_json().str());
  }
  {
    auto index_path {index_parent_path / std::filesystem::path{byte_text_path.filename().string() + ".faster-minuter"}};
    std::ofstream out {index_path};
    std::cout << "serialize faster-minuter to " << std::filesystem::canonical(index_path).string() << "\n";
    index.serialize(out);
    PrintSpace(index_path);
  }
  return;
}

template <>
void ConstructAndSerialize
<CSA::RLCSA>
(
  std::filesystem::path const& byte_text_path,
  std::filesystem::path const& construction_parent_path,
  std::filesystem::path const& index_parent_path
)
{
  std::cout << "--- construct & serialize rlcsa ---\n";
  std::cout << "construct rlcsa of " << std::filesystem::canonical(byte_text_path).string() << "\n";
  auto root {tdc::StatPhase("construction")};
  tdc::StatPhase::wrap
  (
    "rlcsa",
    [&] ()
    {
      sdsl::int_vector<8> byte_text;
      sdsl::load_vector_from_file(byte_text, byte_text_path);
      sdsl::append_zero_symbol(byte_text);
      CSA::RLCSA index((CSA::uchar*)byte_text.data(), std::size(byte_text), 32, 0, 1, false);
    }
  );
  {
    auto construction_path {construction_parent_path / std::filesystem::path{byte_text_path.filename().string() + ".rlcsa.construction"}};
    PrintConstructionTimeAndPeakMemory(construction_path, root.to_json().str());
  }
  {
    auto index_path {index_parent_path / std::filesystem::path{byte_text_path.filename().string() + ".rlcsa.array"}};
    sdsl::int_vector<8> byte_text;
    sdsl::load_vector_from_file(byte_text, byte_text_path);
    sdsl::append_zero_symbol(byte_text);
    CSA::RLCSA index((CSA::uchar*)byte_text.data(), std::size(byte_text), 32, 0, 1, false);
    index.writeTo(byte_text_path.string());
    {
      auto index_from_path {std::filesystem::path{byte_text_path.string() + ".rlcsa.array"}};
      auto index_parameters_from_path {byte_text_path.parent_path() / std::filesystem::path{index_from_path.stem().string() + ".parameters"}};
      auto index_parameters_to_path {index_parent_path / std::filesystem::path{index_path.stem().string() + ".parameters"}};
      auto code {system(("mv " + index_from_path.string() + " " + index_path.string()).c_str())};
      if (WIFSIGNALED(code) && (WTERMSIG(code) == SIGINT || WTERMSIG(code) == SIGQUIT))
      {
        throw std::runtime_error("\033[31mfailed at moving rlcsa");
      }
      code = system(("mv " + index_parameters_from_path.string() + " " + index_parameters_to_path.string()).c_str());
      if (WIFSIGNALED(code) && (WTERMSIG(code) == SIGINT || WTERMSIG(code) == SIGQUIT))
      {
        throw std::runtime_error("\033[31mfailed at moving rlcsa parameters");
      }
      std::cout << "serialize rlcsa to " << std::filesystem::canonical(index_path).string() << "\n";
      std::cout << "write rlcsa parameters to " << std::filesystem::canonical(index_parameters_to_path).string() << "\n";
    }
    PrintSpace(index_path);
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
void Count (sdsl::csa_wt<wt_fbb<>, 0xFFFF'FFFF, 0xFFFF'FFFF>& index, Iterator begin, Iterator end)
{
  sdsl::count(index, begin, end);
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
      if (!time_path.parent_path().empty() && !std::filesystem::exists(time_path.parent_path()))
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
        out << "average counting time per character," << ProperTimeRepresentation(duration) << "\n";
      }
      else
      {
        out << "average counting time per character (ns)," << duration << "\n";
      }
    }
  }
  return;
}
}
