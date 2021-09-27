#pragma once

#include <sdsl/wavelet_trees.hpp>

#include "byte_alphabet.h"
#include "grammar_xbwt_trie.h"
#include "integer_alphabet.h"
#include "sparse_prefix_sum.h"

namespace figiss
{
class Index
{
public:

  static constexpr uint8_t kS {1};
  static constexpr uint8_t kL {0};

  Index () = default;
  Index (std::filesystem::path const& byte_text_path);
  Index (Index const&) = delete;
  Index (Index&&);
  Index& operator= (Index const&) = delete;
  Index& operator= (Index&&);

  void Swap (Index&);

  uint64_t Serialize
  (
    std::filesystem::path const& index_path,
    std::shared_ptr<SpaceNode> parent = nullptr
  );
  void Load (std::filesystem::path const& index_path);

  template <typename Iterator>
  uint64_t Count (Iterator begin, Iterator end);

  friend std::ostream& operator<< (std::ostream& out, Index const& index);

private:

  void CalculateRunLengthText
  (
    sdsl::int_vector<8> const& byte_text,
    sdsl::int_vector<>& run_length_text
  );
  void InsertGrammarRulesIntoGrammarTrie
  (
    sdsl::int_vector<> const& run_length_text,
    GrammarTrie& grammar_trie
  );
  void CalculateRunLengthAlphabetAndLexRankOnGrammarTrie
  (
    IntegerAlphabet& run_length_alphabet,
    GrammarTrie& grammar_trie
  );
  void CalculateColexText
  (
    sdsl::int_vector<> const& run_length_text,
    sdsl::int_vector<>& colex_text
  );
  void CalculateLexRankBucketOffsets (sdsl::int_vector<> const& colex_text);
  void CalculateColexBwt (sdsl::int_vector<> const& colex_text);
  template <typename Iterator>
  bool CalculateRunLengthPattern (Iterator begin, Iterator end, sdsl::int_vector<>& run_length_pattern);
  template <typename Iterator>
  void CalculateSlFactor (Iterator const rend, Iterator& rfirst, Iterator& rlast);
  template <typename Iterator, typename Range>
  Iterator ClassifyPattern (Iterator rbegin, Iterator rend, Range& pattern_range_s, Range& pattern_range_l);
  template <typename Iterator, typename Range>
  void BackwardSearchSuffix (Iterator rbegin, Iterator rend, Range& pattern_range);
  template <typename Iterator, typename Range>
  void BackwardSearchInfix (Iterator rbegin, Iterator rend, Range& pattern_range);
  template <typename Iterator, typename Range>
  void BackwardSearchPrefix (Iterator rbegin, Iterator rend, Range& pattern_range);
  template <typename Iterator, typename Range>
  void BackwardSearchInfix (Iterator rbegin, Iterator rend, Range& pattern_range_s, Range& pattern_range_l);
  template <typename Iterator, typename Range>
  void BackwardSearchPrefix (Iterator rbegin, Iterator rend, Range& pattern_range_s, Range& pattern_range_l);
  template <typename Range>
  uint64_t OrthogonalRangeCount (Range const pattern_range, Range const colex_rank_range);

  template <typename Range> // [,]
  inline bool IsNotEmptyRange (Range const& range) const
  {
    return (std::get<0>(range) <= std::get<1>(range));
  }

  template <typename Range> // [,]
  inline uint64_t RangeSize (Range const& range) const
  {
    if (std::get<0>(range) <= std::get<1>(range))
    {
      return (std::get<1>(range) - std::get<0>(range) + 1);
    }
    return 0;
  }

  ByteAlphabet byte_alphabet_;
  GrammarXbwtTrie grammar_xbwt_trie_;
  SparsePrefixSum lex_rank_bucket_offsets_;
  sdsl::wt_rlmn
  <
    sdsl::sd_vector<>,
    typename sdsl::sd_vector<>::rank_1_type,
    typename sdsl::sd_vector<>::select_1_type,
    sdsl::wt_ap<>
  >
  colex_bwt_;

};

Index::Index (std::filesystem::path const& byte_text_path)
{
  std::cout << "construct index of " << std::filesystem::canonical(byte_text_path) << "\n";
  sdsl::int_vector<> colex_text;
  {
    sdsl::int_vector<> run_length_text;
    {
      sdsl::int_vector<8> byte_text;
      sdsl::load_vector_from_file(byte_text, byte_text_path);
      {
        byte_alphabet_ = decltype(byte_alphabet_)(byte_text);
        // std::cout << byte_alphabet_;
      }
      {
        CalculateRunLengthText(byte_text, run_length_text);
        // {
        //   auto length_width {run_length_text.width() - byte_alphabet_.GetEffectiveAlphabetWidth()};
        //   uint64_t const divisor {1ULL << length_width};
        //   for (auto const& integer : run_length_text)
        //   {
        //     std::cout << "(" << (integer / divisor) << "," << (integer % divisor) << ")";
        //   }
        //   std::cout << "\n";
        // }
      }
    }
    {
      auto length_width {static_cast<uint8_t>(run_length_text.width() - byte_alphabet_.GetEffectiveAlphabetWidth())};
      GrammarTrie grammar_trie {length_width};
      InsertGrammarRulesIntoGrammarTrie(run_length_text, grammar_trie);
      // std::cout << grammar_trie;
      IntegerAlphabet run_length_alphabet;
      CalculateRunLengthAlphabetAndLexRankOnGrammarTrie(run_length_alphabet, grammar_trie);
      // std::cout << grammar_trie;
      // std::cout << run_length_alphabet;
      grammar_xbwt_trie_ = decltype(grammar_xbwt_trie_)(grammar_trie, std::move(run_length_alphabet));
      // std::cout << grammar_xbwt_trie_;
    }
    {
      CalculateColexText(run_length_text, colex_text);
      // Print(colex_text, std::cout);
    }
  }
  {
    CalculateLexRankBucketOffsets(colex_text);
    // std::cout << lex_rank_bucket_offsets_;
  }
  {
    CalculateColexBwt(colex_text);
    // std::cout << colex_bwt_ << "\n";
  }
}

Index::Index (Index&& index)
{
  if (this != &index)
  {
    this->Swap(index);
  }
}

Index& Index::operator= (Index&& index)
{
  if (this != &index)
  {
    Index temp {std::move(index)};
    this->Swap(temp);
  }
  return *this;
}

void Index::Swap (Index& index)
{
  if (this != &index)
  {
    byte_alphabet_.Swap(index.byte_alphabet_);
    grammar_xbwt_trie_.Swap(index.grammar_xbwt_trie_);
  }
  return;
}

uint64_t Index::Serialize
(
  std::filesystem::path const& index_path,
  std::shared_ptr<SpaceNode> root
)
{
  if
  (
    !index_path.parent_path().filename().string().empty()
    && !std::filesystem::exists(index_path.parent_path())
  )
  {
    std::filesystem::create_directories(index_path.parent_path());
  }
  std::fstream out {index_path, std::ios_base::out | std::ios_base::trunc};
  std::cout << "serialize index to" << std::filesystem::canonical(index_path) << "\n";
  uint64_t size_in_bytes {};
  if (!root)
  {
    byte_alphabet_.Serialize(out);
    grammar_xbwt_trie_.Serialize(out);
    lex_rank_bucket_offsets_.Serialize(out);
    colex_bwt_.serialize(out);
  }
  else
  {
    byte_alphabet_.Serialize(out, root, "byte_alphabet_");
    grammar_xbwt_trie_.Serialize(out, root, "grammar_xbwt_trie_");
    lex_rank_bucket_offsets_.Serialize(out, root, "lex_rank_bucket_offsets_");
    root->AddLeaf("colex_bwt_", colex_bwt_.serialize(out));
    size_in_bytes = root->GetSizeInBytes();
  }
  return size_in_bytes;
}

void Index::Load (std::filesystem::path const& index_path)
{
  std::ifstream in {index_path};
  std::cout << "load index from " << std::filesystem::canonical(index_path) << "\n";
  byte_alphabet_.Load(in);
  grammar_xbwt_trie_.Load(in);
  lex_rank_bucket_offsets_.Load(in);
  colex_bwt_.load(in);
  return;
}

template <typename Iterator>
uint64_t Index::Count (Iterator begin, Iterator end)
{
  auto pattern_range_s {std::make_pair<uint64_t, uint64_t>(1, 0)};
  auto pattern_range_l {std::make_pair<uint64_t, uint64_t>(1, 0)};
  sdsl::int_vector<> run_length_pattern;
  if (CalculateRunLengthPattern(begin, end, run_length_pattern))
  {
    // {
    //   uint64_t const divisor {1ULL << grammar_xbwt_trie_.GetLengthWidth()};
    //   for (auto const& integer : run_length_pattern)
    //   {
    //     std::cout << "(" << (integer / divisor) << "," << (integer % divisor) << ")";
    //   }
    //   std::cout << "\n";
    // }
    auto rbegin {std::prev(std::end(run_length_pattern))};
    auto rend {std::prev(std::begin(run_length_pattern))};
    auto rfirst {rbegin};
    auto rlast {rbegin};
    rlast = ClassifyPattern(rbegin, rend, pattern_range_s, pattern_range_l);
    if (rlast != rend)
    {
      while (IsNotEmptyRange(pattern_range_s) || IsNotEmptyRange(pattern_range_l))
      {
        CalculateSlFactor(rend, rfirst, rlast);
        if (rlast != rend)
        {
          BackwardSearchInfix(rfirst, rlast, pattern_range_s, pattern_range_l);
        }
        else
        {
          BackwardSearchPrefix(rfirst, rlast, pattern_range_s, pattern_range_l);
          break;
        }
      }
    }
  }
  return (RangeSize(pattern_range_s) + RangeSize(pattern_range_l));
}

std::ostream& operator<< (std::ostream& out, Index const& index)
{
  out << "byte_alphabet_:\n";
  out << index.byte_alphabet_ << "\n";
  out << "grammar_xbwt_trie_:\n";
  out << index.grammar_xbwt_trie_ << "\n";
  out << "lex_rank_bucket_offsets_:";
  out << index.lex_rank_bucket_offsets_ << "\n";
  out << "colex_bwt_:";
  Print(index.colex_bwt_, out);
  return out;
}

void Index::CalculateRunLengthText
(
  sdsl::int_vector<8> const& byte_text,
  sdsl::int_vector<>& run_length_text
)
{
  uint8_t prev_byte {byte_text[0]};
  uint64_t run_length_text_size {1};
  uint64_t length {1};
  uint64_t max_length {1};
  for (auto const byte : byte_text)
  {
    if (prev_byte != byte)
    {
      if (max_length < length)
      {
        max_length = length;
      }
      prev_byte = byte;
      length = 1;
      ++run_length_text_size;
    }
    else
    {
      ++length;
    }
  }
  if (max_length < length)
  {
    max_length = length;
  }
  {
    auto length_width {sdsl::bits::hi(max_length) + 1};
    run_length_text.width(byte_alphabet_.GetEffectiveAlphabetWidth() + length_width);
    run_length_text.resize(run_length_text_size);
    auto it {std::begin(run_length_text)};
    uint64_t const divisor {1ULL << length_width};
    prev_byte = byte_text[0];
    length = 0;
    for (auto const byte : byte_text)
    {
      if (prev_byte != byte)
      {
        *it++ = (byte_alphabet_[prev_byte] * divisor) + length;
        prev_byte = byte;
        length = 1;
      }
      else
      {
        ++length;
      }
    }
    *it = (byte_alphabet_[prev_byte] * divisor) + length;
  }
}

void Index::InsertGrammarRulesIntoGrammarTrie
(
  sdsl::int_vector<> const& run_length_text,
  GrammarTrie& grammar_trie
)
{
  auto text_rend {std::prev(std::begin(run_length_text))};
  auto text_it {std::prev(std::end(run_length_text), 2)};
  auto next_symbol {*std::prev(std::end(run_length_text))};
  uint8_t sl_type {};
  uint8_t next_sl_type {Index::kL};
  auto sl_factor_end {std::end(run_length_text)};
  while (text_it != text_rend)
  {
    if (*text_it == next_symbol)
    {
      sl_type = next_sl_type;
    }
    else if (*text_it < next_symbol)
    {
      sl_type = Index::kS;
    }
    else
    {
      sl_type = Index::kL;
    }
    if ((sl_type == Index::kL) && (next_sl_type == Index::kS))
    {
      auto sl_factor_begin {std::next(text_it)};
      grammar_trie.Insert(sl_factor_begin, sl_factor_end);
      sl_factor_end = sl_factor_begin;
    }
    next_symbol = *text_it--;
    next_sl_type = sl_type;
  }
  grammar_trie.Insert(std::begin(run_length_text), sl_factor_end);
  return;
}

void Index::CalculateRunLengthAlphabetAndLexRankOnGrammarTrie
(
  IntegerAlphabet& run_length_alphabet,
  GrammarTrie& grammar_trie
)
{
  std::set<uint64_t> alphabet;
  uint64_t lex_rank {};
  std::deque<std::shared_ptr<GrammarTrie::Node>> nodes;
  nodes.emplace_back(grammar_trie.GetRoot());
  while (!nodes.empty())
  {
    auto node {nodes.back()};
    if (node->label)
    {
      alphabet.insert(node->label);
    }
    if (node->lex_rank)
    {
      node->lex_rank = ++lex_rank;
    }
    nodes.pop_back();
    for (auto it {std::rbegin(node->children)}; it != std::rend(node->children); ++it)
    {
      nodes.emplace_back(std::get<1>(*it));
    }
  }
  run_length_alphabet = IntegerAlphabet(alphabet);
  return;
}

void Index::CalculateColexText
(
  sdsl::int_vector<> const& run_length_text,
  sdsl::int_vector<>& colex_text
)
{
  std::deque<uint64_t> temp_colex_text(1, 0);
  auto rend {std::prev(std::begin(run_length_text))};
  auto it {std::prev(std::end(run_length_text), 2)};
  auto next_symbol {*std::prev(std::end(run_length_text))};
  uint8_t sl_type {};
  uint8_t next_sl_type {Index::kL};
  auto sl_factor_end {std::end(run_length_text)};
  while (it != rend)
  {
    if (*it == next_symbol)
    {
      sl_type = next_sl_type;
    }
    else if (*it < next_symbol)
    {
      sl_type = Index::kS;
    }
    else
    {
      sl_type = Index::kL;
    }
    if ((sl_type == Index::kL) && (next_sl_type == Index::kS))
    {
      auto sl_factor_begin {std::next(it)};
      temp_colex_text.push_front(grammar_xbwt_trie_.GetColexRank(sl_factor_begin, sl_factor_end));
      sl_factor_end = sl_factor_begin;
    }
    next_symbol = *it--;
    next_sl_type = sl_type;
  }
  temp_colex_text.push_front(grammar_xbwt_trie_.GetColexRank(std::next(rend), sl_factor_end));
  {
    colex_text.width(grammar_xbwt_trie_.GetColexAlphabetWidth());
    colex_text.resize(std::size(temp_colex_text));
    std::copy(std::begin(temp_colex_text), std::end(temp_colex_text), std::begin(colex_text));
  }
  return;
}

void Index::CalculateLexRankBucketOffsets (sdsl::int_vector<> const& colex_text)
{
  std::map<uint64_t, uint64_t> lex_rank_counts;
  for (auto const& colex_rank : colex_text)
  {
    auto lex_rank {grammar_xbwt_trie_.ColexToLexRank(colex_rank)};
    auto it {lex_rank_counts.find(lex_rank)};
    if (it != lex_rank_counts.end())
    {
      ++std::get<1>(*it);
    }
    else
    {
      lex_rank_counts[lex_rank] = 1;
    }
  }
  std::deque<uint64_t> cumulative_counts;
  for (auto const& pair : lex_rank_counts)
  {
    cumulative_counts.push_back(std::get<1>(pair));
  }
  std::partial_sum(std::begin(cumulative_counts), std::end(cumulative_counts), std::begin(cumulative_counts));
  lex_rank_bucket_offsets_ = decltype(lex_rank_bucket_offsets_)(cumulative_counts);
  return;
}

void Index::CalculateColexBwt (sdsl::int_vector<> const& colex_text)
{
  sdsl::int_vector<> buffer;
  sdsl::qsufsort::construct_sa(buffer, colex_text);
  for (auto it {std::begin(buffer)}; it != std::end(buffer); ++it)
  {
    if (*it != 0)
    {
      *it = colex_text[*it - 1];
    }
  }
  sdsl::construct_im(colex_bwt_, buffer);
  return;
}

template <typename Iterator>
bool Index::CalculateRunLengthPattern (Iterator begin, Iterator end, sdsl::int_vector<>& run_length_pattern)
{
  uint8_t prev_byte {};
  uint64_t run_length_pattern_size {};
  for (auto it {begin}; it != end; ++it)
  {
    if (*it != prev_byte)
    {
      if (byte_alphabet_[*it] == 0)
      {
        return false;
      }
      prev_byte = *it;
      ++run_length_pattern_size;
    }
  }
  run_length_pattern.width(grammar_xbwt_trie_.GetAlphabetWidth());
  run_length_pattern.resize(run_length_pattern_size);
  auto run_length_pattern_it {std::begin(run_length_pattern)};
  {
    uint64_t const divisor {1ULL << grammar_xbwt_trie_.GetLengthWidth()};
    prev_byte = 0;
    uint64_t length {};
    for (auto it {begin}; it != end; ++it)
    {
      if (*it != prev_byte)
      {
        if (prev_byte)
        {
          *run_length_pattern_it++ = byte_alphabet_[prev_byte] * divisor + length;
        }
        prev_byte = *it;
        length = 1;
      }
      else
      {
        ++length;
      }
    }
    *run_length_pattern_it = byte_alphabet_[prev_byte] * divisor + length;
  }
  return true;
}

template <typename Iterator>
void Index::CalculateSlFactor (Iterator const rend, Iterator& rfirst, Iterator& rlast)
{
  auto next_sl_type {Index::kL};
  rfirst = rlast--;
  while ((rlast != rend) && !((next_sl_type == Index::kS) && (*rlast > *std::next(rlast))))
  {
    if((next_sl_type == Index::kL) && (*rlast < *std::next(rlast)))
    {
      next_sl_type = Index::kS;
    }
    --rlast;
  }
  return;
}

template <typename Iterator, typename Range>
Iterator Index::ClassifyPattern (Iterator rbegin, Iterator rend, Range& pattern_range_s, Range& pattern_range_l)
{
  auto rfirst {rbegin};
  auto rlast {rbegin};
  if (rbegin != rend)
  {
    if (std::prev(rbegin) != rend)
    {
      if (*std::prev(rbegin) > *rbegin)
      {
        rlast = std::prev(rbegin);
        BackwardSearchSuffix(rfirst, rlast, pattern_range_s);
        if (IsNotEmptyRange(pattern_range_s))
        {
          CalculateSlFactor(rend, rfirst, rlast);
          if (rlast != rend)
          {
            BackwardSearchInfix(rfirst, rlast, pattern_range_s);
          }
          else
          {
            BackwardSearchPrefix(rfirst, rlast, pattern_range_s);
          }
        }
      }
      else
      {
        CalculateSlFactor(rend, rfirst, rlast);
      }
      if (rlast != rend)
      {
        BackwardSearchSuffix(rfirst, rend, pattern_range_l);
      }
      else
      {
        pattern_range_l = {1, grammar_xbwt_trie_.Count(std::next(rend), std::next(rbegin))};
      }
    }
    else
    {
      pattern_range_l = {1, grammar_xbwt_trie_.Count(std::next(rend), std::next(rbegin))};
    }
  }
  return rlast;
}

template <typename Iterator, typename Range>
void Index::BackwardSearchSuffix (Iterator rbegin, Iterator rend, Range& pattern_range)
{
  auto lex_rank_range {grammar_xbwt_trie_.GetLexRankRange(std::next(rend), std::next(rbegin))};
  if (IsNotEmptyRange(lex_rank_range))
  {
    pattern_range =
    {
      lex_rank_bucket_offsets_[std::get<0>(lex_rank_range)],
      lex_rank_bucket_offsets_[std::get<1>(lex_rank_range) + 1] - 1
    };
  }
  else
  {
    pattern_range = {1, 0};
  }
  return;
}

template <typename Iterator, typename Range>
void Index::BackwardSearchInfix (Iterator rbegin, Iterator rend, Range& pattern_range)
{
  auto colex_rank {grammar_xbwt_trie_.GetColexRank(std::next(rend), std::next(rbegin))};
  if (colex_rank)
  {
    auto begin_offset {lex_rank_bucket_offsets_[grammar_xbwt_trie_.ColexToLexRank(colex_rank)]};
    pattern_range =
    {
      begin_offset + colex_bwt_.rank(std::get<0>(pattern_range), colex_rank),
      begin_offset + colex_bwt_.rank(std::get<1>(pattern_range) + 1, colex_rank) - 1
    };
  }
  else
  {
    pattern_range = {1, 0};
  }
  return;
}

template <typename Iterator, typename Range>
void Index::BackwardSearchPrefix (Iterator rbegin, Iterator rend, Range& pattern_range)
{
  auto colex_rank_range {grammar_xbwt_trie_.GetColexRankRange(std::next(rend), std::next(rbegin))};
  if (IsNotEmptyRange(colex_rank_range))
  {
    pattern_range = {1, OrthogonalRangeCount(pattern_range, colex_rank_range)};
  }
  else
  {
    pattern_range = {1, 0};
  }
  return;
}

template <typename Iterator, typename Range>
void Index::BackwardSearchInfix (Iterator rbegin, Iterator rend, Range& pattern_range_s, Range& pattern_range_l)
{
  auto colex_rank {grammar_xbwt_trie_.GetColexRank(std::next(rend), std::next(rbegin))};
  if (colex_rank)
  {
    auto begin_offset {lex_rank_bucket_offsets_[grammar_xbwt_trie_.ColexToLexRank(colex_rank)]};
    if (IsNotEmptyRange(pattern_range_s))
    {
      pattern_range_s =
      {
        begin_offset + colex_bwt_.rank(std::get<0>(pattern_range_s), colex_rank),
        begin_offset + colex_bwt_.rank(std::get<1>(pattern_range_s) + 1, colex_rank) - 1
      };
    }
    if (IsNotEmptyRange(pattern_range_l))
    {
      pattern_range_l =
      {
        begin_offset + colex_bwt_.rank(std::get<0>(pattern_range_l), colex_rank),
        begin_offset + colex_bwt_.rank(std::get<1>(pattern_range_l) + 1, colex_rank) - 1
      };
    }
  }
  else
  {
    pattern_range_s = pattern_range_l = {1, 0};
  }
  return;
}

template <typename Iterator, typename Range>
void Index::BackwardSearchPrefix (Iterator rbegin, Iterator rend, Range& pattern_range_s, Range& pattern_range_l)
{
  auto colex_rank_range {grammar_xbwt_trie_.GetColexRankRange(std::next(rend), std::next(rbegin))};
  if (IsNotEmptyRange(colex_rank_range))
  {
    if (IsNotEmptyRange(pattern_range_s))
    {
      pattern_range_s = {1, OrthogonalRangeCount(pattern_range_s, colex_rank_range)};
    }
    if (IsNotEmptyRange(pattern_range_l))
    {
      pattern_range_l = {1, OrthogonalRangeCount(pattern_range_l, colex_rank_range)};
    }
  }
  else
  {
    pattern_range_s = pattern_range_l = {1, 0};
  }
  return;
}

template <typename Range>
uint64_t Index::OrthogonalRangeCount (Range const pattern_range, Range const colex_rank_range)
{
  uint64_t count {};
  for (auto colex_rank {std::get<0>(colex_rank_range)}; colex_rank <= std::get<1>(colex_rank_range); ++colex_rank)
  {
    count +=
    (
      colex_bwt_.rank(std::get<1>(pattern_range) + 1, colex_rank)
      - colex_bwt_.rank(std::get<0>(pattern_range), colex_rank)
    );
  }
  return count;
}
}
