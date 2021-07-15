#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>

#include "benchmark.h"
#include "grammar_compressed_index.h"

int main (int argc, char **argv)
{
  if (argc != 2)
  {
    throw std::runtime_error(std::string("usage: ") + argv[0] + " [text path]");
  }
  std::filesystem::path text_path {argv[1]};
  {
    gciis::Index index;
    sdsl::csa_wt<sdsl::wt_rlmn<>, 0xFFFF'FFFF, 0xFFFF'FFFF> rlfm;
    gciis::TestCount(text_path, index);
    gciis::PrintIndexSpace(text_path, index);
    gciis::MeasureCountingTime(text_path, index, rlfm);
  }
  {
    gciis::Index<5> index;
    sdsl::csa_wt<sdsl::wt_rlmn<>, 0xFFFF'FFFF, 0xFFFF'FFFF> rlfm;
    gciis::TestCount(text_path, index);
    gciis::PrintIndexSpace(text_path, index);
    gciis::MeasureCountingTime(text_path, index, rlfm);
  }
  {
    gciis::Index<6> index;
    sdsl::csa_wt<sdsl::wt_rlmn<>, 0xFFFF'FFFF, 0xFFFF'FFFF> rlfm;
    gciis::TestCount(text_path, index);
    gciis::PrintIndexSpace(text_path, index);
    gciis::MeasureCountingTime(text_path, index, rlfm);
  }
  {
    gciis::Index<7> index;
    sdsl::csa_wt<sdsl::wt_rlmn<>, 0xFFFF'FFFF, 0xFFFF'FFFF> rlfm;
    gciis::TestCount(text_path, index);
    gciis::PrintIndexSpace(text_path, index);
    gciis::MeasureCountingTime(text_path, index, rlfm);
  }
  {
    gciis::Index<8> index;
    sdsl::csa_wt<sdsl::wt_rlmn<>, 0xFFFF'FFFF, 0xFFFF'FFFF> rlfm;
    gciis::TestCount(text_path, index);
    gciis::PrintIndexSpace(text_path, index);
    gciis::MeasureCountingTime(text_path, index, rlfm);
  }
  return 0;
}
