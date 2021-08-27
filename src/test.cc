#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>

#include "benchmark.h"
#include "index.h"

void Usage ()
{
  std::cout << "Usage: ./test [k] [text path]\n"
  << "\tconstruct index of text of parameter [k], T, at [text path]\n"
  << "\t[k] can only be integer in [1,8]\n"
  << "\ttest correctness of conunting query on pattern of size being power of 2 <= |T|\n";
  return;
}

int main (int argc, char **argv)
{
  if (argc != 3)
  {
    Usage();
    return -1;
  }
  else
  {
    auto k {std::stoll(argv[1])};
    auto text_path {std::filesystem::path(argv[2])};
    if (k < 1 || k > 8)
    {
      Usage();
      return -1;
    }
    else if (k == 1)
    {
      figiss::Index<1> index;
      figiss::TestCount(text_path, index);
    }
    else if (k == 2)
    {
      figiss::Index<2> index;
      figiss::TestCount(text_path, index);
    }
    else if (k == 3)
    {
      figiss::Index<3> index;
      figiss::TestCount(text_path, index);
    }
    else if (k == 4)
    {
      figiss::Index<> index;
      figiss::TestCount(text_path, index);
    }
    else if (k == 5)
    {
      figiss::Index<5> index;
      figiss::TestCount(text_path, index);
    }
    else if (k == 6)
    {
      figiss::Index<6> index;
      figiss::TestCount(text_path, index);
    }
    else if (k == 7)
    {
      figiss::Index<7> index;
      figiss::TestCount(text_path, index);
    }
    else if (k == 8)
    {
      figiss::Index<8> index;
      figiss::TestCount(text_path, index);
    }
  }
  return 0;
}
