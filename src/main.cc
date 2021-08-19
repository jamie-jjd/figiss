#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>

#include "grammar_compressed_index.h"

void Usage ()
{
  std::cout << "Usage: ./gciis\n"
  << "cs [k] [text path] [index path]\n"
  << "\tconstruct index of text at [text path] and serialize it to [index path]\n"
  << "\tk can only be integer in [1, 8]\n"
  << "lc [k] [index path] [pattern path]\n"
  << "\tload index from [index path] and report number of occurences of pattern at [pattern path]\n"
  << "\t k must be matched with the k of index\n";
  return;
}

int main (int argc, char **argv)
{
  if (argc < 2)
  {
    Usage();
    return -1;
  }
  else
  {
    auto option {std::string(argv[1])};
    if (option == "cs")
    {
      if (argc != 5)
      {
        Usage();
        return -1;
      }
      else
      {
        auto k {std::stoll(argv[2])};
        auto text_path {std::filesystem::path(argv[3])};
        auto index_path {std::filesystem::path(argv[4])};
        if (k < 1 || k > 8)
        {
          Usage();
          return -1;
        }
        else if (k == 1)
        {
          gciis::Index<1> index {text_path};
          index.Serialize(index_path);
        }
        else if (k == 2)
        {
          gciis::Index<2> index {text_path};
          index.Serialize(index_path);
        }
        else if (k == 3)
        {
          gciis::Index<3> index {text_path};
          index.Serialize(index_path);
        }
        else if (k == 4)
        {
          gciis::Index<> index {text_path};
          index.Serialize(index_path);
        }
        else if (k == 5)
        {
          gciis::Index<5> index {text_path};
          index.Serialize(index_path);
        }
        else if (k == 6)
        {
          gciis::Index<6> index {text_path};
          index.Serialize(index_path);
        }
        else if (k == 7)
        {
          gciis::Index<7> index {text_path};
          index.Serialize(index_path);
        }
        else if (k == 8)
        {
          gciis::Index<8> index {text_path};
          index.Serialize(index_path);
        }
      }
    }
    else if (option == "lc")
    {
      if (argc != 5)
      {
        Usage();
        return -1;
      }
      else
      {
        auto k {std::stoll(argv[2])};
        auto index_path {std::filesystem::path(argv[3])};
        auto pattern_path {std::filesystem::path(argv[4])};
        sdsl::int_vector<8> pattern;
        sdsl::load_vector_from_file(pattern, pattern_path);
        if (k < 1 || k > 8)
        {
          Usage();
          return -1;
        }
        if (k == 1)
        {
          gciis::Index<1> index;
          index.Load(index_path);
          std::cout << index.Count(std::begin(pattern), std::end(pattern)) << "\n";
        }
        else if (k == 2)
        {
          gciis::Index<2> index;
          index.Load(index_path);
          std::cout << index.Count(std::begin(pattern), std::end(pattern)) << "\n";
        }
        else if (k == 3)
        {
          gciis::Index<3> index;
          index.Load(index_path);
          std::cout << index.Count(std::begin(pattern), std::end(pattern)) << "\n";
        }
        else if (k == 4)
        {
          gciis::Index<> index;
          index.Load(index_path);
          std::cout << index.Count(std::begin(pattern), std::end(pattern)) << "\n";
        }
        else if (k == 5)
        {
          gciis::Index<5> index;
          index.Load(index_path);
          std::cout << index.Count(std::begin(pattern), std::end(pattern)) << "\n";
        }
        else if (k == 6)
        {
          gciis::Index<6> index;
          index.Load(index_path);
          std::cout << index.Count(std::begin(pattern), std::end(pattern)) << "\n";
        }
        else if (k == 7)
        {
          gciis::Index<7> index;
          index.Load(index_path);
          std::cout << index.Count(std::begin(pattern), std::end(pattern)) << "\n";
        }
        else if (k == 8)
        {
          gciis::Index<8> index;
          index.Load(index_path);
          std::cout << index.Count(std::begin(pattern), std::end(pattern)) << "\n";
        }
      }
    }
    else
    {
      Usage();
      return -1;
    }
  }
  return 0;
}
