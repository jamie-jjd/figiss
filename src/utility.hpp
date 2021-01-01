#ifndef UTILITY_HPP_
#define UTILITY_HPP_

#include <string>

namespace gci
{
namespace util
{
std::string basename (std::string const path)
{
  auto basename_begin_dist {path.find_last_of("/")};
  if (basename_begin_dist != std::string::npos)
  {
    ++basename_begin_dist;
  }
  else
  {
    basename_begin_dist = 0;
  }
  return path.substr(basename_begin_dist);
}
}
}

#endif
