#include <sys/wait.h>

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>

#include "index.h"
#include "test.h"

void DecompressCompressedCorpus
(
  std::filesystem::path const& compressed_corpus_path,
  std::filesystem::path const& corpus_path
)
{
  if (!std::filesystem::exists(corpus_path))
  {
    std::filesystem::create_directories(corpus_path);
  }
  for (auto const& entry : std::filesystem::directory_iterator(compressed_corpus_path))
  {
    if (entry.is_regular_file())
    {
      auto compressed_text_path {entry.path()};
      auto decompressed_text_to_path {corpus_path / entry.path().stem()};
      if (!std::filesystem::exists(decompressed_text_to_path))
      {
        if (compressed_text_path.extension() == ".7z")
        {
          {
            auto code {system(("p7zip -k -d " + compressed_text_path.string()).c_str())};
            if (WIFSIGNALED(code) && (WTERMSIG(code) == SIGINT || WTERMSIG(code) == SIGQUIT))
            {
              throw std::runtime_error
              (
                "\033[31mfailed at decompressing " +
                std::filesystem::canonical(compressed_text_path).string() +
                "\033[0m\n"
              );
            }
          }
          {
            auto decompressed_text_from_path {decompressed_text_to_path.filename()};
            auto code {system(("mv " + decompressed_text_from_path.string() + " " + decompressed_text_to_path.string()).c_str())};
            if (WIFSIGNALED(code) && (WTERMSIG(code) == SIGINT || WTERMSIG(code) == SIGQUIT))
            {
              throw std::runtime_error("\033[31mfailed at moving decompressed text");
            }
          }
        }
        else if (compressed_text_path.extension() == ".xz")
        {
          {
            auto code {system(("xz -kdv -T0 " + compressed_text_path.string()).c_str())};
            if (WIFSIGNALED(code) && (WTERMSIG(code) == SIGINT || WTERMSIG(code) == SIGQUIT))
            {
              throw std::runtime_error
              (
                "\033[31mfailed at decompressing " +
                std::filesystem::canonical(compressed_text_path).string() +
                "\033[0m\n"
              );
            }
            {
              auto decompressed_text_from_path {compressed_corpus_path / compressed_text_path.stem()};
              auto code {system(("mv " + decompressed_text_from_path.string() + " " + decompressed_text_to_path.string()).c_str())};
              if (WIFSIGNALED(code) && (WTERMSIG(code) == SIGINT || WTERMSIG(code) == SIGQUIT))
              {
                throw std::runtime_error("\033[31mfailed at moving decompressed text");
              }
            }
          }
        }
      }
    }
  }
  return;
}

int main ()
{
  auto compressed_corpus_path {std::filesystem::path{"../data/compressed_corpus"}};
  auto corpus_path {std::filesystem::path{"../data/corpus"}};
  DecompressCompressedCorpus(compressed_corpus_path, corpus_path);
  std::cout << "test counting of index built from text in corpus (" << std::filesystem::canonical(corpus_path).string() << ")\n";
  for (auto const& entry : std::filesystem::directory_iterator(corpus_path))
  {
    if (entry.is_regular_file())
    {
      auto byte_text_path {entry.path()};
      auto index_path {std::filesystem::path{"../data/index/ours/" + entry.path().filename().string() + ".i"}};
      auto pattern_parent_path {std::filesystem::path{"../data/pattern"}};
      figiss::Index<> index {byte_text_path};
      index.Serialize(index_path);
      figiss::TestCounting(byte_text_path, index_path, pattern_parent_path);
    }
  }
  return 0;
}
