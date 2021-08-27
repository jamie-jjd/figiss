# 🗂️  FIGISS: FM-Indexing Grammars Induced by Suffix Sorting for Long Patterns

The FIGISS index is a full-text self-index capable of counting the occurrences of patterns.
The index is the FM index built upon the Burrows-Wheeler transform on the right hand side of the start symbol of a grammar.
In our case, this grammar (see the references) is based on the SA-IS algorithm, but it can be exchanged with other grammars that have locally consistent properties.

## Usage

### 🚀 Complete Test Run

A complete test-run of our index is done with the following lines:

```bash:
cd figiss/src
make
./figiss cs 4 ../data/text.txt ../data/index
./figiss lc 4 ../data/index ../data/pattern.txt
```

The output should be `21`, which is the number of occurrences of `../data/pattern.txt` in `../data/text.txt`.


### ✔️ Prerequisites

To compile this project, you need the following tools:

 - `make`
 - a recent `g++` with C++17 support
 - [sdsl-lite](https://github.com/simongog/sdsl-lite)


### ⚙️ Compilation

The FIGISS index consists in a single executable called `figiss`.
It can be compiled with:

```bash:
cd figiss/src
make
```

Counting query testing can be done by an executable called `test`.
It can be compiled with:
```bash:
cd src
make test
```

### 🏗️ Index Construction

Syntax:

```bash:
./figiss cs [k] [text path] [index path]
```

 - constructs the index of the text at `[text path]` and serializes it to `[index path]`
 - `k` must be an integer in [1..8] as described in the paper. Good choices are `4` and `7`.

Example:

```bash:
./figiss cs 4 ../data/text.txt ../data/index
```

Output:

```bash:
construct index of "/figiss/data/text.txt"
serialize index to "/figiss/data/index
```

### 🔎 Counting Query

```bash:
./figiss lc [k] [index path] [pattern path]
```

 - loads the index from `[index path]` and report number of occurrences of pattern at `[pattern path]`
 - ⚠️ the same `k` as during the construction must be used
 - ⚠️ all characters of the file `[pattern path]` are considered, even if the file ends with a new line character (which is commonly appended by text editors!)

Example:

```bash:
./figiss lc 4 ../data/index ../data/pattern.txt
```

Output:

```bash:
load index from ".../figiss/data/index"
21
```

### Testing Counting Query

```bash:
./test [k] [text path]
```

- constructs index of text at `[text path]` with parameter `k` and tests correctness of conunting query on pattern of size being power of 2 <= text size.
- `k` must be an integer in [1..8] as described in the paper.

## 📚 References

- Daniel Saad Nogueira Nunes, Felipe A. Louza, Simon Gog, Mauricio Ayala-Rincón, Gonzalo Navarro: [A Grammar Compression Algorithm Based on Induced Suffix Sorting. DCC 2018: 42-51](https://doi.org/10.1109/DCC.2018.00012)


## Using as a Library

You can also use FIGISS as a library in the following way:

### Sample Code for Constructing Index

```c++:
#include "index.h"
int main (int argc, char **argv)
{  
  if (argc == 2)
  {
    auto text_path {std::filesystem::path{argv[1]}};
    figiss::Index<> index {text_path}; // figiss::Index<?>, ? can be replaced by [1..8] (by default 4)  
  }
  return 0;
}
```

### Sample Code for Serializing Index

```c++:
#include "index.h"
int main (int argc, char **argv)
{
  if (argc == 3)
  {
    auto text_path {std::filesystem::path{argv[1]}};
    auto index_path {std::filesystem::path{argv[2]}};
    figiss::Index<> index {text_path};
    index.Serialize(index_path);
  }
  return 0;
}
```

### Sample Code for Loading Index

```c++:
#include "index.h"
int main (int argc, char **argv)
{
  if (argc == 2)
  {
    auto index_path {std::filesystem::path{argv[1]}};
    figiss::Index<> index; // parameter should be matched
    index.Load(index_path);
  }
  return 0;
}
```

### Sample Code for Counting

```c++:
#include "index.h"
int main (int argc, char **argv)
{
  if (argc == 2)
  {
    auto text_path {std::filesystem::path{argv[1]}};
    figiss::Index<> index {text_path};
    std::string pattern {"your favorite pattern"};
    std::cout << index.Count(std::begin(pattern), std::end(pattern)) << "\n";
  }
  return 0;
}
```

### Sample Code for Testing Counting

```c++:
#include "benchmark.h"
#include "index.h"
int main (int argc, char **argv)
{
  if (argc == 2)
  {
    auto text_path {std::filesystem::path{argv[1]}};
    figiss::Index<> index; // figiss::Index<?>, ? can be replaced by [1..8] (by default 4)
    figiss::TestCount(text_path, index);
  }
  return 0;
}
```
