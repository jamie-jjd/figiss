# 🗂️  GCIIS: FM-Indexing Grammars Induced by Suffix Sorting for Long Patterns

The GCIIS index is a full-text self-index capable of counting the occurrences of patterns. 
The index is the FM index built upon the Burrows-Wheeler transform on the right hand side of the start symbol of a grammar.
In our case, this grammar (see the refernces) is based on the SA-IS algorithm, but it can be exchanged with other grammars that have locally consistent properties.

## Usage

### 🚀 Prequisites

To compile this project, you need the following tools:

 - `make`
 - a recent `g++` with C++17 support
 - [sdsl-lite](https://github.com/simongog/sdsl-lite)


###  ⚙️ Compilation

The GCIIS index consists in a single executable called `gciis`.
It can be compiled with:

```bash:
git clone https://github.com/jamie-jjd/gciis
cd src
make
```

### 🏗️ Index Construction


Syntax:

```bash:
./gciis cs [k] [text path] [index path]
```

 - constructs the index of the text at `[text path]` and serializes it to `[index path]`
 - `k` must be an integer in [1..8] as described in the paper. Good choices are `4` and `7`.

### 🔎 Counting Query

```bash:
./gciis lc [k] [index path] [pattern path]
```

 - loads the index from `[index path]` and report number of occurences of pattern at `[pattern path]`
 - ⚠️ the same `k` as during the construction must be used
 - ⚠️ all characters of the file `[pattern path]` are considered, even if the file ends with a new line character (which is commonly appended by text editors!)

## 📚 References

- Daniel Saad Nogueira Nunes, Felipe A. Louza, Simon Gog, Mauricio Ayala-Rincón, Gonzalo Navarro: [A Grammar Compression Algorithm Based on Induced Suffix Sorting. DCC 2018: 42-51](https://doi.org/10.1109/DCC.2018.00012)


## Using as a Library

You can also use GCIIS as a library in the following way:

### Sample Code for Constructing Index

```c++:
#include "grammar_compressed_index.h"
int main (int argc, char **argv)
{  
  if (argc == 2) 
  {
    auto text_path {std::filesystem::path{argv[1]}};
    gciis::Index<> index {text_path}; // gciis::Index<?>, ? can be replaced by 1 ~ 8 (by default 4)  
  }
  return 0;
}
```

### Sample Code for Serializing Index

```c++:
#include "grammar_compressed_index.h"
int main (int argc, char **argv)
{
  if (argc == 3)
  {
    auto text_path {std::filesystem::path{argv[1]}};
    auto index_path {std::filesystem::path{argv[2]}};
    gciis::Index<> index {text_path};
    index.Serialize(index_path);
  }
  return 0;
}
```

### Sample Code for Loading Index

```c++:
#include "grammar_compressed_index.h"
int main (int argc, char **argv)
{
  if (argc == 2)
  {
    auto index_path {std::filesystem::path{argv[1]}};
    gciis::Index<> index; // parameter should be matched
    index.Load(index_path);
  }
  return 0;
}
```

### Sample Code for Counting

```c++:
#include "grammar_compressed_index.h"
int main (int argc, char **argv)
{
  if (argc == 2)
  {
    auto text_path {std::filesystem::path{argv[1]}};
    gciis::Index<> index {text_path};
    std::string pattern {"your favorite pattern"};
    std::cout << index.Count(std::begin(pattern), std::end(pattern)) << "\n";
  }
  return 0;
}
```
