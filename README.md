# 🗂️  FIGISS: FM-Indexing Grammars Induced by Suffix Sorting for Long Patterns

The FIGISS index is a full-text self-index capable of counting the occurrences of patterns.
The index is the FM-index built upon the Burrows-Wheeler transform on the right hand side of the start symbol of a grammar.
In our case, this grammar (see the references) is based on the SA-IS algorithm, but it can be exchanged with other grammars that have locally consistent properties.

## Usage

### 🚀 Complete Test Run

A complete test run of our index can be done with the following lines:

```bash:
cd figiss/src
make test
./test 1 8 1024 0 8
```

### 🏎️ Play With Toy Example

A full set of construction, serialization, load and counting can be done with the following line:

```bash:
cd figiss/src
make
./figiss cs 4 ../data/text.txt ../data/text.figiss
./figiss lc ../data/text.figiss ../data/pattern.txt
```

The output should be `5`, which is the number of occurrences of `../data/pattern.txt` in `../data/text.txt`.

### ✔️ Prerequisites

To compile this project and run testing and experiments, you need the following tools:

 - `make`
 - a recent `g++` with C++17 support
 - [sdsl-lite](https://github.com/simongog/sdsl-lite)
 - `p7zip`
 - `xz`
 - 
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
cd figiss/src
make test
```

Experimental result presented in our papar is processed from the information generated by an executable call `experiment`
It can be compiled with:
```bash:
cd figiss/src
make experiment
```

### 🏗️ Index Construction

Syntax:

```bash:
./figiss cs λ "byte-text path" "index path"
```

 - construct index of (ASCII-encoded) text at `"byte-text path"`
	- with max factor size being `λ`, which can only be an integer in `[1..8]`
	- and serialize it to `"index path"`
	-`4` and `7` are preferred choices.  

Example:

```bash:
./figiss cs 4 ../data/text.txt ../data/text.figiss
```

Output:

```bash:
construct index of "/figiss/data/text.txt"
serialize index to "/figiss/data/text.figiss
```

### 🔎 Counting Query

```bash:
./figiss lc "index path" "byte-pattern path"
```

 - load index from `"index path"` and report number of occurences of (ACSII-encoded) pattern at `"pattern path"`
 - ⚠️ every byte in the file `"byte-pattern path"` are considered, even if the file ends with a new line character (which is commonly appended by text editors!)

Example:

```bash:
./figiss lc ../data/text.figiss ../data/pattern.txt
```

Output:

```bash:
load index from ".../figiss/data/text.figiss"
5
```

### Testing Counting Query

```bash:
make test
./test λ_l λ_r b p_l p_r
```

- test counting of figiss, where λ is from `λ_l` to `λ_r` and 1 <= `λ_l` <= `λ_r` <= 8
- for each p from `p_l` to `p_r`, pattern collection of `b` random patterns of length 2^p is generated 

### Re-run Experiments

```bash:
make experiment
./experiment λ_l λ_r b p_l p_r
```
- run experiments of figiss, where λ is from `λ_l` to `λ_r` and 1 <= `λ_l` <= `λ_r` <= 8
- for each p from `p_l` to `p_r`, pattern collection of `b` random patterns of length 2^p is generated 

## 📚 References

- Daniel Saad Nogueira Nunes, Felipe A. Louza, Simon Gog, Mauricio Ayala-Rincón, Gonzalo Navarro: [A Grammar Compression Algorithm Based on Induced Suffix Sorting. DCC 2018: 42-51](https://doi.org/10.1109/DCC.2018.00012)

## Using as a Library

You can also use FIGISS as a library in the following way:

### Sample Code for Constructing Index

```c++:
#include "index.h"
int main (int argc, char** argv)
{
  if (argc == 3)
  {
    auto byte_text_path {std::filesystem::path{argv[1]}};
    auto max_factor_size {static_cast<uint8_t>(std::stoull(argv[2]))}; // 1 <= max_factor_size <= 8
    figiss::Index index {byte_text_path, max_factor_size};
  }
  return 0;
}
```

### Sample Code for Serializing Index

```c++:
#include "index.h"
int main (int argc, char** argv)
{
  if (argc == 4)
  {
    auto byte_text_path {std::filesystem::path{argv[1]}};
    auto max_factor_size {static_cast<uint8_t>(std::stoull(argv[2]))}; // 1 <= max_factor_size <= 8
    auto index_path {std::filesystem::path{argv[3]}};
    figiss::Index index {byte_text_path, max_factor_size};
    index.Serialize(index_path);
  }
  return 0;
}
```

### Sample Code for Loading Index

```c++:
#include "index.h"
int main (int argc, char** argv)
{
  if (argc == 2)
  {
    auto index_path {std::filesystem::path{argv[1]}};
    figiss::Index index;
    index.Load(index_path);
  }
  return 0;
}
```

### Sample Code for Counting

```c++:
#include "index.h"
int main (int argc, char** argv)
{
  if (argc == 2)
  {
    auto index_path {std::filesystem::path{argv[1]}};
    figiss::Index index;
    index.Load(index_path);
    std::string pattern {"acgt"};
    std::cout << index.Count(std::begin(pattern), std::end(pattern)) << "\n";
  }
  return 0;
}
```
