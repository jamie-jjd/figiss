### Compilation

```bash:
make
```

This would generete executable `gciis`.

### Execution

```bash:
# construct index of text at [text path] and serialize it to [index path]
# k can only be integer in [1, 8]
./gciis cs [k] [text path] [index path]
```

```bash:
# load index from [index path] and report number of occurences of pattern at [pattern path]
# k must be matched with the k of index
./gciis lc [k] [text path] [index path]  
```

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
