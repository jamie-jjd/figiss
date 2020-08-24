CC									:= g++
WARNING_ERROR_FLAGS	:= -Wall -Werror -Wextra -Wpedantic #-Wfatal-errors
CXXFLAGS						:= $(WARNING_ERROR_FLAGS) -std=c++17
CXXHEADERS					:= lib/sl_types.hpp lib/text_index_buckets.hpp lib/trie.hpp lib/utility.hpp
OBJECTS							:=
EXECUTABLE					:= main

all: $(EXECUTABLE)

$(EXECUTABLE): main.cpp $(CXXHEADERS) $(OBJECTS)
	$(CC) $(CXXFLAGS) $^ -o $@

clean:
	rm -f $(EXECUTABLE) $(OBJECTS)

.PHONY: all clean
