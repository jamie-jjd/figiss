CC									:= g++
WARNING_ERROR_FLAGS	:= -Wall -Wfatal-errors #-Wextra -Wpedantic -Werror
CXX_FLAGS						:= -O3 -std=c++17
LIBRARIES						:= -ldivsufsort -ldivsufsort64 -lsdsl

all: main.cpp
	$(CC) -DNDEBUG $(WARNING_ERROR_FLAGS) $(CPP_FLAGS) $(CXX_FLAGS) $^ -o main $(LIBRARIES)

test: test.cpp
	$(CC) -DNDEBUG $(WARNING_ERROR_FLAGS) $(CXX_FLAGS) $^ -o $@ $(LIBRARIES)

clean:
	rm -f main test

.PHONY: all clean
