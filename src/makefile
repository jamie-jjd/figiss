CC									:= g++
WARNING_ERROR_FLAGS	:= -Wall -Wextra -Wpedantic -Werror -Wfatal-errors
CXX_FLAGS						:= -g -O3 -std=c++17 -I /home/sada/sdsl/include/ -I /home/sada/sdsl/external/libdivsufsort/include/ -L /home/sada/sdsl/lib/ -L /home/sada/sdsl/external/libdivsufsort/lib/
LIBRARIES						:= -ldivsufsort -ldivsufsort64 -lsdsl

all: main.cpp
	$(CC) -DNDEBUG $(WARNING_ERROR_FLAGS) $(CXX_FLAGS) $^ -o main $(LIBRARIES)

test: test.cpp
	$(CC) -DNDEBUG $(WARNING_ERROR_FLAGS) $(CXX_FLAGS) $^ -o test $(LIBRARIES)

clean:
	rm -f main test

.PHONY: all clean
