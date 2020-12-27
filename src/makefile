CC 								:= g++
WARN_ERROR_FLAGS	:= -Wall -Wextra -Wpedantic -Werror -Wfatal-errors
CXX_FLAGS 				:= -O3 -std=c++17
LIBS 							:= -ldivsufsort -ldivsufsort64 -lsdsl

all: main.cpp
	$(CC) -DNDEBUG $(WARN_ERROR_FLAGS) $(CXX_FLAGS) $^ -o main $(LIBS)

clean:
	rm -f main

.PHONY: all clean
