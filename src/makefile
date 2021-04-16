CC 									:= g++
CXX_FLAGS 					:= -std=c++17 -O3
DEBUG_FLAGS					:= -DNDEBUG
WARNING_ERROR_FLAGS	:= -Wall -Wextra -Wpedantic -Werror -Wfatal-errors
LIBS 								:= -ldivsufsort -ldivsufsort64 -lm -lsdsl

all: main.cc
	${CC} ${CXX_FLAGS} ${DEBUG_FLAGS} ${WARNING_ERROR_FLAGS} $^ -o main ${LIBS}

clean:
	rm -f main

.PHONY: all clean
