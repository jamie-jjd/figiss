CC 									:= g++
CXX_FLAGS 					:= -std=c++17 -O3
INCLUDE_PATH				:= -isystem./include/others -I./include
LIB_PATH						:= -L./lib/others -L./lib
DEBUG_FLAGS					:= -DNDEBUG
WARNING_ERROR_FLAGS	:= -Wall -Wextra -Wpedantic -Werror -Wfatal-errors
LIBS 								:= -ldivsufsort -ldivsufsort64 -lm -lsdsl

all: main.cc
	${CC} ${CXX_FLAGS} ${INCLUDE_PATH} ${LIB_PATH} ${DEBUG_FLAGS} ${WARNING_ERROR_FLAGS} $^ -o main ${LIBS}

clean:
	rm -f main

.PHONY: all clean
