CC 									:= g++
CXX_FLAGS 					:= -std=c++17 -O3
INCLUDE_PATH				:= -isystem./include/others -I./include
LIB_PATH						:= -L./lib/others -L./lib
DEBUG_FLAGS					:= -DNDEBUG
PARALLEL_FLAGS			:= -D_GLIBCXX_PARALLEL -fopenmp
WARNING_ERROR_FLAGS	:= -Wall -Wextra -Wpedantic -Werror -Wfatal-errors
LIBS 								:= -ldivsufsort -ldivsufsort64 -lm -lrlcsa -lsdsl -lstdc++fs -ltudocomp_stat

all: main.cc
	${CC} ${CXX_FLAGS} ${INCLUDE_PATH} ${LIB_PATH} ${DEBUG_FLAGS} ${WARNING_ERROR_FLAGS} $^ -o figiss ${LIBS}

test: test.cc
	${CC} ${CXX_FLAGS} ${INCLUDE_PATH} ${LIB_PATH} ${DEBUG_FLAGS} ${WARNING_ERROR_FLAGS} $^ -o test ${LIBS}

experiment: experiment.cc
	${CC} ${CXX_FLAGS} ${INCLUDE_PATH} ${LIB_PATH} ${DEBUG_FLAGS} ${PARALLEL_FLAGS} ${WARNING_ERROR_FLAGS} $^ -o experiment ${LIBS}

clean:
	rm -f figiss test experiment
	# rm -rf ../data/corpus ../data/parameter ../data/construction ../data/index ../data/index_space  ../data/pattern ../data/counting_time
.PHONY: all clean
