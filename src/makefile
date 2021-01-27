CC 								:= g++
WARN_ERROR_FLAGS	:= -Wall -Wextra -Wpedantic -Werror -Wfatal-errors
CXX_FLAGS 				:= -O3 -std=c++17
LIBS 							:= -ldivsufsort -ldivsufsort64 -lsdsl -ltudocomp_stat
LOG_FLAGS 				:= -DLOG
# LOG_FLAGS					+= -DLOG_GCI_CONSTRUCT -DLOG_FMI_CONSTRUCT
# LOG_FLAGS					+= -DLOG_GCI_SERIALIZE -DLOG_FMI_SERIALIZE
LOG_FLAGS					+= -DLOG_GCI_LOAD -DLOG_FMI_LOAD
# LOG_FLAGS					+= -DLOG_GCI_COUNT -DLOG_FMI_COUNT
LOG_FLAGS					+= -DLOG_VERBOSE_GCI_COUNT
DEBUG_FLAGS				:= -DNDEBUG
# DEBUG_FLAGS				+= -DDEBUG_GCI_COUNT

all: main.cpp
	$(CC) $(DEBUG_FLAGS) $(LOG_FLAGS) $(WARN_ERROR_FLAGS) $(CXX_FLAGS) $^ -o main $(LIBS)

clean:
	rm -f main

.PHONY: all clean
