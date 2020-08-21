CC									:= g++
WARNING_ERROR_FLAGS	:= -Wall -Werror -Wextra -Wpedantic -Wfatal-errors
CXXFLAGS						:= $(WARNING_ERROR_FLAGS) -std=c++17
CXXHEADERS					:= lib/sl_types.hpp lib/text_index_buckets.hpp lib/utility.hpp
OBJECTS							:= utility.o
EXECUTABLE					:= main

all: $(EXECUTABLE)

$(EXECUTABLE): main.cpp $(CXXHEADERS) $(OBJECTS)
	$(CC) $(CXXFLAGS) $^ -o $@

utility.o: lib/utility.cpp lib/utility.hpp
	$(CC) $(CXXFLAGS) -c $<

clean:
	rm -f $(EXECUTABLE) $(OBJECTS)

.PHONY: all clean
