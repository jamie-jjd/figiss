CC									:= g++
WARNING_ERROR_FLAGS	:= -Wall -Werror -Wextra -Wpedantic #-Wfatal-errors
CXXFLAGS						:= $(WARNING_ERROR_FLAGS) -std=c++17
OBJECTS							:= utility.o main.o
EXECUTABLE					:= main

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(CXXFLAGS) $? -o $@

main.o: main.cpp lib/utility.hpp
	$(CC) $(CXXFLAGS) -c $<

utility.o: lib/utility.cpp lib/utility.hpp
	$(CC) $(CXXFLAGS) -c $<

clean:
	rm -f $(EXECUTABLE) $(OBJECTS)

.PHONY: all clean
