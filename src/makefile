CC									:= g++
WARNING_ERROR_FLAGS	:= -Wall -Werror -Wextra -Wpedantic -Wfatal-errors
CXXFLAGS						:= $(WARNING_ERROR_FLAGS) -std=c++17

all:
	$(CC) $(CXXFLAGS) -I$(CURDIR)/include/ main.cpp -o main

test:
	$(CC) $(CXXFLAGS) -I$(CURDIR)/include/ test.cpp -o test

clean:
	rm -f main test

.PHONY: all clean
