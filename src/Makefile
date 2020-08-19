CC									:= g++
WARNING_ERROR_FLAGS	:= -Wall -Werror -Wextra -Wpedantic #-Wfatal-errors
CXXFLAGS						:= $(WARNING_ERROR_FLAGS) -std=c++17
OBJECTS							:= main.o
EXECUTABLE					:= main

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(CFLAGS) $(CXXFLAGS) $(OBJECTS) -o $@

clean:
	rm -f $(OBJECTS) main

.PHONY: all clean
