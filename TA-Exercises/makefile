
CC = gcc
CFLAGS = -Wall -Wextra -O2 -std=c11
TOOLS_DIR = ../../C-introduction/C4/library
CPPFLAGS = -I$(TOOLS_DIR)/include
LDFLAGS = -lgsl -lgslcblas -lm

TARGET = simulation
SRC = main.c verlet.c unit_conversion.c
OBJ = $(SRC:.c=.o) tools.o

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CC) $(OBJ) -o $(TARGET) $(LDFLAGS)

%.o: %.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

tools.o: $(TOOLS_DIR)/src/tools.c $(TOOLS_DIR)/include/tools.h
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

run: $(TARGET)
	./$(TARGET)

clean:
	rm -f $(OBJ) $(TARGET)

.PHONY: all run clean
