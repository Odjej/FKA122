CC ?= gcc
CFLAGS ?= -O2 -Wall -Wextra -std=c99
CPPFLAGS += -I$(HOME)/.local/include
LDFLAGS += -L$(HOME)/.local/lib -Wl,-rpath,$(HOME)/.local/lib
LDLIBS += -lgsl -lgslcblas -lm

OBJECTS = main.o linalg.o gsl_random.o
PROGRAM = program

%.o: %.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c -o $@ $<

all: $(PROGRAM)

$(PROGRAM): $(OBJECTS)
	$(CC) $(LDFLAGS) -o $@ $^ $(LDLIBS)

clean:
	$(RM) $(OBJECTS) $(PROGRAM)

.PHONY: all clean
