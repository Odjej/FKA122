CC = gcc

LOCAL = ${HOME}/.local
CFLAGS = \
	-pedantic \
	-Wall \
	-Werror \
	-g \
	-Iinclude \
	-MMD \
	-I${LOCAL}/include
CFLAGS_OPT = \
	-O2 \
	-march=native
LDFLAGS = \
	-Wl,-rpath,${LOCAL}/lib \
	-L${LOCAL}/lib
LIB = \
	-lm \
	-lgsl \
	-lgslcblas

OBJ = \
	obj/tools.o \
	obj/unit_conversion.o \
	obj/verlet.o
MAIN = obj/main.o

program: $(MAIN) $(OBJ)
	$(CC) $(LDFLAGS) $^ -o $@ $(LIB)

obj/%.o: src/%.c | obj
	$(CC) -c $(CFLAGS) $(CFLAGS_OPT) $< -o $@

obj:
	mkdir -p obj

-include $(OBJ:.o=.d) $(MAIN:.o=.d)

clean:
	rm -rf obj program

.PHONY: clean