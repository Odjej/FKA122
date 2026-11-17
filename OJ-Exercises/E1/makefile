CC = gcc

LOCAL = ${HOME}/.local
CFLAGS = \
	 -pedantic \
	 -Wall \
	 -Werror \
	 -g \
	 -Iinclude \
	 -MMD \
	 -Wl,-rpath,${LOCAL}/lib \
	 -L${LOCAL}/lib \
	 -I${LOCAL}/include
CFLAGS_OPT = \
	     -O2 \
	     -march=native

LIB = \
      -lm \
      -lgsl \
      -lgslcblas
OBJ = \
      obj/tools.o\
	  obj/unit_conversion.o\
	  obj/verlet.o
MAIN = \
       obj/main.o

program: obj _program

-include $(OBJ:.o=.d) $(MAIN:.o=.d)
_program: $(MAIN) $(OBJ)
	$(CC) $(CFLAGS) $^ -o program $(LIB)

obj/%.o: src/%.c
	$(CC) -c $(CFLAGS) $< -o $@ 

obj:
	mkdir -p obj

clean:
	find -iname "*.o" -exec rm {} \;
	find -iname "*.d" -exec rm {} \;


.PHONY: clean