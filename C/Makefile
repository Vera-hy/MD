SRC=MD.c control.c #util.c
OBJ=$(SRC:.c=.o)
CC=icc 
CFLAGS= -restrict  -Ofast -xCORE-AVX2 -no-vec
#-g -O0 -check=uninit -check-pointers:rw -qopenmp-simd -pg -qopt-report=5


all: MD

MD: $(OBJ)
	$(CC) $(CFLAGS) -o $@  $(OBJ) -lm


output.dat: MD input.dat
	./MD


clean:
	rm -f MD $(OBJ) output.dat*

$(OBJ) : coord.h Makefile


