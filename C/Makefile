SRC=MD.c control.c util.c 
OBJ=$(SRC:.c=.o)
CC=icc 
CFLAGS= #-g -O0 -check=uninit -check-pointers:rw -no-vec
Profile_FLAG= #-pg
OPTIMIZATION_FLAG= -Ofast
ISA_FLAG= -xBROADWELL


all: MD

MD: $(OBJ)
	$(CC) $(CFLAGS) $(Profile_FLAG) $(OPTIMIZATION_FLAG) $(ISA_FLAG) -o $@  $(OBJ) -lm


output.dat: MD input.dat
	./MD


clean:
	rm -f MD $(OBJ) output.dat*

$(OBJ) : coord.h Makefile


