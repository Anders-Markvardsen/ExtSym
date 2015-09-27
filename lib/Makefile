CC=gcc
#CFLAGS=-O6 -c -w
CFLAGS=-c -static -w
#FFLAGS=-O2 -c -w
FFLAGS=-c -w

OBJ=nrutil.o\
	jacobi.o\
	gasdev.o ran1.o\
	ffactrl.o\
	inverse.o choldc.o\
	ludcmp.o lubksb.o inverse_jac.o\
	erfcc.o \
	log_erfcc.o parcylu_my.o

.c.o :
	$(CC) $(CFLAGS) $<
libnr.a : $(OBJ)
	ar rv libnr.a $(OBJ)
clean :
	rm -f *.o *~
