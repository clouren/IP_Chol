default: all

include ./SuiteSparse_config/SuiteSparse_config.mk

# CF = $(CFLAGS) $(CPPFLAGS) $(TARGET_ARCH) -O
I = -I./Include -I./Source -I./SLIP_LU-master/SuiteSparse_config -I./SLIP_LU-master/COLAMD/Include -I./SLIP_LU-master/AMD/Include -I./SLIP_LU-master/SLIP_LU/Include -I./SLIP_LU-master/SLIP_LU/Source -I./SLIP_LU-master/SLIP_LU/Lib

# LDFLAGS = -L../../lib

LDLIBS += -lm -lgmp -lmpfr -lcolamd -lamd
CS = ./SLIP_LU-master/SLIP_LU/Lib/libslip.a ./Lib/libipchol.a $(LDLIBS)


all: lib SLIP_Chol Up_Chol SLIP_Chol_debug Up_Chol_debug
	- ./SLIP_Chol
	- ./Up_Chol

lib:
	( cd ./Lib ; $(MAKE) )

SLIP_Chol_debug: lib SLIP_Chol.c Makefile
	$(CC) $(LDFLAGS) $(CF) $(I) -g -o SLIP_Chol_debug SLIP_Chol.c $(CS)

SLIP_Chol: lib SLIP_Chol.c Makefile
	$(CC) $(LDFLAGS) $(CF) $(I) -o SLIP_Chol SLIP_Chol.c $(CS)

Up_Chol_debug: lib Up_Chol.c Makefile
	$(CC) $(LDFLAGS) $(CF) $(I) -g -o Up_Chol_debug Up_Chol.c $(CS)
	
Up_Chol: lib Up_Chol.c Makefile
	$(CC) $(LDFLAGS) $(CF) $(I) -o Up_Chol Up_Chol.c $(CS)

#lu: lib SLIPLU.c demos.h demos.c Makefile
#	$(CC) $(LDFLAGS) $(CF) $(I) -o SLIPLU SLIPLU.c demos.c $(CS)

clean:
	- $(RM) *.o

purge: distclean

distclean: clean
	- $(RM) -r SLIP_Chol Up_Chol *.a *.dSYM *.obj *.dll
