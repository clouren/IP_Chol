default: all

include ./SLIP_LU-master/SuiteSparse_config/SuiteSparse_config.mk


#LIB_DIRS=./SLIP_LU-master/lib/
PATH_TO_SLIP_LU=-Wl,-rpath,"./SLIP_LU-master/lib/"
CF = $(CFLAGS) $(CPPFLAGS) $(TARGET_ARCH) $(PATH_TO_SLIP_LU) -O
I = -I./Include -I./Source -I./SLIP_LU-master/SuiteSparse_config -I./SLIP_LU-master/COLAMD/Include -I./SLIP_LU-master/AMD/Include -I./SLIP_LU-master/SLIP_LU/Include -I./SLIP_LU-master/SLIP_LU/Source -I./SLIP_LU-master/SLIP_LU/Lib -I./SLIP_LU-master/Include/



LDLIBS += -lm -lgmp -lmpfr -lcolamd -lamd -L./SLIP_LU-master/lib/ -lsliplu
CS = ./Lib/libipchol.a  $(LDLIBS)


all: lib IP_Chol IP_Chol_debug
	- ./IP_Chol

lib:
	( cd ./Lib ; $(MAKE) )

IP_Chol_debug: lib IP_Chol.c Makefile
	$(CC) $(LDFLAGS) $(CF) $(I) -g -o IP_Chol_debug IP_Chol.c $(CS) 


IP_Chol: lib IP_Chol.c Makefile
	$(CC) $(LDFLAGS) $(CF) $(I) -o IP_Chol IP_Chol.c $(CS) 

clean:
	- $(RM) *.o

purge: distclean

distclean: clean
	- $(RM) -r IP_Chol IP_Chol_debug *.a *.dSYM *.obj *.dll
