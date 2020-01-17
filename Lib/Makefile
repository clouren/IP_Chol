# On the Mac, you need to first install Xcode.  Then
# type this command in the Terminal:
#       xcode-select --install
# to ensure that you have the command-line tools installed.
#
# To just compile the library, do:
#       make
#
# To run a demo using the library 
#       cd ../Demo ; make
#


LIBRARY = libipchol
CF = $(CFLAGS) $(CPPFLAGS) $(TARGET_ARCH) -O
        
I = -I../Include -I../Source -I../SLIP_LU-master/SuiteSparse_config -I../SLIP_LU-master/COLAMD/Include -I../SLIP_LU-master/AMD/Include -I../SLIP_LU-master/SLIP_LU/Include -I../SLIP_LU-master/SLIP_LU/Source 
RANLIB = ranlib
ARCHIVE = $(AR) $(ARFLAGS)
CP = cp -f

AR_TARGET = $(LIBRARY).a

all: install

CS = IP_array_div.o IP_array_mul.o IP_Chol_counts.o IP_Chol_ereach.o IP_Chol_etree.o \
     IP_Chol_leaf.o IP_Chol_ltsolve.o IP_Chol_permute_A.o IP_Chol_post.o \
     IP_Chol_tdfs.o IP_create_mpz_array2.o IP_cumsum_chol.o IP_forward_sub.o IP_get_column.o \
     IP_Left_Chol_Factor.o IP_Left_Chol_triangular_solve.o IP_Pre_Left_Factor.o  \
     IP_process_command_line.o IP_reset_int_array.o IP_reset_int_array2.o IP_reset_mpz_array.o \
     IP_Solve.o IP_sparse_alloc.o IP_sparse_alloc2.o IP_Sym_chol_free.o IP_transpose.o \
     IP_tripread_double.o IP_Up_Chol_Factor.o IP_Up_Chol_triangular_solve.o IP_Up_get_column.o \

$(CS): ../Include/IP-Chol.h Makefile

%.o: ../Source/%.c ../Include/IP-Chol.h 
	$(CC) $(CF) $(I) -c $<

static: $(AR_TARGET)

$(AR_TARGET): $(CS)
	$(ARCHIVE)  $@ $^
	- $(RANLIB) $@

# install SLIP LU in this directory
install: $(AR_TARGET)

# uninstall SLIP LU
uninstall: purge

clean:
	- $(RM) *.o

purge: distclean

distclean: clean
	- $(RM) *.a *.obj *.dll *.dylib *.so *.so.*

