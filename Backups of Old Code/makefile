#-------------------------------------------------------------------
# Compiler selection and flags
#--------------------------------------------------------------------
CC  = gcc
CCC = g++
debug = -g -o
optimize = -O2 -o


#--------------------------------------------------------------------
#Library flags
#--------------------------------------------------------------------
CCLSTD = -std=c++11
CCLGMP =   -lgmp -lgmpxx
CCLMPFR = -lmpfr
CCLNCOLS = -lcolamd -lamd
CCLNUMF = -lumfpack
#
#

#--------------------------------------------------------------------
# Applications
#--------------------------------------------------------------------- 
APPS = SLIP_Chol Up_Chol
Debugs = debugChol debugUp

all: $(APPS) $(Debugs)
debug: $(Debugs)
	
SLIP_Chol: SLIP_Chol.cpp
	$(CCC) SLIP_Chol.cpp $(optimize) SLIP_Chol.exe $(CCLGMP) $(CCLMPFR) $(CCLNCOLS) $(CCLSTD) $(CCLNUMF)

debugChol: SLIP_Chol.cpp
	$(CCC) SLIP_Chol.cpp $(debug) debug_SLIP_Chol.exe $(CCLGMP) $(CCLMPFR) $(CCLNCOLS) $(CCLSTD) $(CCLNUMF)
	
Up_Chol: UP_Chol.cpp
	$(CCC) UP_Chol.cpp $(optimize) Up_Chol.exe $(CCLGMP) $(CCLMPFR) $(CCLNCOLS) $(CCLSTD) $(CCLNUMF)

debugUp: UP_Chol.cpp
	$(CCC) UP_Chol.cpp $(debug) debug_Up_Chol.exe $(CCLGMP) $(CCLMPFR) $(CCLNCOLS) $(CCLSTD) $(CCLNUMF)

#SLIP_LU_mg: SLIP_LU_mg.cpp
#	$(CCC) SLIP_LU_mg.cpp $(optimize) SLIP_LU_mg.exe $(CCLGMP) $(CCLMPFR) $(CCLNCOLS) $(CCLSTD) $(CCLNUMF)

