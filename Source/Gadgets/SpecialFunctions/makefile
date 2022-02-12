############################################################
#            Makefile for Fishertest                       #
#                     Yang-Yu Liu                          #
#                     04/25/14                             #
############################################################
# home directory
BASE = ./

INCLALL= -I$(BASE) -I/usr/local/include 


# executable file
EXE = fishertest

# macro for name of my compiler: 
CC   = g++

# macro for compilation options:
CFLAGS = -O3

#LIBS = -lm  -lgsl  #-lglpk -lgslcblas
LIBS = -lm -L/usr/local/lib -lgsl


# main file
MAIN = test_hypergeometric.cpp
OBJMAIN = test_hypergeometric.o

# tool files
TOOL1 = hypergeometric.cpp
OBJTOOL1 = hypergeometric.o

OBJS = $(OBJMAIN) $(OBJTOOL1)

$(EXE): $(OBJS) 
	@echo "linking ..."
	$(CC) $(CFLAGS) -o $(EXE) $(OBJS) $(LIBS)

$(OBJMAIN): $(MAIN)  
	@echo "compiling ..."
	$(CC) $(CFLAGS) -c $(MAIN) $(INCLALL)

$(OBJTOOL1): $(TOOL1)  
	@echo "compiling ..."
	$(CC) $(CFLAGS) -c $(TOOL1) $(INCLALL)


clean:	
	rm -f $(OBJS) $(EXE)


