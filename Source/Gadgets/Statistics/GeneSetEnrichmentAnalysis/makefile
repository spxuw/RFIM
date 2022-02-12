############################################################
#            Makefile for GetHistogram.cpp                 #
#                     Yang-Yu Liu                          #
#                     08/25/13                             #
############################################################
# home directory
BASE = ./

INCLALL= -I$(BASE) -I../../Random

# executable file
EXE = gsea

# macro for name of my compiler: 
CC   = g++

# macro for compilation options:
CFLAGS = -O3

LIBS = -lm  -lgsl  #-lglpk -lgslcblas

# main file
MAIN = main.cpp
OBJMAIN = main.o

# tool files
TOOL1 = gsea.cpp
OBJTOOL1 = gsea.o

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


