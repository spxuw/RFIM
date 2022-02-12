############################################################
#            Makefile for test.cpp                         #
#                     Yang Liu                             #
#                     09/05/04                             #
############################################################
# home directory
BASE = /home/yangliu/work/C++ing/NR/cpp

# directory for include files
INCLS = -I$(BASE)/other 

# head file
HEAD = $(BASE)/other/nr.h

############################################################
# main file
MAIN = test.cpp
OBJMAIN = test.o 

# numerical recipe file  
RECIPES = $(BASE)/recipes/airy.cpp    \
          $(BASE)/recipes/bessik.cpp  \
          $(BASE)/recipes/bessjy.cpp  \
          $(BASE)/recipes/beschb.cpp  \
          $(BASE)/recipes/chebev.cpp  \

OBJRECIPES = airy.o   \
             bessik.o \
             bessjy.o \
             beschb.o \
             chebev.o \
############################################################

# object file 
OBJS = $(OBJMAIN) $(OBJRECIPES)

# executable file
EXE = run

# macro for name of my compiler: 
CC   = g++

# macro for compilation options:
CFLAG = -O

$(EXE): $(OBJS) 
	@echo "linking ..."
	$(CC) $(CFLAG) -o $(EXE) $(OBJS)

$(OBJMAIN): $(MAIN) $(HEAD) 
	$(CC) $(CFLAG) -c $*.cpp $(INCLS)

$(OBJRECIPES): $(RECIPES) $(HEAD)
	$(CC) $(CFLAG) -c $(RECIPES) $(INCLS)

clean:	
	rm -f $(OBJS) $(EXE)


