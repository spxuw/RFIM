############################################################
#            Makefile for test2.cpp                         #
#                     Yang Liu                             #
#                     09/07/04                             #
############################################################
# home directory
BASE = /home/yangliu/work/C++ing/NR/cpp

# directory for include files
INCLS = -I$(BASE)/other 

# head file
HEAD = $(BASE)/other/nr.h

############################################################
# main file
MAIN = xscrsho.cpp
OBJMAIN = xscrsho.o 

# numerical recipe file  
RECIPES = $(BASE)/recipes/scrsho.cpp \
          $(BASE)/recipes/bessj0.cpp \

OBJRECIPES = scrsho.o   \
             bessj0.o \
            
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


