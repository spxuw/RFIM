############################################################
#        Makefile for CheckReloadedOperators_ly.cpp        #
#                     Yang Liu                             #
#                     10/23/04                             #
############################################################
# home directory
BASE = /home/yangliu/Tools/C++/Library/NR/cpp/

# directory for include files
INCLS = -I$(BASE)/other 

# head file
HEAD = $(BASE)/other/nr.h

############################################################
# main file
#MAIN = CheckReloadedOperators_ly.cpp
#OBJMAIN = CheckReloadedOperators_ly.o 

MAIN = xludcmp.cpp
OBJMAIN = xludcmp.o 

# numerical recipe file  
RECIPES = $(BASE)/recipes/ludcmp.cpp 
OBJRECIPES = ludcmp.o 


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
	$(CC) $(CFLAG) -c $(MAIN) $(INCLS)

$(OBJRECIPES): $(RECIPES) $(HEAD)
	$(CC) $(CFLAG) -c $(RECIPES) $(INCLS)

clean:	
	rm -f $(OBJS) $(EXE)


