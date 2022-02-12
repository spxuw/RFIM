############################################################
#            Makefile for parser.cpp                       #
#                     Yang-Yu Liu                          #
#                     09/05/13                             #
############################################################
# home directory
LIBRARY = /Users/spyli/Gadgets/
INCLALL   =  -I$(LIBRARY)/Random \
             -I$(LIBRARY)/Timer \
             -I$(LIBRARY)/Find \
             -I$(LIBRARY)/Color \
             -I$(LIBRARY)/Directory 


# executable file
EXE = parser

# macro for name of my compiler: 
CC   = g++

# macro for compilation options:
CFLAGS = -O

LIBS = -lm  #-lgsl  #-lglpk -lgslcblas

# main file
MAIN = main.cpp
OBJMAIN = main.o

# tool files
TOOL1 = Parser.cpp
OBJTOOL1 = Parser.o

TOOL2 = $(LIBRARY)/Find/Find.cpp
OBJTOOL2 = Find.o

OBJS = $(OBJMAIN) $(OBJTOOL1) $(OBJTOOL2) 

$(EXE): $(OBJS) 
	@echo "linking ..."
	$(CC) $(CFLAGS) -o $(EXE) $(OBJS) $(LIBS)

$(OBJMAIN): $(MAIN)  
	@echo "compiling ..."
	$(CC) $(CFLAGS) -c $(MAIN) $(INCLALL)

$(OBJTOOL1): $(TOOL1)  
	@echo "compiling ..."
	$(CC) $(CFLAGS) -c $(TOOL1) $(INCLALL)

$(OBJTOOL2): $(TOOL2)  
	@echo "compiling ..."
	$(CC) $(CFLAGS) -c $(TOOL2) $(INCLALL)

clean:	
	rm -f $(OBJS) $(EXE)


