############################################################
#            Makefile for GetHistogram.cpp                 #
#                     Yang-Yu Liu                          #
#                     08/25/13                             #
############################################################
# home directory
BASE = ./

INCLALL= -I$(BASE) -I/usr/local/include 


# executable file
EXE = histogram

# macro for name of my compiler: 
CC   = g++

# macro for compilation options:
CFLAGS = -O3

#LIBS = -lm  -lgsl  #-lglpk -lgslcblas
LIBS = -lm  -L/usr/local/lib -lgsl -stdlib=libstdc++ #-lgslcblas  


# main file
MAIN = main.cpp
OBJMAIN = main.o

# tool files
TOOL1 = Histogram.cpp
OBJTOOL1 = Histogram.o

TOOL2 = Statistics.cpp
OBJTOOL2 = Statistics.o

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


