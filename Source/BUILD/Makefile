############################################################
#                 Makefile for RFIM                        #
############################################################
# home directory
#BASEBGL= /Users/xu-wenwang/Dropbox/Github/RFIM/RFIM/Source/boost_1_33_1
BASE   = ../
LIBRARY = /Users/xu-wenwang/Dropbox/Github/RFIM/RFIM/Source/Gadgets/

INCLALL   =  -I$(LIBRARY)/Random \
             -I$(LIBRARY)/Timer \
             -I$(LIBRARY)/Splitter \
	     -I$(LIBRARY)/SpecialFunctions \
             -I$(LIBRARY)/Find \
             -I$(LIBRARY)/Color \
             -I$(LIBRARY)/Statistics/Histogram \
	     -I$(BASE)/WINDOW \
             -I$(BASE)/STATE \
             -I$(BASE)/HIPR  \
             -I$(BASE)/MAIN \
             -I$(BASE)/TOOL \
             -I/usr/local/include 


############################################################
# main files
MAIN0 = $(BASE)/MAIN/main_preprocess.cpp 
OBJMAIN0 = main_preprocess.o 

MAIN4 = $(BASE)/MAIN/main_Network_Real.cpp 
OBJMAIN4 = main_Network_Real.o 

MAIN6 = $(BASE)/MAIN/main_permute.cpp 
OBJMAIN6 = main_permute.o 


# solver file
SOLVER = $(BASE)/HIPR/calculate_GS.cpp
OBJSOLVER = calculate_GS.o

# class file 1
#CLASS1 = $(BASE)/WINDOW/DropletAnalyse.cpp  
#OBJCLASS1 = DropletAnalyse.o 

# class file 2
CLASS2 = $(BASE)/STATE/State.cpp  
OBJCLASS2 = State.o 

CLASS2a = $(BASE)/STATE/State_constructor.cpp  
OBJCLASS2a = State_constructor.o 

CLASS2b = $(BASE)/STATE/State_initialize.cpp  
OBJCLASS2b = State_initialize.o 

CLASS2c = $(BASE)/STATE/State_setRandomField.cpp  
OBJCLASS2c = State_setRandomField.o 

CLASS2d = $(BASE)/STATE/State_energy.cpp  
OBJCLASS2d = State_energy.o 

CLASS2e = $(BASE)/STATE/State_tools.cpp  
OBJCLASS2e = State_tools.o 

CLASS2f = $(BASE)/STATE/State_mapping.cpp  
OBJCLASS2f = State_mapping.o 

CLASS2g = $(BASE)/STATE/State_advanced_mapping.cpp  
OBJCLASS2g = State_advanced_mapping.o 

CLASS2h = $(BASE)/STATE/State_avalanche.cpp  
OBJCLASS2h = State_avalanche.o 

CLASS2i = $(BASE)/STATE/State_Score.cpp  
OBJCLASS2i = State_Score.o 

CLASS2j = $(BASE)/STATE/State_Write.cpp  
OBJCLASS2j = State_Write.o 

CLASS2k = $(BASE)/STATE/State_EnrichmentAnalysis.cpp  
OBJCLASS2k = State_EnrichmentAnalysis.o 

CLASS2l = $(BASE)/STATE/State_DFS.cpp  
OBJCLASS2l = State_DFS.o 

CLASS2m = $(BASE)/STATE/State_SingleSpinFlipDynamics.cpp  
OBJCLASS2m = State_SingleSpinFlipDynamics.o 

CLASS2n = $(BASE)/STATE/State_jActive.cpp  
OBJCLASS2n = State_jActive.o 

CLASS2o = $(BASE)/STATE/State_dmGWAS.cpp  
OBJCLASS2o = State_dmGWAS.o 



OBJCLASSES = $(OBJCLASS1) $(OBJCLASS2) $(OBJCLASS2a) $(OBJCLASS2b) $(OBJCLASS2c) $(OBJCLASS2d) $(OBJCLASS2e) $(OBJCLASS2f) $(OBJCLASS2g) $(OBJCLASS2h) $(OBJCLASS2i) $(OBJCLASS2j) $(OBJCLASS2k) $(OBJCLASS2l) $(OBJCLASS2m) $(OBJCLASS2n) $(OBJCLASS2o)

# tool files 
TOOL0 = $(LIBRARY)/Find/find.cpp
OBJTOOL0 = find.o 

TOOL1 = $(LIBRARY)/SpecialFunctions/inversenormalCDF.cpp
OBJTOOL1 = inversenormalCDF.o 

TOOL2 = $(LIBRARY)/SpecialFunctions/hypergeometric.cpp
OBJTOOL2 = hypergeometric.o 

TOOL3 = $(LIBRARY)/Statistics/Histogram/Histogram.cpp
OBJTOOL3 = Histogram.o 

TOOL4 = $(LIBRARY)/Statistics/Histogram/Statistics.cpp
OBJTOOL4 = Statistics.o 

TOOL5 = $(BASE)/TOOL/MakeDirectory.cpp  
OBJTOOL5 = MakeDirectory.o 

OBJTOOLS = $(OBJTOOL0) $(OBJTOOL1) $(OBJTOOL2) $(OBJTOOL3) $(OBJTOOL4) $(OBJTOOL5)



OBJRECIPES = $(OBJRECIPE1)  $(OBJRECIPE2) 

# object file 
OBJ0 = $(OBJMAIN0) $(OBJTOOLS)
OBJ1 = $(OBJMAIN1) $(OBJSOLVER) $(OBJRECIPES) $(OBJCLASSES) $(OBJTOOLS)
OBJ2 = $(OBJMAIN2) $(OBJSOLVER) $(OBJRECIPES) $(OBJCLASSES) $(OBJTOOLS)
OBJ3 = $(OBJMAIN3) $(OBJSOLVER) $(OBJRECIPES) $(OBJCLASSES) $(OBJTOOLS)
OBJ4 = $(OBJMAIN4) $(OBJSOLVER) $(OBJRECIPES) $(OBJCLASSES) $(OBJTOOLS)
OBJ5 = $(OBJMAIN5) $(OBJSOLVER) $(OBJRECIPES) $(OBJCLASSES) $(OBJTOOLS)
OBJ6 = $(OBJMAIN6) $(OBJTOOLS) 

OBJS = $(OBJMAIN0) $(OBJMAIN1) $(OBJMAIN2) $(OBJMAIN3) $(OBJMAIN4) $(OBJMAIN5) $(OBJMAIN6) \
       $(OBJSOLVER) $(OBJRECIPES) $(OBJCLASSES) $(OBJTOOLS)


############################################################
# executable file
EXE0 = HMS_PreProcess
EXE4 = HMS_Network_Real
EXE6 = HMS_Permute

EXES = $(EXE0) $(EXE1) $(EXE2) $(EXE3) $(EXE4) $(EXE6)   
all : HMS_PreProcess HMS_Network_Real HMS_Permute


# macro for name of my compiler: 
CC   = g++ 

# compiler flags:

# PRINT_TIME         show the computing time, added by YL 05/29/06
#
# CHECK_SOLUTION     check if at the end the sink is unreachable 
#                    from the source
# CUT_ONLY           define to compute min-cut but not max-flow
#
# PRINT_STAT         define to print detailed operation statistics
# PRINT_FLOW         define to print out flow value
# PRINT_CUT          define to print out the sink size of the cut
# WAVE_INIT          wave-style initialization (initial global update,
#                    then a wave of flow pushing). Better on some problems
# OLD_INIT           define to use old-style initialization (large excess
#                    at the source, global update)
# INIT_UPDATE        perform initial global update
# EXCESS_TYPE_LONG   set excessType to long, otherwise it is long long int
#                    if "long long int" not supported, change to "double"
#                    (in types.h)

# macro for compilation options:

#CFLAG  = -g -Wall
#CFLAGS = -g -DPRINT_TIME -DPRINT_CUT -DCUT_ONLY -Wall

CFLAG  = -O3 -Wall -Wno-write-strings
CFLAGS = -O3 -DPRINT_TIME -DPRINT_CUT -DCUT_ONLY -Wall -Wno-write-strings

# YYL 04/16/2014
# -O4 sometimes cause linking error: ld: lto: could not merge in main.o because Invalid ALLOCA record for architecture x86_64

LIBS = -lm -L/usr/local/lib -lgsl

############################################################
$(EXE0): $(OBJ0) 
	@echo "linking .........................."
	$(CC) $(CFLAGS) $(LIBS) -o $(EXE0) $(OBJ0)

$(EXE1): $(OBJ1) 
	@echo "linking .........................."
	$(CC) $(CFLAGS) $(LIBS) -o $(EXE1) $(OBJ1)

$(EXE2): $(OBJ2) 
	@echo "linking .........................."
	$(CC) $(CFLAGS) $(LIBS) -o $(EXE2) $(OBJ2)

$(EXE3): $(OBJ3) 
	@echo "linking .........................."
	$(CC) $(CFLAGS) $(LIBS) -o $(EXE3) $(OBJ3)

$(EXE4): $(OBJ4) 
	@echo "linking .........................."
	$(CC) $(CFLAGS) $(LIBS) -o $(EXE4) $(OBJ4)
	@echo "Done !!! :)"

$(EXE5): $(OBJ5) 
	@echo "linking .........................."
	$(CC) $(CFLAGS) $(LIBS) -o $(EXE5) $(OBJ5)
	@echo "Done !!! :)"

$(EXE6): $(OBJ6) 
	@echo "linking .........................."
	$(CC) $(CFLAGS) $(LIBS) -o $(EXE6) $(OBJ6)
	@echo "Done !!! :)"

$(OBJMAIN0): $(MAIN0) 
	@echo "compling main file0 ..............."
	$(CC) $(CFLAG) -c $(MAIN0) $(INCLALL) 

$(OBJMAIN1): $(MAIN1) 
	@echo "compling main file1 ..............."
	$(CC) $(CFLAG) -c $(MAIN1) $(INCLALL) 

$(OBJMAIN2): $(MAIN2) 
	@echo "compling main file2 ..............."
	$(CC) $(CFLAG) -c $(MAIN2) $(INCLALL) 

$(OBJMAIN3): $(MAIN3) 
	@echo "compling main file3 ..............."
	$(CC) $(CFLAG) -c $(MAIN3) $(INCLALL) 

$(OBJMAIN4): $(MAIN4) 
	@echo "compling main file4 ..............."
	$(CC) $(CFLAG) -c $(MAIN4) $(INCLALL) 

$(OBJMAIN5): $(MAIN5) 
	@echo "compling main file5 ..............."
	$(CC) $(CFLAG) -c $(MAIN5) $(INCLALL) 

$(OBJMAIN6): $(MAIN6) 
	@echo "compling main file6 ..............."
	$(CC) $(CFLAG) -c $(MAIN6) $(INCLALL) 

$(OBJSOLVER): $(SOLVER) 
	@echo "compling solver file ..."
	$(CC) $(CFLAGS) -c $(SOLVER) $(INCLALL) 

$(OBJRECIPE1): $(RECIPE1) 
	@echo "compling numerical recipe file 1..."
	$(CC) $(CFLAG) -c $(RECIPE1) $(INCLALL)

$(OBJRECIPE2): $(RECIPE2) 
	@echo "compling numerical recipe file 2..."
	$(CC) $(CFLAG) -c $(RECIPE2) $(INCLALL)



$(OBJCLASS1): $(CLASS1) 
	@echo "compling class file 1..."
	$(CC) $(CFLAG) -c $(CLASS1) $(INCLALL)

$(OBJCLASS2): $(CLASS2) 
	@echo "compling class file 2..."
	$(CC) $(CFLAG) -c $(CLASS2) $(INCLALL)

$(OBJCLASS2a): $(CLASS2a)
	@echo "compling class file 2..."
	$(CC) $(CFLAG) -c $(CLASS2a) $(INCLALL)
$(OBJCLASS2b): $(CLASS2b)
	@echo "compling class file 2..."
	$(CC) $(CFLAG) -c $(CLASS2b) $(INCLALL)
$(OBJCLASS2c): $(CLASS2c)
	@echo "compling class file 2..."
	$(CC) $(CFLAG) -c $(CLASS2c) $(INCLALL)
$(OBJCLASS2d): $(CLASS2d)
	@echo "compling class file 2..."
	$(CC) $(CFLAG) -c $(CLASS2d) $(INCLALL)
$(OBJCLASS2e): $(CLASS2e)
	@echo "compling class file 2..."
	$(CC) $(CFLAG) -c $(CLASS2e) $(INCLALL)
$(OBJCLASS2f): $(CLASS2f)
	@echo "compling class file 2..."
	$(CC) $(CFLAG) -c $(CLASS2f) $(INCLALL)
$(OBJCLASS2g): $(CLASS2g)
	@echo "compling class file 2..."
	$(CC) $(CFLAG) -c $(CLASS2g) $(INCLALL)
$(OBJCLASS2h): $(CLASS2h)
	@echo "compling class file 2..."
	$(CC) $(CFLAG) -c $(CLASS2h) $(INCLALL)
$(OBJCLASS2i): $(CLASS2i)
	@echo "compling class file 2..."
	$(CC) $(CFLAG) -c $(CLASS2i) $(INCLALL)
$(OBJCLASS2j): $(CLASS2j)
	@echo "compling class file 2..."
	$(CC) $(CFLAG) -c $(CLASS2j) $(INCLALL)
$(OBJCLASS2k): $(CLASS2k)
	@echo "compling class file 2..."
	$(CC) $(CFLAG) -c $(CLASS2k) $(INCLALL)
$(OBJCLASS2l): $(CLASS2l)
	@echo "compling class file 2..."
	$(CC) $(CFLAG) -c $(CLASS2l) $(INCLALL)
$(OBJCLASS2m): $(CLASS2m)
	@echo "compling class file 2..."
	$(CC) $(CFLAG) -c $(CLASS2m) $(INCLALL)
$(OBJCLASS2n): $(CLASS2n)
	@echo "compling class file 2..."
	$(CC) $(CFLAG) -c $(CLASS2n) $(INCLALL)
$(OBJCLASS2o): $(CLASS2o)
	@echo "compling class file 2..."
	$(CC) $(CFLAG) -c $(CLASS2o) $(INCLALL)

$(OBJTOOL0): $(TOOL0) 
	@echo "compling tool file 0..."
	$(CC) $(CFLAG) -c $(TOOL0) $(INCLALL)

$(OBJTOOL1): $(TOOL1) 
	@echo "compling tool file 1..."
	$(CC) $(CFLAG) -c $(TOOL1) $(INCLALL)

$(OBJTOOL2): $(TOOL2) 
	@echo "compling tool file 2..."
	$(CC) $(CFLAG) -c $(TOOL2) $(INCLALL)

$(OBJTOOL3): $(TOOL3) 
	@echo "compling tool file 3..."
	$(CC) $(CFLAG) -c $(TOOL3) $(INCLALL)

$(OBJTOOL4): $(TOOL4) 
	@echo "compling tool file 4..."
	$(CC) $(CFLAG) -c $(TOOL4) $(INCLALL)

$(OBJTOOL5): $(TOOL5) 
	@echo "compling tool file 5..."
	$(CC) $(CFLAG) -c $(TOOL5) $(INCLALL)

clean:
	@echo "delete .o files .........."	
	rm -f $(OBJS) $(EXES)



