cflag = -O3  
#lib   = -I$(HOME)/VF/include
#lflag = -L$(HOME)/VF/lib -lvf -lstdc++ -lm

com   = common.h
uts   = occur.h
cc    = g++

all:  FFSM queryDB
FFSM    :  sym.o gBase.o adj_matrix.o Ullman.o pattern.o induPatt.o \
          myTimer.o register.o  $(com) 
	$(cc) sym.o gBase.o myTimer.o pattern.o Ullman.o induPatt.o\
        adj_matrix.o register.o  -o sym $(cflag) $(lflag)
getBF   : getBF.o gBase.o adj_matrix.o Ullman.o common.h
	$(cc) getBF.o Ullman.o  gBase.o adj_matrix.o \
        -o getBF $(cflag)  $(lflag)
queryDB   : queryDB.o gBase.o adj_matrix.o Ullman.o common.h
	$(cc) queryDB.o Ullman.o  gBase.o adj_matrix.o \
        -o queryDB $(cflag)  $(lflag)
sortFE    : sortFE.o gBase.o adj_matrix.o Ullman.o common.h
	$(cc) sortFE.o gBase.o adj_matrix.o Ullman.o -o sortFE $(cflag)  $(lflag)
register.o: register.cpp register.h ${com}
	$(cc)  -c register.cpp ${cflag} 
myTimer.o : myTimer.h myTimer.cpp   ${com}
	$(cc) -c myTimer.cpp ${cflag}
adj_matrix.o : adj_matrix.cpp cost.h adj_matrix.h attribute.h ${com}
	$(cc) -c adj_matrix.cpp $(lib) $(cflag)
gBase.o      : gBase.cpp gBase.h occur.h ${com} ${uts}
	$(cc) -c gBase.cpp $(cflag) $(lib)
pattern.o : pattern.cpp  ${com} ${uts}
	$(cc) -c pattern.cpp $(cflag) $(lib)
induPatt.o : induPatt.cpp induPatt.h
	$(cc) -c induPatt.cpp $(cflag) $(lib)
sym.o     :  sym.cpp ${com} ${uts}  
	$(cc) -c sym.cpp $(cflag) $(lib)
Ullman.o   : Ullman.h Ullman.cpp ${com}
	$(cc) -c Ullman.cpp $(cflag) $(lib)
getBF.o   : getBF.cpp ${com} ${uts}
	$(cc) -c getBF.cpp ${cflag} ${lib}
queryDB.o :  queryDB.cpp ${com} ${uts}
	$(cc) -c queryDB.cpp ${cflag} ${lib}
sortFE.o  :  sortFE.cpp 
	$(cc) -c sortFE.cpp ${cflag} ${lib}
clean:  
	rm *.o ;  rm sym
