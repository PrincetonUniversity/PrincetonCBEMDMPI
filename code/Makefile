atom = atom.o
bond = bond.o
misc = misc.o
system = system.o
main = main.o

# pull
# git merge origin/master
# push
CXXFLAGS = -Wall -O3 -lstdc++ -limf -Wno-deprecated -std=c++0x

LDFLAGS = -lm 

#MPICPP=/usr/local/openmpi/1.4.3/intel111/x86_64/bin/mpic++
#MPICC=/usr/local/openmpi/1.4.3/intel111/x86_64/bin/mpicc
MPICC=mpicc
MPICPP=mpic++
OMPFLAGS = -fopenmp

MPICXXFLAGS = -Wall -O3 -lstdc++ -limf -Wno-deprecated -std=c++0x
all: driver

%.o : %.cpp
	$(MPICPP) $(MPICXXFLAGS) -c $< -o $*.o

#atom.o : atom.cpp
#	$(MPICC) $(CXXFLAGS) -c atom.cpp

#atom: $(atom)
#	$(MPICC) $(MPICXXFLAGS) -o $@ $(LDFLAGS) $(CXXFLAGS) $^

#bond.o : bond.cpp
#	$(MPICC) $(CXXFLAGS) -c $<

#bond: $(bond)
#	$(MPICC) $(CXXFLAGS) -o $@ $(LDFLAGS) $(CXXFLAGS) $^

#misc.o : misc.cpp
#	$(MPICC) $(CXXFLAGS) -c $<

#misc: $(misc)
#	$(MPICC) $(CXXFLAGS) -o $@ $(LDFLAGS) $(CXXFLAGS) $^


#main.o : main.cpp 
#	$(MPICPP) $(CXXFLAGS) -c $< 

driver : $(main) $(atom) $(bond) $(misc) $(system)
	 $(MPICPP) $(MPICXXFLAGS) -o $@ $(LDFLAGS) $(CXXFLAGS) $^

clean:
	$(RM) *.o
	$(RM) .depend
