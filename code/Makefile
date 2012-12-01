atom = atom.o
bond = bond.o
misc = misc.o

# pull
# git merge origin/master
# push
CXXFLAGS = -Wall -O3 -lstdc++ -limf -Wno-deprecated

LDFLAGS = -lm 

MPICPP=/usr/local/openmpi/1.4.3/intel111/x86_64/bin/mpic++
OMPFLAGS = -fopenmp

MPICXXFLAGS = -Wall -O3 -lstdc++ -limf -Wno-deprecated
all: atom bond misc

%.o : %.c
	$(MPICPP) $(CXXFLAGS) -c $<

atom.o : atom.cpp
	$(MPICPP) $(CXXFLAGS) -c atom.cpp

atom: $(atom)
	$(MPICPP) $(MPICXXFLAGS) -o $@ $(LDFLAGS) $(CXXFLAGS) $^

bond.o : bond.cpp
	$(CXX) $(CXXFLAGS) -c $<

bond: $(bond)
	$(CXX) $(CXXFLAGS) -o $@ $(LDFLAGS) $(CXXFLAGS) $^

misc.o : misc.cpp
	$(CXX) $(CXXFLAGS) -c $<

misc: $(misc)
	$(CXX) $(CXXFLAGS) -o $@ $(LDFLAGS) $(CXXFLAGS) $^

clean:
	$(RM) *.o
	$(RM) .depend
