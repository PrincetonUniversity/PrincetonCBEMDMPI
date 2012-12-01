atom = atom.o

CXXFLAGS = -Wall -O3 -lstdc++ -limf -Wno-deprecated

LDFLAGS = -lm 

MPICPP=/usr/local/openmpi/1.4.3/intel111/x86_64/bin/mpic++
OMPFLAGS = -fopenmp

MPICXXFLAGS = -Wall -O3 -lstdc++ -limf -Wno-deprecated
all: atom

%.o : %.c
	$(MPICPP) $(CXXFLAGS) -c $<

atom.o : atom.cpp
	$(MPICPP) $(CXXFLAGS) -c atom.cpp

atom: $(atom)
	$(MPICPP) $(MPICXXFLAGS) -o $@ $(LDFLAGS) $(CXXFLAGS) $^

clean:
	$(RM) *.o
	$(RM) .depend
