LDFLAGS = -lm 

MPICXX=mpic++

MPICXXFLAGS = -O3 -g -Wall -pedantic

all: cbemd

%.o : %.cpp
	$(MPICXX) $(MPICXXFLAGS) -c $< 

cbemd : main_mpi.o force_calc.o initialize.o read_xml.o read_interaction.o system.o atom.o misc.o domain_decomp.o integrator.o interaction.o  
	$(MPICXX) -o $@ $^ $(GTESTFLAGS) $(MPICXXFLAGS) 

clean:
	$(RM) *.o
	$(RM) .depend
	$(RM) cbemd
