LDFLAGS = -lm 

MPICXX=mpic++

MPICXXFLAGS = -O3 -Wno-deprecated -I /home/gkhoury/boost_1_52_0

all: cbemd andersen

%.o : %.cpp
	$(MPICXX) $(MPICXXFLAGS) -c $< 

cbemd : main_mpi.o force_calc.o initialize.o read_xml.o read_interaction.o system.o atom.o misc.o integrator.o interaction.o  
	$(MPICXX) -o $@ $^ $(GTESTFLAGS) $(MPICXXFLAGS) 

andersen : main_anderson_mpi.o force_calc.o initialize.o read_xml.o read_interaction.o system.o atom.o misc.o integrator.o interaction.o
	$(MPICXX) -o $@ $^ $(GTESTFLAGS) $(MPICXXFLAGS)

clean:
	$(RM) *.o
	$(RM) .depend
	$(RM) cbemd
