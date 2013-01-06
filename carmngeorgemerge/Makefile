LDFLAGS = -lm 

MPICXX=mpic++

MPICXXFLAGS = -O3 -Wno-deprecated -Wall

all: cbemd andersen

%.o : %.cpp
	$(MPICXX) $(MPICXXFLAGS) -c $< 

cbemd : verlet.o force_calc.o initialize.o read_xml.o read_interaction.o system.o atom.o misc.o integrator.o interaction.o  
	$(MPICXX) -o $@ $^ $(GTESTFLAGS) $(MPICXXFLAGS) 

andersen : andersen.o force_calc.o initialize.o read_xml.o read_interaction.o system.o atom.o misc.o integrator.o interaction.o
	$(MPICXX) -o $@ $^ $(GTESTFLAGS) $(MPICXXFLAGS)

clean:
	$(RM) *.o
	$(RM) .depend
	$(RM) cbemd, andersen
