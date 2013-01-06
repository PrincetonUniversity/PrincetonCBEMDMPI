LDFLAGS = -lm 

MPICXX=mpic++

PATHTOBOOST = /home/gkhoury/boost_1_52_0

MPICXXFLAGS = -O3 -Wno-deprecated -I $(PATHTOBOOST)

all: verlet andersen

%.o : %.cpp
	$(MPICXX) $(MPICXXFLAGS) -c $< 

verlet : verlet.o force_calc.o initialize.o read_xml.o read_interaction.o system.o atom.o misc.o integrator.o interaction.o  
	$(MPICXX) -o $@ $^ $(GTESTFLAGS) $(MPICXXFLAGS) 

andersen : andersen.o force_calc.o initialize.o read_xml.o read_interaction.o system.o atom.o misc.o integrator.o interaction.o
	$(MPICXX) -o $@ $^ $(GTESTFLAGS) $(MPICXXFLAGS)

clean:
	$(RM) *.o
	$(RM) .depend
	$(RM) verlet
	$(RM) andersen
