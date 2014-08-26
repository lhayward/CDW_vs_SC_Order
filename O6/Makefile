CC = c++ -Wall -funroll-loops #-pg #-xHOST # -O3

PROG = o6
HDRS = MersenneTwister.h FileReading.h Hypercube.h Lattice.h Measure.h Model.h SimParameters.h
SRCS = Main.cpp FileReading.cpp Hypercube.cpp Lattice.cpp Measure.cpp Model.cpp SimParameters.cpp
OBJS = Main.o FileReading.o Hypercube.o Lattice.o Measure.o Model.o SimParameters.o

$(PROG): $(OBJS)
	$(CC) $(OBJS) -o $(PROG)
	
FileReading.o: FileReading.cpp $(HDRS)
	$(CC) -c FileReading.cpp -o FileReading.o
	
Hypercube.o: Hypercube.cpp $(HDRS)
	$(CC) -c Hypercube.cpp -o Hypercube.o
	
Lattice.o: Lattice.cpp $(HDRS)
	$(CC) -c Lattice.cpp -o Lattice.o
	
Measure.o: Measure.cpp $(HDRS)
	$(CC) -c Measure.cpp -o Measure.o
	
Model.o: Model.cpp $(HDRS)
	$(CC) -c Model.cpp -o Model.o
	
SimParameters.o: SimParameters.cpp $(HDRS)
	$(CC) -c SimParameters.cpp -o SimParameters.o
	
Main.o: Main.cpp $(HDRS)
	$(CC) -c Main.cpp -o Main.o 
	
clean:
	rm -f $(PROG) $(OBJS) 