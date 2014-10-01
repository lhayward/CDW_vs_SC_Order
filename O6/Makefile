CC = c++ 
CFLAGS = -std=c++0x -Wall -funroll-loops -O3 #-pg #-xHOST

PROG = o6
HDRS = MersenneTwister.h FileReading.h Hyperrectangle.h Lattice.h Measure.h Model.h O6_Model.h SimParameters.h Vector_NDim.h VectorSpins.h
SRCS = Main.cpp FileReading.cpp Hyperrectangle.cpp Lattice.cpp Measure.cpp Model.cpp O6_Model.cpp SimParameters.cpp Vector_NDim.cpp VectorSpins.cpp
OBJS = Main.o FileReading.o Hyperrectangle.o Lattice.o Measure.o Model.o O6_Model.o SimParameters.o Vector_NDim.o VectorSpins.o

$(PROG): $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o $(PROG)
	
FileReading.o: FileReading.cpp $(HDRS)
	$(CC) $(CFLAGS) -c FileReading.cpp -o FileReading.o
	
Hyperrectangle.o: Hyperrectangle.cpp $(HDRS)
	$(CC) $(CFLAGS) -c Hyperrectangle.cpp -o Hyperrectangle.o
	
Lattice.o: Lattice.cpp $(HDRS)
	$(CC) $(CFLAGS) -c Lattice.cpp -o Lattice.o
	
Measure.o: Measure.cpp $(HDRS)
	$(CC) $(CFLAGS) -c Measure.cpp -o Measure.o
	
Model.o: Model.cpp $(HDRS)
	$(CC) $(CFLAGS) -c Model.cpp -o Model.o
	
O6_Model.o: O6_Model.cpp $(HDRS)
	$(CC) $(CFLAGS) -c O6_Model.cpp -o O6_Model.o
	
SimParameters.o: SimParameters.cpp $(HDRS)
	$(CC) $(CFLAGS) -c SimParameters.cpp -o SimParameters.o
	
Vector_NDim.o: Vector_NDim.cpp $(HDRS)
	$(CC) $(CFLAGS) -c Vector_NDim.cpp -o Vector_NDim.o
	
VectorSpins.o: VectorSpins.cpp $(HDRS)
	$(CC) $(CFLAGS) -c VectorSpins.cpp -o VectorSpins.o
	
Main.o: Main.cpp $(HDRS)
	$(CC) $(CFLAGS) -c Main.cpp -o Main.o 
	
clean:
	rm -f $(PROG) $(OBJS) 
