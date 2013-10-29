CC = c++ -Wall -funroll-loops #-pg #-xHOST # -O3

PROG = gaussian
HDRS = MersenneTwister.h Simulation.h VecND.h
SRCS = Main.cpp Simulation.cpp VecND.cpp
OBJS = Main.o Simulation.o VecND.o

$(PROG): $(OBJS)
	$(CC) $(OBJS) -o $(PROG)
	
Simulation.o: Simulation.cpp $(HDRS)
	$(CC) -c Simulation.cpp -o Simulation.o
  
VecND.o: VecND.cpp $(HDRS)
	$(CC) -c VecND.cpp -o VecND.o

Main.o: Main.cpp $(HDRS)
	$(CC) -c Main.cpp -o Main.o 

clean:
	rm -f $(PROG) $(OBJS) 
