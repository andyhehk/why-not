CC      = g++
CFLAGS  = -g 

SRC = Whynot.cpp DominantGraph.cpp Structure.cpp config.cpp 
HDR = Whynot.h DominantGraph.h Structure.h config.h
OBJ = Whynot.o DominantGraph.o Structure.o config.o

#MATLAB_INCLUDE_DIR=/usr/local/MATLAB/R2010b/extern/include
#MATLAB_LIB_DIR=/usr/local/MATLAB/R2010b/extern/lib

#INCLUDE_DIR=./lib

all: whynot

whynot: $(OBJ)
	$(CC) -o $@ $^  

#OBJ: $(SRC) $(HDR)
#	$(CC) -c -I$(MATLAB_INCLUDE_DIR) $(CFLAGS) $<

Whynot.o : Whynot.cpp Whynot.h DominantGraph.h Structure.h 
	$(CC) -c $(CFLAGS) Whynot.cpp	

#libqp_splx.o : libqp_splx.c libqp.h
#	gcc -c libqp_splx.c

DominantGraph.o : DominantGraph.cpp DominantGraph.h Structure.h
	$(CC) -c $(CFLAGS) DominantGraph.cpp

Structure.o : Structure.cpp Structure.h
	$(CC) -c $(CFLAGS) Structure.cpp

.PHONY: clean cleanest

clean:
	rm *.o

cleanest: clean
	rm whynot
