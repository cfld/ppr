CFLAGS += -mtune=native -march=native -Wall -O3 -DNDEBUG -g -fopenmp
LIBS += -lrt

all: ppr

ppr: ppr.c
	gcc $(CFLAGS) -o bin/ppr ppr.c

clean:
	rm -rf bin/* *.o
