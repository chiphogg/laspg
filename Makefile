CC=./gfilt
OBJS=g_of_r.o lattice.o
CFLAGS=-ggdb

all: g_of_r

g_of_r: $(OBJS)
	$(CC) $(LFLAGS) $(CFLAGS) $(OBJS) -o g_of_r

g_of_r.o : g_of_r.cc lattice.h
	$(CC) $(LFLAGS) $(CFLAGS) -c g_of_r.cc

lattice.o : lattice.cc lattice.h
	$(CC) $(LFLAGS) $(CFLAGS) -c lattice.cc

clean:
	rm -Rf *.o
