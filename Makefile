CC=g++
CFLAGS=-O2 -std=c++11 -march=native -Wall

EIGENDIR=/usr/include/eigen3
LIBINTDIR=/opt/libint

INCLUDE=-isystem$(EIGENDIR) -isystem$(LIBINTDIR)/include -I./include
LDFLAGS=-L$(LIBINTDIR)/lib -Wl,-Bstatic -lint2 -Wl,-Bdynamic -lpthread 

OBJ=obj/main.o obj/properties.o obj/qm_residue.o obj/grid.o obj/algebra.o

all: tresp

tresp: $(OBJ)
	$(CC) $^  $(LDFLAGS) -o $@

tresp_bu: $(OBJ)
	$(CC) $^  $(LDFLAGS) -o $@

obj/%.o: src/%.cc 
	$(CC) -c $(CFLAGS) $(INCLUDE) $< -o $@

clean:
	rm -rf $(OBJ)
	rm -rf tresp

