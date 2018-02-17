CC=g++
CFLAGS=-O2 -std=c++11 -march=native

EIGENDIR=/usr/include/eigen3
LIBINTDIR=/opt/libint

INCLUDE=-I$(EIGENDIR) -I$(LIBINTDIR)/include -I./include
LDFLAGS=-L$(LIBINTDIR)/lib -Wl,-Bstatic -lint2 -Wl,-Bdynamic -lpthread 

OBJ=obj/main.o obj/properties.o obj/qm_residue.o

all: tresp

tresp: $(OBJ)
	$(CC) $^  $(LDFLAGS) -o $@
	
obj/%.o: src/%.cc 
	$(CC) -c $(CFLAGS) $(INCLUDE) $< -o $@

clean:
	rm -rf $(OBJ)
	rm -rf tresp
