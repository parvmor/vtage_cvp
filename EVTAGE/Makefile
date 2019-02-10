CC = g++
OPT = -O3 -g
LIBS = -lz -lcvp -lz
FLAGS = -std=c++11 -L../cvp_kit $(LIBS) $(OPT) -no-pie

OBJ = mypredictor.o
DEPS = ../cvp_kit/cvp.h mypredictor.h

all: cvp

cvp: $(OBJ)
	$(CC) $(FLAGS) -o $@ $^

%.o: %.cc $(DEPS)
	$(CC) $(FLAGS) -c -o $@ $<

.PHONY: clean

clean:
	rm -f *.o cvp
