CC = g++
OPT_DBG = -O3 -g -Wall -Wextra -pedantic -O3 -Wshadow -Wformat=2 -Wfloat-equal -Wconversion -Wlogical-op -Wfatal-errors -Wshift-overflow=2 -Wduplicated-cond -Wcast-qual -Wcast-align -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC -D_FORTIFY_SOURCE=2 -fsanitize=address -fsanitize=undefined -fno-sanitize-recover -fstack-protector
OPT = -O3 -g
LIBS = -lz -lcvp -lz
FLAGS = -std=c++14 -L../cvp_kit $(LIBS) $(OPT) -no-pie

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
