all: a.out
SRC=$(shell ls *.cpp ../*.cpp)

CC=g++
CPPFLAGS= -O3 -march=native

a.out: $(SRC)
	$(CC) $(CPPFLAGS) $^ -o $@

clean:
	rm -f a.out
