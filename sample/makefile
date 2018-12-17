all: ansmaker_test minlex_test

CC=g++
CPPFLAGS= -O3 -march=native

ansmaker_test: ansmaker_test.cpp ../grid.cpp
	$(CC) $(CPPFLAGS) $^ -o $@

minlex_test: minlex_test.cpp ../grid.cpp
	$(CC) $(CPPFLAGS) $^ -o $@

test: ansmaker_test minlex_test
	./ansmaker_test > ansmaker_test.test
	diff ansmaker_test.orig ansmaker_test.test
	time ./minlex_test sample.txt > sample.test
	diff sample.minlex sample.test

.PHONY: clean

clean:
	rm -f ansmaker_test minlex_test
	rm -f ansmaker_test.test minlex_test.test
