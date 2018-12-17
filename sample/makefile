TARGET= ansmaker_test minlex_test solve_test
all: $(TARGET)

CC=g++
CPPFLAGS= -O3 -march=native

%: %.cpp ../grid.cpp
	$(CC) $(CPPFLAGS) $^ -o $@

test: $(TARGET)
	./ansmaker_test > ansmaker.test
	diff ansmaker.answer ansmaker.test
	time ./minlex_test minlex.input > minlex.test
	diff minlex.answer minlex.test
	time ./solve_test solve.input > solve.test
	diff solve.answer solve.test

.PHONY: clean

clean:
	rm -f $(TARGET)
	rm -f ansmaker.test minlex.test solve.test
