OPTS = -std=c++11 -O3 -fstack-protector-all -g -W -Wall -Wextra -Wunused -Wcast-align -Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security -Wmissing-format-attribute -Wformat=2 -Wwrite-strings -Wcast-align -Wno-long-long -Woverloaded-virtual -Wnon-virtual-dtor -Wcast-qual -Wno-property-attribute-mismatch -pthread
all: clear a.out
a.out: main.o functions.o matrix.o solver.o
	g++ $(OPTS) $^ -o a.out
main.o: main.cpp functions.h matrix.h solver.h approximation.h common.h
	g++ -c $(OPTS) $<
functions.o: functions.cpp functions.h
	g++ -c $(OPTS) $<
solver.o: solver.cpp solver.h
	g++ -c $(OPTS) $<
matrix.o: matrix.cpp matrix.h
	g++ -c $(OPTS) $<
clear:
	rm -f *.o
clean:
	rm -f *.out *.o *.bak
