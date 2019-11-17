TARGET = pdeSolver
LIBS = -lm -llikwid
CC = gcc
CFLAGS = -O3 -D LIKWID_PERFMON -I /home/soft/likwid/include -L /home/soft/likwid/lib -mavx -march=native 

.PHONY: default all clean

default: $(TARGET)
all: default

OBJECTS = $(patsubst %.c, %.o, $(wildcard *.c))
HEADERS = $(wildcard *.h)

%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

.PRECIOUS: $(TARGET) $(OBJECTS)

$(TARGET): $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) $(LIBS) -o $@

clean:
	-rm -f *.o core $(TARGET) *.gch *.dat *.txt
	-rm -rf html/ latex/

doc:
	doxygen Doxyfile
