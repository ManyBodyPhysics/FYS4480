CC= clang-omp++
CFLAGS= -Wall -Wextra -pedantic -std=c++11 -fopenmp -g -I/usr/local/opt/openblas/include
LIBS= -lopenblas -lm  -L/usr/local/opt/openblas/lib

EXECUTABLES = main.exe

O_FILES = $(patsubst %.cpp, %.o, $(wildcard *.cpp))

$(EXECUTABLES): $(O_FILES)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

%.o: %.cpp  
	$(CC) $(CFLAGS) -c -o $@ $< 

clean:
	rm *.o *~ *# $(EXECUTABLES) 
 
