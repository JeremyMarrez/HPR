CXX=g++ -std=c++11
CFLAGS=-Wall -Wextra -Ofast -march=native
HEADERS=Params.h HPR.h
OBJECTS=HPR.o test.o
LIBRARIES=-lntl -lgmp -lpthread

test: $(OBJECTS)
	$(CXX) -o $@ $^ $(CFLAGS) $(LIBRARIES)

%.o: %.cpp $(HEADERS)
	$(CXX) -c -o $@ $(CFLAGS) $<
