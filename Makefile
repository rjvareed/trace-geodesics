CC=g++
LDFLAGS=-lgsl -lginac -lcln
CFLAGS=-g -Wall `pkg-config --cflags --libs gtkmm-4.0` -std=c++17

all: main.cpp
	$(CC) -c DrawingArea.cpp $(CFLAGS)
	$(CC) -c main.cpp $(CFLAGS) $(LDFLAGS)
	$(CC) main.o DrawingArea.o -o app $(CFLAGS) $(LDFLAGS)
	rm -f *.o

clean:
	rm -f *.o
	rm app
