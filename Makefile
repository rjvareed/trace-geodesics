CC=g++
LDFLAGS=-lgsl -lginac -lcln

all: main.cpp
	$(CC) main.cpp -o app $(LDFLAGS)

clean:
	rm app
