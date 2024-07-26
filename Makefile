CC=gcc
LDFLAGS=-lgsl

all: main.c
	$(CC) main.c -o app $(LDFLAGS)

clean:
	rm app
