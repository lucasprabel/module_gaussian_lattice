CC = gcc
CFLAGS = -g -Wall

EXEC = multiplicative_order primes

all: $(EXEC)

multiplicative_order: multiplicative_order.c
	$(CC) $(CFLAGS) -o $@ $^

primes: primes.c
	$(CC) $(CFLAGS) -o $@ $^ -lgmp

clean:
	rm -f $(EXEC)
