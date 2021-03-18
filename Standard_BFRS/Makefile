
#CC = gcc
CC = clang
CFLAGS = -Wall -O3 -mtune=native
CVECFLAGS := $(CFLAGS) -mavx2 -ftree-vectorize #-fopt-info-vec-optimized
CFLAGS += -fno-tree-vectorize
NFL_AVX = -DNFL_OPTIMIZED=ON -DNTT_AVX2

EXEC = signature_tests sampling_tests timing main_signature

OBJ = arithmetic.o random.o crt_trees.o

HDR = common.h

all: sampling_tests signature_tests

timing: timing.o crt_trees.o signature.o sampling.o random.o arithmetic.o cpucycles.o
	$(CC) $(CFLAGS) -o $@ $^ -lm

signature_tests: signature_tests.o signature.o sampling.o random.o arithmetic.o crt_trees.o
	$(CC) $(CFLAGS) -o $@ $^ -lm

main_signature: main_signature.o signature.o sampling.o random.o $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^ -lm

sampling_tests: sampling_tests.o sampling.o random.o arithmetic.o crt_trees.o
	$(CC) $(CFLAGS) -o $@ $^ -lm


arithmetic.o: arithmetic.c random.o $(HDR)
	$(CC) $(CVECFLAGS) -c -o $@ $<

sampling.o: sampling.c random.o $(HDR)
	$(CC) $(CVECFLAGS) -c -o $@ $<

random.o: random.c
	$(CC) $(CVECFLAGS) -maes -c -o $@ $<

%.o: %.c $(HDR)
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm -f $(EXEC) *.o *.s
