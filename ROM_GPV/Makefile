
#CC = gcc
CC = clang
CFLAGS = -Wall -O3 -march=native -mtune=native -fomit-frame-pointer -ftree-vectorize -funsafe-math-optimizations -mfpmath=sse
CVECFLAGS := $(CFLAGS) -mavx2 -ftree-vectorize #-fopt-info-vec-optimized



EXEC = signature_tests sampling_tests timing main_signature
OBJ = arithmetic.o random.o
HDR = common.h


all: signature_tests timing

timing: timing.o hash.o signature.o sampling.o random.o arithmetic.o cpucycles.o
	$(CC) $(CFLAGS) -o $@ $^ -lm

signature_tests: signature_tests.o hash.o signature.o sampling.o random.o arithmetic.o
	$(CC) $(CFLAGS) -o $@ $^ -lm

main_signature: main_signature.o signature.o sampling.o random.o $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^ -lm

sampling_tests: sampling_tests.o sampling.o random.o arithmetic.o
	$(CC) $(CFLAGS) -o $@ $^ -lm

arithmetic.o: arithmetic.c random.o $(HDR)
	$(CC) $(CVECFLAGS) -c -o $@ $<

sampling.o: sampling.c random.o $(HDR)
	$(CC) $(CVECFLAGS) -c -o $@ $<

random.o: random.c $(HDR)
	$(CC) $(CVECFLAGS) -maes -c -o $@ $<

timing.o : timing.c common.h
	$(CC) $(CFLAGS) -c -o $@ $<

%.o: %.c $(HDR)
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm -f $(EXEC) *.o *.s
