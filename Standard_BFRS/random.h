#include <inttypes.h>
#include <x86intrin.h>

#include "common.h"

uint32_t uniform_int_distribution(uint32_t n);

scalar uniform_mod_q(void);

/*
	Code from random_aesni.c
*/

//public API
void random_bytes_init(void);

void random_bytes(uint8_t * restrict data);

/*
	Code from exp_aes.cpp
*/

double algorithm_EA(uint64_t * n);

/*
	Code from algoF_aes.cpp
*/

int algorithmF(double mu, double sigma);
