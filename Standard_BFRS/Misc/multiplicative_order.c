#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>

// a^x mod m
uint64_t mod_exp(uint64_t a, uint64_t x, uint64_t m)
	{
	uint64_t r = 1;
	
	for(int i=0 ; i < 64 ; i++)
		{
		if(x & 1)
			{
			r = (((r * a) % m) + m) % m;
			}
		a = (((a * a) % m) + m) % m;
		x = x >> 1;
		}
	
	return r;
	}

// multiplicative order of a mod 2^e
uint64_t order(uint64_t a, uint64_t e)
	{
	uint64_t ones = (1 << e) - 1;
	
	for(int i=0 ; i < e ; i++)
		{
		if(a == 1)
			{
			return (1 << i);
			}
		a = (a * a) & ones;
		}
	
	return 0;
	}

int main(int argc, char **argv)
	{
	if(argc < 4)
		{
		fprintf(stderr, "Usage : %s a e d\n", argv[0]);
		exit(EXIT_FAILURE);
		}
	
	uint64_t a, e, a_i, d;
	
	sscanf(argv[1], "%" SCNu64, &a);
	sscanf(argv[2], "%" SCNu64, &e);
	sscanf(argv[3], "%" SCNu64, &d);
	
	uint64_t ones = (1 << e) - 1;
	
	for(int i=1 ; i <= d ; i++)
		{
		a_i = mod_exp(a, i, (1 << e));
		printf("%" PRIu64 " = %" PRIu64 "^%d has multiplicative order %" PRIu64 " mod 2^%" PRIu64 "\n", a_i, a, i, order(a_i, e), e);
		}
	
	return EXIT_SUCCESS;
	}
