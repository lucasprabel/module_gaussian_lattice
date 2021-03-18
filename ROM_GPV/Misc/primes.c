#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include <gmp.h>

int main(int argc, char **argv)
	{
	mpz_t N;
	uint64_t n;
	
	mpz_init(N);
	
	printf("Finding sparse primes congruent to 17 mod 32...\n");
	
	printf("Finding primes of the form 2^i + 2^4 + 1...\n");
	for(uint64_t i = 5 ; i < 64 ; i++)
		{
		n = (((uint64_t) 1)<<i) + 17;
		mpz_set_ui(N, n);
		switch(mpz_probab_prime_p(N, 20))
			{
			case 1:
				printf("i = %" PRIu64 " : %" PRIu64 " is probably prime\n", i, n);
				break;
			
			case 2:
				printf("i = %" PRIu64 " : %" PRIu64 " is prime\n", i, n);
				break;
			}
		}
	
	printf("Finding primes of the form 2^i - 2^4 + 1...\n");
	for(uint64_t i = 5 ; i < 64 ; i++)
		{
		n = (((uint64_t) 1)<<i) - 15;
		mpz_set_ui(N, n);
		switch(mpz_probab_prime_p(N, 20))
			{
			case 1:
				printf("i = %" PRIu64 " : %" PRIu64 " is probably prime\n", i, n);
				break;
			
			case 2:
				printf("i = %" PRIu64 " : %" PRIu64 " is prime\n", i, n);
				break;
			}
		}
	
	return EXIT_SUCCESS;
	}
