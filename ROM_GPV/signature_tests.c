#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include <time.h>
#include <string.h>

#include "common.h"
#include "random.h"
#include "sampling.h"
#include "signature.h"
#include "arithmetic.h"
#include "hash.h"

#define MESSAGE_BYTES 512 // keep it divisible by 16

void test_KeyGen(int n)
	{
	scalar *A_coeffs = malloc(PARAM_N * PARAM_D * (PARAM_M - PARAM_D) * sizeof(scalar));
	scalar *IA_coeffs = malloc(PARAM_N * PARAM_D * PARAM_M * sizeof(scalar));
	scalar *TI_coeffs = malloc(PARAM_N * PARAM_M * PARAM_D * PARAM_K * sizeof(scalar));
	scalar *AtimesTI_coeffs = malloc(PARAM_N * PARAM_D * PARAM_D * PARAM_K * sizeof(scalar));
	cplx *sch_comp_coeffs = malloc(PARAM_N * PARAM_D * (2 * PARAM_D + 1) * sizeof(cplx));
	cplx *cplx_T_coeffs = malloc(PARAM_N * 2 * PARAM_D * PARAM_D * PARAM_K * sizeof(cplx));
	
	poly_matrix A = A_coeffs, IA = IA_coeffs, TI = TI_coeffs, T = TI_coeffs, I_dk = poly_matrix_element(TI, PARAM_D * PARAM_K, 2 * PARAM_D, 0), AtimesTI = AtimesTI_coeffs;
	cplx_poly_matrix sch_comp = sch_comp_coeffs, cplx_T = cplx_T_coeffs;
	
	// build I_dk in the CRT domain
	memset(I_dk, 0, PARAM_N * PARAM_D * PARAM_K * PARAM_D * PARAM_K * sizeof(scalar));
	for(int i = 0 ; i < PARAM_D * PARAM_K ; ++i)
		{
		poly I_dk_ii = poly_matrix_element(I_dk, PARAM_D * PARAM_K, i, i);
		
		I_dk_ii[0] = 1;
		
		ntt(I_dk_ii);
		}
	
	// build I_d in the CRT domain
	for(int i = 0 ; i < PARAM_D ; ++i)
		{
		poly_matrix IA_i = poly_matrix_element(IA, PARAM_M, i, 0);
		poly_matrix IA_ii = poly_matrix_element(IA, PARAM_M, i, i);
		
		memset(IA_i, 0, PARAM_N * PARAM_D * sizeof(scalar));
		IA_ii[0] = 1;
		
		ntt(IA_ii);
		}
	
	printf("Generating and verifying keys\n");
	
	for(int i = 0 ; i < n ; ++i)
		{
		// Generate the keys
		KeyGen(A, T, cplx_T, sch_comp);
		
		// Test the keys
		for(int j = 0 ; j < PARAM_D ; ++j)
			{
			poly_matrix IA_jd = poly_matrix_element(IA, PARAM_M, j, PARAM_D);
			poly_matrix A_j = poly_matrix_element(A, PARAM_M - PARAM_D, j, 0);
			
			memcpy(IA_jd, A_j, PARAM_N * (PARAM_M - PARAM_D) * sizeof(scalar));
			}
		
		mul_crt_poly_matrix(AtimesTI, IA, TI, PARAM_D, PARAM_M, PARAM_D * PARAM_K);
		
		if(!is_zero_poly(AtimesTI, PARAM_N * PARAM_D * PARAM_D * PARAM_K - 1))
			{
			printf("AtimesTI\n");
			print_poly_matrix(AtimesTI, PARAM_D, PARAM_D * PARAM_K);
			}
		else
			{
			printf("ok @ %d\n", i);
			}
		}
	
	free(A);
	free(IA);
	free(TI);
	free(AtimesTI);
	free(sch_comp);
	free(cplx_T);
	}

int time_KeyGen1(int n)
	{
	clock_t t0, t1;
	
	scalar *A_coeffs = malloc(PARAM_N * PARAM_D * (PARAM_M - PARAM_D) * sizeof(scalar)), *T_coeffs = malloc(PARAM_N * 2 * PARAM_D * PARAM_D * PARAM_K * sizeof(scalar));
	cplx *sch_comp_coeffs = malloc(PARAM_N * PARAM_D * (2 * PARAM_D + 1) * sizeof(cplx)), *cplx_T_coeffs = malloc(PARAM_N * 2 * PARAM_D * PARAM_D * PARAM_K * sizeof(cplx));
	poly_matrix A = A_coeffs, T = T_coeffs;
	cplx_poly_matrix sch_comp = sch_comp_coeffs, cplx_T = cplx_T_coeffs;
	
	// Generate the keys
	printf("Generating keys (A, T, cplx_T, sch_comp)... ");
	fflush(stdout);
	
	t0 = clock();
	for(int i = 0 ; i < n ; ++i)
		{
		KeyGen(A, T, cplx_T, sch_comp);
		}
	t1 = clock() - t0;
	
	printf("Done : %ld\n\n", t1);
	
	free(A);
	free(T);
	free(sch_comp);
	free(cplx_T);
	
	return t1;
	}

int time_KeyGen2(int n)
	{
	clock_t t0, t1;
	
	scalar A_coeffs[PARAM_N * PARAM_D * (PARAM_M - PARAM_D)], T_coeffs[PARAM_N * 2 * PARAM_D * PARAM_D * PARAM_K];
	poly_matrix A = A_coeffs, T = T_coeffs;
	
	printf("Generating keys (A,T)...\n");
	
	t0 = clock();
	for(int i = 0 ; i < n ; ++i)
		{
		TrapGen(A, T);
		}
	t1 = clock() - t0;
	
	printf("Done : %ld\n\n", t1);
	
	free(A);
	free(T);
	
	return t1;
	}

void test_Sign_and_Verify(int n)
	{
	scalar *A_coeffs = malloc(PARAM_N * PARAM_D * (PARAM_M - PARAM_D) * sizeof(scalar)), *T_coeffs = malloc(PARAM_N * 2 * PARAM_D * PARAM_D * PARAM_K * sizeof(scalar));
	cplx *sch_comp_coeffs = malloc(PARAM_N * PARAM_D * (2 * PARAM_D + 1) * sizeof(cplx)), *cplx_T_coeffs = malloc(PARAM_N * 2 * PARAM_D * PARAM_D * PARAM_K * sizeof(cplx));
	poly_matrix A = A_coeffs, T = T_coeffs;
	cplx_poly_matrix sch_comp = sch_comp_coeffs, cplx_T = cplx_T_coeffs;
	
	// Generate the keys
	printf("Generating keys (A, T, cplx_T, sch_comp)...\n");
	
	KeyGen(A, T, cplx_T, sch_comp);
	
	printf("Done\n");
	
	
	// Generate a message (an array of MESSAGE_BYTES bytes)
	printf("Generating a message m...\n");
	
	uint8_t m[MESSAGE_BYTES];
	
	for(int i = 0 ; i < MESSAGE_BYTES / 16 ; ++i)
		{
		random_bytes(&m[16 * i]);
		}

	printf("Done\n");
	
	// Compute a signature (nu, r)
	uint8_t r[SALT_BYTES];
	scalar nu_coeffs[PARAM_N * PARAM_M];
	poly_matrix nu = nu_coeffs;
	
	printf("Signing...\n");
	
	for(int i = 0 ; i < n ; ++i)
		{
		Sign(nu, r, A, T, cplx_T, sch_comp, m, MESSAGE_BYTES);
		
		if(Verify(nu, r, A, m, MESSAGE_BYTES) == false)
			{
			printf("verification failed @ %d\n", i);
			}
		else
			{
			printf("verification successful @ %d\n", i);
			}
		}
	
	free(A);
	free(T);
	free(sch_comp);
	free(cplx_T);
	}

void test_H(int n)
	{
	int message_bytes = 512; // keep it divisible by 16
	uint8_t m[message_bytes];
	scalar t_coeffs[PARAM_N * PARAM_D];
	poly_matrix t = t_coeffs;
	
	printf("Generating a random message...");
	fflush(stdout);
	
	// Generate a random message
	for(int i = 0 ; i < message_bytes / 16 ; ++i)
		{
		random_bytes(&m[16 * i]);
		}
	
	printf("Done.\n");
	
	printf("m\n[");
	for(int i = 0 ; i < message_bytes ; ++i)
		{
		printf("%02" PRIX8 ", ", m[i]);
		}
	printf("]\n\n");
	
	printf("\"Hashing\" it multiple times...\n\n");
	
	uint8_t r[SALT_BYTES];
	
	for(int i = 0 ; i < n ; ++i)
		{
		printf("i = %d\n", i);
		
		salt(r);
		H(t, m, message_bytes, r);
		
		printf("r\n[");
		for(int j = 0 ; j < SALT_BYTES ; ++j)
			{
			printf("%02" PRIX8 ", ", r[j]);
			}
		printf("]\n\n");
		printf("t\n");
		print_poly_matrix(t, PARAM_D, 1);
		printf("\n\n");
		}
	}

int time_Sign1(int n)
	{
	clock_t t0, t1;
	
	scalar *A_coeffs = malloc(PARAM_N * PARAM_D * (PARAM_M - PARAM_D) * sizeof(scalar)), *T_coeffs = malloc(PARAM_N * 2 * PARAM_D * PARAM_D * PARAM_K * sizeof(scalar));
	cplx *sch_comp_coeffs = malloc(PARAM_N * PARAM_D * (2 * PARAM_D + 1) * sizeof(cplx)), *cplx_T_coeffs = malloc(PARAM_N * 2 * PARAM_D * PARAM_D * PARAM_K * sizeof(cplx));
	poly_matrix A = A_coeffs, T = T_coeffs;
	cplx_poly_matrix sch_comp = sch_comp_coeffs, cplx_T = cplx_T_coeffs;
	
	// Generate the keys
	printf("Generating keys (A, T, cplx_T, sch_comp)... ");
	fflush(stdout);
	
	KeyGen(A, T, cplx_T, sch_comp);
	
	printf("Done.\n");
	
	// Generate a message (an array of MESSAGE_BYTES bytes)
	printf("Generating a message m... ");
	fflush(stdout);
	
	uint8_t m[MESSAGE_BYTES];
	
	for(int i = 0 ; i < MESSAGE_BYTES / 16 ; ++i)
		{
		random_bytes(&m[16 * i]);
		}
	
	printf("Done.\n");
	
	// Compute a signature (nu, r)
	uint8_t r[SALT_BYTES];
	scalar nu_coeffs[PARAM_N * PARAM_M];
	poly_matrix nu = nu_coeffs;
	
	printf("Signing the same message multiple times... ");
	fflush(stdout);
	
	t0 = clock();
	for(int i = 0 ; i < n ; ++i)
		{
		Sign(nu, r, A, T, cplx_T, sch_comp, m, MESSAGE_BYTES);
		}
	t1 = clock() - t0;
	
	printf("Done : %ld\n\n", t1);
	
	free(A);
	free(T);
	free(sch_comp);
	free(cplx_T);
	
	return t1;
	}

int time_Sign2(int n)
	{
	clock_t t0, t1;
	
	scalar *A_coeffs = malloc(PARAM_N * PARAM_D * (PARAM_M - PARAM_D) * sizeof(scalar)), *T_coeffs = malloc(PARAM_N * 2 * PARAM_D * PARAM_D * PARAM_K * sizeof(scalar));
	cplx *sch_comp_coeffs = malloc(PARAM_N * PARAM_D * (2 * PARAM_D + 1) * sizeof(cplx)), *cplx_T_coeffs = malloc(PARAM_N * 2 * PARAM_D * PARAM_D * PARAM_K * sizeof(cplx));
	poly_matrix A = A_coeffs, T = T_coeffs;
	cplx_poly_matrix sch_comp = sch_comp_coeffs, cplx_T = cplx_T_coeffs;
	
	// Generate the keys
	printf("Generating keys (A, T, cplx_T, sch_comp)... ");
	fflush(stdout);
	
	KeyGen(A, T, cplx_T, sch_comp);
	
	printf("Done.\n");
	
	uint8_t m[MESSAGE_BYTES];
	uint8_t r[SALT_BYTES];
	scalar nu_coeffs[PARAM_N * PARAM_M];
	poly_matrix nu = nu_coeffs;
	
	printf("Signing different messages... ");
	fflush(stdout);
	
	t1 = 0;
	for(int i = 0 ; i < n ; ++i)
		{
		// Generate a message (an array of MESSAGE_BYTES bytes)
		for(int i = 0 ; i < MESSAGE_BYTES / 16 ; ++i)
			{
			random_bytes(&m[16 * i]);
			}
		
		t0 = clock();
		
		// Sign it
		Sign(nu, r, A, T, cplx_T, sch_comp, m, MESSAGE_BYTES);
		
		t1 += clock() - t0;
		}
	
	printf("Done : %ld\n\n", t1);
	
	free(A);
	free(T);
	free(sch_comp);
	free(cplx_T);
	
	return t1;
	}

int time_Verify1(int n)
	{
	clock_t t0, t1;
	
	scalar *A_coeffs = malloc(PARAM_N * PARAM_D * (PARAM_M - PARAM_D) * sizeof(scalar)), *T_coeffs = malloc(PARAM_N * 2 * PARAM_D * PARAM_D * PARAM_K * sizeof(scalar));
	cplx *sch_comp_coeffs = malloc(PARAM_N * PARAM_D * (2 * PARAM_D + 1) * sizeof(cplx)), *cplx_T_coeffs = malloc(PARAM_N * 2 * PARAM_D * PARAM_D * PARAM_K * sizeof(cplx));
	poly_matrix A = A_coeffs, T = T_coeffs;
	cplx_poly_matrix sch_comp = sch_comp_coeffs, cplx_T = cplx_T_coeffs;
	
	// Generate the keys
	printf("Generating keys (A, T, cplx_T, sch_comp)... ");
	fflush(stdout);
	
	KeyGen(A, T, cplx_T, sch_comp);
	
	printf("Done.\n");
	
	// Generate a message (an array of MESSAGE_BYTES bytes)
	printf("Generating a message m... ");
	fflush(stdout);
	
	uint8_t m[MESSAGE_BYTES];
	
	for(int i = 0 ; i < MESSAGE_BYTES / 16 ; ++i)
		{
		random_bytes(&m[16 * i]);
		}
	
	printf("Done.\n");
	
	// Compute a signature (nu, r)
	uint8_t r[SALT_BYTES];
	scalar nu_coeffs[PARAM_N * PARAM_M];
	poly_matrix nu = nu_coeffs;
	
	printf("Signing the message... ");
	fflush(stdout);
	
	Sign(nu, r, A, T, cplx_T, sch_comp, m, MESSAGE_BYTES);
	
	printf("Done.\n");
	
	// Verify the signature
	printf("Verifying the same signature multiple times... ");
	fflush(stdout);
	
	t0 = clock();
	for(int i = 0 ; i < n ; ++i)
		{
		Verify(nu, r, A, m, MESSAGE_BYTES);
		}
	t1 = clock() - t0;
	
	printf("Done : %ld\n\n", t1);
	
	free(A);
	free(T);
	free(sch_comp);
	free(cplx_T);
	
	return t1;
	}

int time_Verify2(int n)
	{
	clock_t t0, t1;
	
	scalar *A_coeffs = malloc(PARAM_N * PARAM_D * (PARAM_M - PARAM_D) * sizeof(scalar)), *T_coeffs = malloc(PARAM_N * 2 * PARAM_D * PARAM_D * PARAM_K * sizeof(scalar));
	cplx *sch_comp_coeffs = malloc(PARAM_N * PARAM_D * (2 * PARAM_D + 1) * sizeof(cplx)), *cplx_T_coeffs = malloc(PARAM_N * 2 * PARAM_D * PARAM_D * PARAM_K * sizeof(cplx));
	poly_matrix A = A_coeffs, T = T_coeffs;
	cplx_poly_matrix sch_comp = sch_comp_coeffs, cplx_T = cplx_T_coeffs;
	
	// Generate the keys
	printf("Generating keys (A, T, cplx_T, sch_comp)... ");
	fflush(stdout);
	
	KeyGen(A, T, cplx_T, sch_comp);
	
	printf("Done.\n");
	
	uint8_t m[MESSAGE_BYTES];
	uint8_t r[SALT_BYTES];
	scalar nu_coeffs[PARAM_N * PARAM_M];
	poly_matrix nu = nu_coeffs;
	
	printf("Verifying different signatures for different messages... ");
	fflush(stdout);
	
	t1 = 0;
	for(int i = 0 ; i < n ; ++i)
		{
		// Generate a message (an array of MESSAGE_BYTES bytes)
		for(int i = 0 ; i < MESSAGE_BYTES / 16 ; ++i)
			{
			random_bytes(&m[16 * i]);
			}
		
		// Sign it
		Sign(nu, r, A, T, cplx_T, sch_comp, m, MESSAGE_BYTES);
		
		// Verify the signature
		t0 = clock();
		
		Verify(nu, r, A, m, MESSAGE_BYTES);
		
		t1 += clock() - t0;
		}
	
	printf("Done : %ld\n\n", t1);
	
	free(A);
	free(T);
	free(sch_comp);
	free(cplx_T);
	
	return t1;
	}

// Time sample_pre and its subroutines (sample_perturb, sample_G, etc.)
void time_sample_pre_subroutines(int n)
	{
	clock_t t0, t1;
	
	// Generate the keys
	scalar *A_coeffs = malloc(PARAM_N * PARAM_D * (PARAM_M - PARAM_D) * sizeof(scalar)), *T_coeffs = malloc(PARAM_N * 2 * PARAM_D * PARAM_D * PARAM_K * sizeof(scalar));
	cplx *sch_comp_coeffs = malloc(PARAM_N * PARAM_D * (2 * PARAM_D + 1) * sizeof(cplx)), *cplx_T_coeffs = malloc(PARAM_N * 2 * PARAM_D * PARAM_D * PARAM_K * sizeof(cplx));
	poly_matrix A = A_coeffs, T = T_coeffs;
	cplx_poly_matrix sch_comp = sch_comp_coeffs, cplx_T = cplx_T_coeffs;
	
	printf("Generating keys (A, T, cplx_T, sch_comp)... ");
	fflush(stdout);
	
	KeyGen(A, T, cplx_T, sch_comp);
	
	printf("Done.\n\n");
	
	// Generate a dummy target for the sampling algorithms
	scalar v_coeffs[PARAM_N * PARAM_M], t_coeffs[PARAM_N * PARAM_D];
	poly_matrix v = v_coeffs, t = t_coeffs;
	
	printf("Generating a dummy target t... ");
	fflush(stdout);
	
	random_poly(t, PARAM_N * PARAM_D - 1);
	
	printf("Done.\n\n");
	
	// Time the G-sampling
	printf("G-sampling... ");
	fflush(stdout);
	
	t0 = clock();
	for(int i = 0 ; i < n ; ++i)
		{
		module_sample_G((signed_poly_matrix) v, t);
		}
	t1 = clock() - t0;
	
	printf("Done : %ld\n\n", t1);
	
	// Time the perturbation sampling
	printf("Sampling perturbations... ");
	fflush(stdout);
	
	t0 = clock();
	for(int i = 0 ; i < n ; ++i)
		{
		sample_perturb((signed_poly_matrix) v, cplx_T, sch_comp);
		}
	t1 = clock() - t0;
	
	printf("Done : %ld\n\n", t1);
	
	// Time the complete preimage sampling
	printf("Sampling preimages... ");
	fflush(stdout);
	
	t0 = clock();
	for(int i = 0 ; i < n ; ++i)
		{
		sample_pre(v, A, T, cplx_T, sch_comp, t);
		}
	t1 = clock() - t0;
	
	printf("Done : %ld\n\n", t1);
	
	
	}

unsigned long long timing_overhead;
unsigned long long timing_sampleZ_G = 0;
unsigned long long timing_sampleZ_P = 0;
unsigned long long timing_sampleG = 0;
unsigned long long timing_samplePerturb = 0;
unsigned long long timing_sampleArith = 0;

unsigned long long timing_sign = 0;

unsigned long long timing_sampleZ_KG = 0;
unsigned long long timing_precomp_KG = 0;
unsigned long long timing_arith_KG = 0;
unsigned long long timing_keygen = 0;

int main(int argc, char **argv)
	{
	if(argc < 2)
		{
		fprintf(stderr, "Usage : %s n\n", argv[0]);
		exit(EXIT_FAILURE);
		}
	
	int n = atoi(argv[1]);
	
	printf("=========== Parameters ===========\n");
	printf("\tn = %d\n", PARAM_N);
	printf("\tq = %d\n", PARAM_Q);
	printf("\tk = %d\n", PARAM_K);
	printf("\td = %d\n", PARAM_D);
	printf("\tsigma = %f\n", PARAM_SIGMA);
	printf("\tzeta = %f\n", PARAM_ZETA);
	printf("\talpha = %f\n", PARAM_ALPHA);
	printf("==================================\n\n");
	
	init_D_lattice_coeffs();
	init_cplx_roots_of_unity();
	
	//random_bytes_init();
	
	//test_KeyGen(n);
	test_Sign_and_Verify(n);
	//test_H(n);
	
	//time_KeyGen1(n);
	//time_Sign1(n);
	//time_Sign2(n);
	//time_Verify1(n);
	//time_Verify2(n);
	
	//time_sample_pre_subroutines(n);
	
	return 0;
	}
