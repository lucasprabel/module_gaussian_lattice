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
		
		crt_representation(I_dk_ii, LOG_R);
		}
	
	// build I_d in the CRT domain
	for(int i = 0 ; i < PARAM_D ; ++i)
		{
		poly_matrix IA_i = poly_matrix_element(IA, PARAM_M, i, 0);
		poly_matrix IA_ii = poly_matrix_element(IA, PARAM_M, i, i);
		
		memset(IA_i, 0, PARAM_N * PARAM_D * sizeof(scalar));
		IA_ii[0] = 1;
		
		crt_representation(IA_ii, LOG_R);
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
		
		mul_crt_poly_matrix(AtimesTI, IA, TI, PARAM_D, PARAM_M, PARAM_D * PARAM_K, LOG_R);
		
		if(!is_zero_poly(AtimesTI, PARAM_N * PARAM_D * PARAM_D * PARAM_K - 1))
			{
			printf("AtimesTI\n");
			print_poly_matrix(AtimesTI, PARAM_D * PARAM_K, PARAM_D);
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

void test_construct_and_deconstruct_A_m(int n)
	{
	scalar *A_coeffs = malloc(PARAM_N * PARAM_D * (PARAM_M - PARAM_D) * sizeof(scalar));
	scalar *A_copy = malloc(PARAM_N * PARAM_D * (PARAM_M - PARAM_D) * sizeof(scalar));
	scalar m_coeffs[SMALL_DEGREE];
	
	poly A = A_coeffs, m = m_coeffs;
	
	// Generate some "public key"-like A and message-like m
	printf("Generating A and m...\n");
	
	random_poly(A, PARAM_N * PARAM_D * (PARAM_M - PARAM_D) - 1);
	random_poly(m, SMALL_DEGREE - 1);
	memcpy(A_copy, A, PARAM_N * PARAM_D * (PARAM_M - PARAM_D) * sizeof(scalar));
	
	printf("Done\n");
	
	// Construct and deconstruct A_m, and check that it hasn't changed
	for(int i = 0 ; i < n ; ++i)
		{
		printf("%d\n", i);
		construct_A_m(A, m);
		deconstruct_A_m(A, m);
		
		for(int j = 0 ; j < PARAM_D ; ++j)
			{
			for(int k = 0 ; k < PARAM_M - PARAM_D ; ++k)
				{
				for(int l = 0 ; l < PARAM_N ; ++l)
					{
					poly A_jk = poly_matrix_element(A, PARAM_M - PARAM_D, j, k);
					poly A_copy_jk = poly_matrix_element(A_copy, PARAM_M - PARAM_D, j, k);
					
					if(A_copy_jk[l] != A_jk[l])
						{
						printf("(%d, %d, %d) : %" PRIu32 " != %" PRIu32 "\n", j, k, l, A_copy_jk[l], A_jk[l]);
						}
					}
				}
			}
		}
	
	free(A);
	free(A_copy);
	}

void test_Sign_and_Verify(int n)
	{
	scalar *A_coeffs = malloc(PARAM_N * PARAM_D * (PARAM_M - PARAM_D) * sizeof(scalar)), *T_coeffs = malloc(PARAM_N * 2 * PARAM_D * PARAM_D * PARAM_K * sizeof(scalar));
	cplx *sch_comp_coeffs = malloc(PARAM_N * PARAM_D * (2 * PARAM_D + 1) * sizeof(cplx)), *cplx_T_coeffs = malloc(PARAM_N * 2 * PARAM_D * PARAM_D * PARAM_K * sizeof(cplx));
	poly_matrix A = A_coeffs, T = T_coeffs;
	cplx_poly_matrix sch_comp = sch_comp_coeffs, cplx_T = cplx_T_coeffs;
	
	
	scalar *A_copy = malloc(PARAM_N * PARAM_D * (PARAM_M - PARAM_D) * sizeof(scalar)), *T_copy = malloc(PARAM_N * 2 * PARAM_D * PARAM_D * PARAM_K * sizeof(scalar));
	cplx *sch_comp_copy = malloc(PARAM_N * PARAM_D * (2 * PARAM_D + 1) * sizeof(cplx)), *cplx_T_copy = malloc(PARAM_N * 2 * PARAM_D * PARAM_D * PARAM_K * sizeof(cplx));
	
	// Generate the keys
	printf("Generating keys (A, T, cplx_T, sch_comp)...\n");
	
	KeyGen(A, T, cplx_T, sch_comp);
	
	printf("Done\n");
	
	
	// Copy them to later check if they have changed
	memcpy(A_copy, A, PARAM_N * PARAM_D * (PARAM_M - PARAM_D) * sizeof(scalar));
	memcpy(T_copy, T, PARAM_N * 2 * PARAM_D * PARAM_D * PARAM_K * sizeof(scalar));
	memcpy(cplx_T_copy, cplx_T, PARAM_N * 2 * PARAM_D * PARAM_D * PARAM_K * sizeof(cplx));
	memcpy(sch_comp_copy, sch_comp, PARAM_N * PARAM_D * (2 * PARAM_D + 1) * sizeof(cplx));

	
	// Generate a message (a polynomial of degree < n/r)
	printf("Generating a message m...\n");
	
	scalar m_coeffs[PARAM_N] = {0};
	poly m = m_coeffs;
	
	random_poly(m, SMALL_DEGREE - 1);
	
	printf("Done\n");
	
	// Sample a signature nu
	scalar nu_coeffs[PARAM_N * PARAM_M];
	poly_matrix nu = nu_coeffs;
	
	printf("Signing...\n");
	
	for(int i = 0 ; i < n ; ++i)
		{
		Sign(nu, A, T, cplx_T, sch_comp, m);
		
		for(int j = 0 ; j < PARAM_D ; ++j)
			{
			for(int k = 0 ; k < PARAM_M - PARAM_D ; ++k)
				{
				for(int l = 0 ; l < PARAM_N ; ++l)
					{
					poly A_jk = poly_matrix_element(A, PARAM_M - PARAM_D, j, k);
					poly A_copy_jk = poly_matrix_element(A_copy, PARAM_M - PARAM_D, j, k);
					
					if(A_copy_jk[l] != A_jk[l])
						{
						printf("(%d, %d, %d) : %" PRIu32 " != %" PRIu32 "\n", j, k, l, A_copy_jk[l], A_jk[l]);
						}
					}
				}
			}

		if(Verify(nu, A, m) == false)
			{
			printf("verification failed @ %d\n", i);
			}
		else
			{
			printf("verification successful @ %d\n", i);
			}
		
		if(memcmp(A, A_copy, PARAM_N * PARAM_D * (PARAM_M - PARAM_D) * sizeof(scalar)) != 0)
			{
			printf("A changed after verifying\n");
			}

		}
	
	free(A);
	free(T);
	free(sch_comp);
	free(cplx_T);
	
	free(A_copy);
	free(T_copy);
	free(sch_comp_copy);
	free(cplx_T_copy);
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
	
	// Generate a message (a polynomial of degree < n/r)
	printf("Generating a message m... ");
	fflush(stdout);
	
	scalar m_coeffs[PARAM_N] = {0};
	poly m = m_coeffs;
	
	random_poly(m, SMALL_DEGREE - 1);
	
	printf("Done.\n");
	
	// Sample a signature nu
	scalar nu_coeffs[PARAM_N * PARAM_M];
	poly_matrix nu = nu_coeffs;
	
	printf("Signing the same message multiple times... ");
	fflush(stdout);
	
	t0 = clock();
	for(int i = 0 ; i < n ; ++i)
		{
		Sign(nu, A, T, cplx_T, sch_comp, m);
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
	
	// Generate a message (a polynomial of degree < n/r)
	printf("Generating a message m... ");
	fflush(stdout);
	
	scalar m_coeffs[PARAM_N] = {0};
	poly m = m_coeffs;
	
	random_poly(m, SMALL_DEGREE - 1);
	
	printf("Done.\n");
	
	// Sample a signature nu
	scalar nu_coeffs[PARAM_N * PARAM_M];
	poly_matrix nu = nu_coeffs;
	
	printf("Signing different messages... ");
	fflush(stdout);
	
	t1 = 0;
	for(int i = 0 ; i < n ; ++i)
		{
		// Generate a message (a polynomial of degree < n/r)
		random_poly(m, SMALL_DEGREE - 1);
		
		t0 = clock();
		
		// Sign it
		Sign(nu, A, T, cplx_T, sch_comp, m);
		
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
	
	// Generate a message (a polynomial of degree < n/r)
	printf("Generating a message m... ");
	fflush(stdout);
	
	scalar m_coeffs[PARAM_N] = {0};
	poly m = m_coeffs;
	
	random_poly(m, SMALL_DEGREE - 1);
	
	printf("Done.\n");
	
	// Sample a signature nu
	scalar nu_coeffs[PARAM_N * PARAM_M];
	poly_matrix nu = nu_coeffs;
	
	printf("Signing the message... ");
	fflush(stdout);
	
	Sign(nu, A, T, cplx_T, sch_comp, m);
	
	printf("Done.\n");
	
	// Verify the signature
	printf("Verifying the same signature multiple times... ");
	fflush(stdout);
	
	t0 = clock();
	for(int i = 0 ; i < n ; ++i)
		{
		Verify(nu, A, m);
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
	
	printf("Verifying different signatures for different messages... ");
	fflush(stdout);
	
	t1 = 0;
	for(int i = 0 ; i < n ; ++i)
		{
		// Generate a message
		scalar m_coeffs[PARAM_N] = {0};
		poly m = m_coeffs;
		
		random_poly(m, SMALL_DEGREE - 1);
		
		// Sign it
		scalar nu_coeffs[PARAM_N * PARAM_M];
		poly_matrix nu = nu_coeffs;
		
		Sign(nu, A, T, cplx_T, sch_comp, m);
		
		// Verify the signature
		t0 = clock();
		Verify(nu, A, m);
		t1 += clock() - t0;
		}
	
	printf("Done : %ld\n\n", t1);
	
	free(A);
	free(T);
	free(sch_comp);
	free(cplx_T);
	
	return t1;
	}

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
	printf("\tr = %d\n", PARAM_R);
	printf("\tq = %d\n", PARAM_Q);
	printf("\tk = %d\n", PARAM_K);
	printf("\td = %d\n", PARAM_D);
	printf("\tsigma = %f\n", PARAM_SIGMA);
	printf("\tzeta = %f\n", PARAM_ZETA);
	printf("\talpha = %f\n", PARAM_ALPHA);
	printf("==================================\n\n");
	
	init_crt_trees();
	init_D_lattice_coeffs();
	init_cplx_roots_of_unity();
	
	//random_bytes_init();
	
	//test_KeyGen(n);
	test_Sign_and_Verify(n);
	//test_construct_and_deconstruct_A_m(n);
	
	//time_KeyGen1(n);
	//time_Sign1(n);
	//time_Sign2(n);
	//time_Verify1(n);
	//time_Verify2(n);
	
	return 0;
	}
