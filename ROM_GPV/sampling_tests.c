#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include <time.h>
#include <string.h>

#include "common.h"
#include "random.h"
#include "sampling.h"
#include "arithmetic.h"

void test_random_bytes(int n)
	{
	uint8_t r[16];
	
	for(int i = 0 ; i < n ; ++i)
		{
		random_bytes(r);
		for(int j = 0 ; j < 16 ; ++j)
			{
			printf("%d ", r[j]);
			}
		printf("\n");
		}
	}

void time_SampleZ(int n)
	{
	clock_t t0, t1;
	real c = 0.5, sigma = PARAM_SIGMA;
	
	printf("Sampling from the discrete Gaussian over Z with center %f and parameter %f...\n", c, sigma);
	
	t0 = clock();
	for(int i = 0 ; i < n ; ++i)
		{
		SampleZ(c, sigma);
		}
	t1 = clock() - t0;
	
	printf("Done : %ld\n", t1);
	}

void test_SampleZ(int n)
	{
	real c = -10, sigma = PARAM_SIGMA;
	
	printf("Sampling from the discrete Gaussian over Z with center %f and parameter %f...\n", c, sigma);
	
	for(int i = 0 ; i < n ; ++i)
		{
		printf("%d\n", SampleZ(c, sigma));
		}
	}

void test_sample_D(int n)
	{
	signed_scalar z[PARAM_K];
	real c[PARAM_K], sigma = PARAM_ALPHA / 3;
	
	for(int i = 0 ; i < PARAM_K ; ++i)
		{
		c[i] = 0;
		}
	
	for(int i = 0 ; i < n ; ++i)
		{
		sample_D(z, c, sigma);
		print_signed_poly(z, PARAM_K - 1);
		printf("\n");
		}
	}

void test_scalar_sample_G(int n)
	{
	signed_scalar z[PARAM_K];
	scalar u, prod;
	
	for(int i = 0 ; i < n ; ++i)
		{
		random_poly(&u, 0);
		scalar_sample_G(z, u);
		
		//printf("\n");
		print_signed_poly(z, PARAM_K - 1);
		
		multiply_by_scalar_gadget_vector(&prod, z);
		
		printf("(%" PRIu32 ", %" PRIu32 ")\n", u, prod);
		}
	}

void test_ring_sample_G(int n)
	{
	signed_scalar z_coeffs[PARAM_K * PARAM_N];
	scalar prod_coeffs[PARAM_N], u_coeffs[PARAM_N];
	signed_poly_matrix z = z_coeffs;
	poly prod = prod_coeffs, u = u_coeffs;
	
	for(int i = 0 ; i < n ; ++i)
		{
		random_poly(u, PARAM_N - 1);
		// Sample z such that <g,z> = 0 mod q
		ring_sample_G(z, u);
		
		
		// Check that <g,z> = u mod q
		multiply_by_ring_gadget_vector(prod, z);
		for(int j = 0 ; j < PARAM_N ; ++j)
			{
			if(prod[j] != u[j])
				{
				printf("%d %d : %" PRIu32 "\n", i, j, prod[j]);
				}
			}
		}
	}

void test_module_sample_G(int n)
	{
	signed_scalar z_coeffs[PARAM_D * PARAM_K * PARAM_N];
	scalar prod_coeffs[PARAM_D * PARAM_N], u_coeffs[PARAM_D * PARAM_N];
	signed_poly_matrix z = z_coeffs;
	poly_matrix prod = prod_coeffs, u = u_coeffs;
	
	for(int i = 0 ; i < n ; ++i)
		{
		random_poly(u, PARAM_D * PARAM_N - 1);
		// Sample z such that Gz = 0 mod q
		module_sample_G(z, u);
		
		// Check that Gz = u mod q
		multiply_by_module_gadget_matrix(prod, z);
		for(int j = 0 ; j < PARAM_D * PARAM_N ; ++j)
			{
			if(prod[j] != u[j])
				{
				printf("%d %d : %" PRIu32 "\n", i, j, prod[j]);
				}
			}

		}
	}

void time_module_sample_G(int n)
	{
	signed_scalar z_coeffs[PARAM_D * PARAM_K * PARAM_N];
	signed_poly_matrix z = z_coeffs;
	
	printf("Presampling from the module gadget lattice... (not really)");
	fflush(stdout);
	
	clock_t t0 = clock();
	for(int i = 0 ; i < n ; ++i)
		{
		// Sample z such that Gz = 0 mod q
		//module_sample_G(z);
		}
	clock_t t1 = clock() - t0;
	
	printf("Done : %ld\n", t1);
	}


/*
	Checks whether the l by l triangular matrix M has only positive coefficients on its diagonal
*/
bool is_diagonally_positive(cplx_poly_matrix M, int l)
	{
	bool diag_pos = true;
	
	for(int i = 0 ; i < l ; ++i)
		{
		cplx_poly M_ii = triangular_poly_matrix_element(M, i, i);
		
		for(int j = 0 ; j < PARAM_N ; ++j)
			{
			if((fabs(cimag(M_ii[j])) != 0.0) || (creal(M_ii[j]) <= 0.0))
				{
				printf("M[%d,%d][%d] = %f %+f * I\n", i, i, j, creal(M_ii[j]), cimag(M_ii[j]));
				diag_pos = false;
				}
			}
		}
	
	return diag_pos;
	}

void test_sch_comp_computations(void)
	{
	scalar *T_coeffs = malloc(PARAM_N * 2 * PARAM_D * PARAM_D * PARAM_K * sizeof(scalar));
	cplx *sch_comp_coeffs = malloc(PARAM_N * PARAM_D * (2 * PARAM_D + 1) * sizeof(cplx)), *cplx_T_coeffs = malloc(PARAM_N * 2 * PARAM_D * PARAM_D * PARAM_K * sizeof(cplx));
	
	poly_matrix T = T_coeffs;
	cplx_poly_matrix sch_comp = sch_comp_coeffs, cplx_T = cplx_T_coeffs;
	
	// Generate T
	printf("Generating gaussian T...\n");
	
	SampleR_matrix_centered((signed_poly_matrix) T, 2*PARAM_D, PARAM_D * PARAM_K, PARAM_SIGMA);
	
	printf("Done\n");
	
	// Compute cplx_T (CRT domain)
	printf("Computing cplx_T in the CRT domain...\n");
	
	for(int i = 0 ; i < PARAM_N * 2 * PARAM_D * PARAM_D * PARAM_K ; ++i)
		{
		cplx_T[i] = (signed_scalar) T[i];
		}
	
	matrix_cplx_crt_representation(cplx_T, 2 * PARAM_D, PARAM_D * PARAM_K);
	
	printf("Done\n");
	
	// Compute the Schur complements and check if they are diagonally positive
	printf("Computing the Schur complements...\n");
	
	construct_T_schur_complement(sch_comp, cplx_T);
	
	printf("first sch_comp done\n");
	
	if(!is_diagonally_positive(sch_comp, 2 * PARAM_D))
		{
		print_triangular_cplx_poly_matrix(sch_comp, 2 * PARAM_D);
		return;
		}
	
	for(int i = 2 * PARAM_D - 1 ; i > 1 ; --i)
		{
		construct_schur_complement(sch_comp, i);
		
		printf("sch_comp %d done\n", i);
		
		if(!is_diagonally_positive(sch_comp, i))
			{
			print_triangular_cplx_poly_matrix(sch_comp, i);
			return;
			}
		}
	}

int main(int argc, char **argv)
	{
	if(argc < 2)
		{
		fprintf(stderr, "Usage : %s n\n", argv[0]);
		exit(EXIT_FAILURE);
		}
	
	int n = atoi(argv[1]);
	
	init_crt_trees();
	init_D_lattice_coeffs();
	init_cplx_roots_of_unity();
	
	random_bytes_init();
	
	time_sample_pre_subroutines(n);
	
	return 0;
	}
