#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include <time.h>
#include <string.h>

#include "common.h"
#include "arithmetic.h"
#include "sampling.h"
#include "signature.h"

void check_trapgen(void)
	{
	scalar A_coeffs[PARAM_D * (PARAM_M - PARAM_D) * PARAM_N], TI_coeffs[PARAM_M * PARAM_D * PARAM_K * PARAM_N], prod_coeffs[PARAM_D * PARAM_D * PARAM_K * PARAM_N], AwithI_coeffs[PARAM_D * PARAM_M * PARAM_N];
	poly_matrix A = A_coeffs, T = TI_coeffs, TI = TI_coeffs, Idk = poly_matrix_element(TI, PARAM_D * PARAM_K, 2*PARAM_D, 0), prod = prod_coeffs, AwithI = AwithI_coeffs;
	
	// Generate A and T
	TrapGen(A, T);
	
	// Generate Idk (in the CRT domain) and append it to T vertically
	memset(Idk, 0, PARAM_D * PARAM_K * PARAM_D * PARAM_K * PARAM_N * sizeof(scalar));
	for(int i = 0 ; i < PARAM_D * PARAM_K ; ++i)
		{
		poly p_i = poly_matrix_element(Idk, PARAM_D * PARAM_K, i, i);
		for(int j = 0 ; j < PARAM_R ; ++j)
			{
			poly p_ij = crt_poly_component(p_i, SMALL_DEGREE, j);
			p_ij[0] = 1;
			}
		}
	
	// Prepend Id (in the CRT domain) to A horizontally
	for(int i = 0 ; i < PARAM_D ; ++i)
		{
		poly_matrix AwithI_i = poly_matrix_element(AwithI, PARAM_M, i, 0);
		poly_matrix A_i = poly_matrix_element(A, PARAM_M - PARAM_D, i, 0);
		
		// build Id
		memset(AwithI_i, 0, PARAM_D * PARAM_N * sizeof(scalar));
		poly p_i = poly_matrix_element(AwithI, PARAM_M, i, i);
		for(int j = 0 ; j < PARAM_R ; ++j)
			{
			poly p_ij = crt_poly_component(p_i, SMALL_DEGREE, j);
			p_ij[0] = 1;
			}
		
		// copy A
		memcpy(poly_matrix_element(AwithI_i, PARAM_M, 0, PARAM_D), A_i, (PARAM_M - PARAM_D) * PARAM_N * sizeof(scalar));
		}
	
	
	// prod = AwithI * TI
	mul_crt_poly_matrix(prod, AwithI, TI, PARAM_D, PARAM_M, PARAM_D * PARAM_K, LOG_R);
	
	
	for(int i = 0 ; i < PARAM_D * PARAM_D * PARAM_K * PARAM_N ; ++i)
		{
		if(prod[i] != 0)
			{
			printf("%d : %" PRIu32 "\n", i, prod[i]);
			}
		else
			{
			//printf("yay\n");
			}
		}
	}

int main(int argc, char **argv)
	{
	init_crt_trees();
	init_D_lattice_coeffs();
	init_cplx_roots_of_unity();
	
	srand48(0);
	srand(0);
	
	check_trapgen();
	
	return 0;
	}
