#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include <math.h>
#include <assert.h>
#include <string.h>

#include "common.h"
#include "arithmetic.h"
#include "sampling.h"

/*
	Generates a signing key (T, cplx_T, sch_comp) and an associated verification key A
*/
void KeyGen(poly_matrix A, poly_matrix T, cplx_poly_matrix cplx_T, cplx_poly_matrix sch_comp)
	{
	//scalar A_hat_coeffs[PARAM_D * PARAM_D * PARAM_N], AprimeT_coeffs[PARAM_D * PARAM_D * PARAM_K * PARAM_N];
	scalar *A_hat_coeffs = malloc(PARAM_D * PARAM_D * PARAM_N * sizeof(scalar)), *AprimeT_coeffs = malloc(PARAM_D * PARAM_D * PARAM_K * PARAM_N * sizeof(scalar));
	poly_matrix A_hat = A_hat_coeffs, AprimeT = AprimeT_coeffs;
	
	// A_hat <- U(R_q^{d,d}) is considered to be in the CRT domain
	random_poly(A_hat, PARAM_N * PARAM_D * PARAM_D - 1);
		
	// T <- D_{R^{2d,dk},sigma}
	SampleR_matrix_centered((signed_poly_matrix) T, 2*PARAM_D, PARAM_D * PARAM_K, PARAM_SIGMA);

	// Compute the Schur complements
	construct_complex_private_key(cplx_T, sch_comp, T);
	
	// Add q to each component of T and put it in the CRT domain
	for(int i = 0 ; i < PARAM_N * 2 * PARAM_D * PARAM_D * PARAM_K ; ++i)
		{
		T[i] += PARAM_Q;
		}

	matrix_crt_representation(T, 2*PARAM_D, PARAM_D * PARAM_K, LOG_R);

	// AprimeT = A_hat * T2 + T1, where T1 and T2 are the upper and lower half of T
	poly_matrix T1 = T, T2 = poly_matrix_element(T, PARAM_D * PARAM_K, PARAM_D, 0);
	
	mul_crt_poly_matrix(AprimeT, A_hat, T2, PARAM_D, PARAM_D, PARAM_D * PARAM_K, LOG_R);
	add_to_poly_matrix(AprimeT, T1, PARAM_D, PARAM_D * PARAM_K);
	
	// A = (A_hat | -A'T) ( = (I | A_hat | -A'T) implicitly)
	for(int i = 0 ; i < PARAM_D ; ++i)
		{
		poly_matrix A_i0 = poly_matrix_element(A, PARAM_M - PARAM_D, i, 0);
		poly_matrix A_hat_i = poly_matrix_element(A_hat, PARAM_D, i, 0);
		
		memcpy(A_i0, A_hat_i, PARAM_D * PARAM_N * sizeof(scalar));
		}
	for(int i = 0 ; i < PARAM_D ; ++i)
		{
		poly_matrix A_i1 = poly_matrix_element(A, PARAM_M - PARAM_D, i, PARAM_D);
		poly_matrix AprimeT_i = poly_matrix_element(AprimeT, PARAM_D * PARAM_K, i, 0);
		for(int j = 0 ; j < PARAM_D * PARAM_K * PARAM_N ; ++j)
			{
			A_i1[j] = 2*PARAM_Q - AprimeT_i[j];
			}
		}
	
	// Reduce A's coefficients mod q
	freeze_poly(A, PARAM_N * PARAM_D * (PARAM_M - PARAM_D) - 1);
	
	free(A_hat);
	free(AprimeT);
	}

/*
	Signs a message m using the signing key (T, cplx_T, sch_comp) and the verification key A
*/
void Sign(poly_matrix nu, poly_matrix A, poly_matrix T, cplx_poly_matrix cplx_T, cplx_poly_matrix sch_comp, scalar *m)
	{
	// Compute m's inverse and put it in the CRT domain
	scalar m_inv_coeffs[PARAM_N];
	poly m_inv = m_inv_coeffs;
	
	invert_poly(m_inv, m, PARAM_N, 1);
	crt_representation(m_inv, LOG_R);
	
	// Use m to construct A_m
	construct_A_m(A, m);
	
	// Sample nu
	sample_pre(nu, A, T, cplx_T, sch_comp, m_inv);
	
	// Deconstruct A_m
	deconstruct_A_m(A, m);
	}

/*
	Checks is the signature nu is valid for the message m, given the verification key A
*/
bool Verify(poly_matrix nu, poly_matrix A, scalar *m)
	{
	// Verify that A * nu = m mod q
	scalar prod_coeffs[PARAM_N * PARAM_D];
	poly_matrix prod = prod_coeffs;
	
	construct_A_m(A, m);
	multiply_by_A(prod, A, (poly_matrix) nu);
	deconstruct_A_m(A, m);
	
	
	if(!is_zero_poly(prod, PARAM_N * PARAM_D - 1))
		{
		
		printf("prod (CRT)\n");
		print_poly_matrix(prod, PARAM_D, 1);
		
		return false;
		}
	
	
	// Verify that nu has a small norm
	matrix_coeffs_representation(nu, PARAM_M, 1, LOG_R);
	
	double_scalar norm_nu_squared = norm_squared((poly_matrix) nu, PARAM_M);
	double bound_squared = PARAM_T * PARAM_T * PARAM_ZETA * PARAM_ZETA * PARAM_N * PARAM_M;

	matrix_crt_representation(nu, PARAM_M, 1, LOG_R);
	
	if(norm_nu_squared >= bound_squared)
		{
		printf("norm^2 = %" PRIu64 " %c bound^2 = %f\n", norm_nu_squared, (norm_nu_squared < bound_squared)?'<':'>', bound_squared);
		
		return false;
		}
	
	return true;
	}
