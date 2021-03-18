#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include <math.h>
#include <assert.h>
#include <string.h>

#include "common.h"
#include "arithmetic.h"
#include "sampling.h"
#include "random.h"
#include "hash.h"

#include "cpucycles.h"


extern unsigned long long timing_sampleZ_KG;
extern unsigned long long timing_precomp_KG;
extern unsigned long long timing_arith_KG;

/*
	Generates a signing key (T, cplx_T, sch_comp) and an associated verification key A
*/
void KeyGen(poly_matrix A, poly_matrix T, cplx_poly_matrix cplx_T, cplx_poly_matrix sch_comp)
	{
	unsigned long long begin_timing = 0;
	unsigned long long end_timing   = 0;
	
	//scalar A_hat_coeffs[PARAM_D * PARAM_D * PARAM_N], AprimeT_coeffs[PARAM_D * PARAM_D * PARAM_K * PARAM_N];
	scalar *A_hat_coeffs = malloc(PARAM_D * PARAM_D * PARAM_N * sizeof(scalar)), *AprimeT_coeffs = malloc(PARAM_D * PARAM_D * PARAM_K * PARAM_N * sizeof(scalar));
	poly_matrix A_hat = A_hat_coeffs, AprimeT = AprimeT_coeffs;
	
	// A_hat <- U(R_q^{d,d}) is considered to be in the NTT_op domain
	random_poly(A_hat, PARAM_N * PARAM_D * PARAM_D - 1);
	
	
	// T <- D_{R^{2d,dk},sigma}
	begin_timing = cpucycles_start();
	
	SampleR_matrix_centered((signed_poly_matrix) T, 2*PARAM_D, PARAM_D * PARAM_K, PARAM_SIGMA);
	
	end_timing = cpucycles_stop();
	timing_sampleZ_KG += (end_timing - begin_timing);

	// Compute the Schur complements
	begin_timing = cpucycles_start();
	
	construct_complex_private_key(cplx_T, sch_comp, T);
	
	end_timing = cpucycles_stop();
	timing_precomp_KG += (end_timing - begin_timing);
	
	
	begin_timing = cpucycles_start();
	
	// Add q to each component of T and put it in the CRT domain
	for(int i = 0 ; i < PARAM_N * 2 * PARAM_D * PARAM_D * PARAM_K ; ++i)
		{
		T[i] += PARAM_Q;
		}

	matrix_ntt(T, 2*PARAM_D, PARAM_D * PARAM_K);
	freeze_poly(T, PARAM_N * 2 * PARAM_D * PARAM_D * PARAM_K - 1);

	// T1 and T2 are the upper and lower half of T
	poly_matrix T1 = T, T2 = poly_matrix_element(T, PARAM_D * PARAM_K, PARAM_D, 0);
	
	// AprimeT <- A_hat * T2
	mul_crt_poly_matrix(AprimeT, A_hat, T2, PARAM_D, PARAM_D, PARAM_D * PARAM_K);
	
	// Convert AprimeT from the NTT_res domain into the NTT_op domain
	multiply_by_2pow32(AprimeT, PARAM_N * PARAM_D * PARAM_D * PARAM_K - 1);
	
	// AprimeT <- AprimeT + T1
	add_to_poly_matrix(AprimeT, T1, PARAM_D, PARAM_D * PARAM_K);
	
	freeze_poly(AprimeT, PARAM_N * PARAM_D * PARAM_D * PARAM_K - 1);
	
	
	// AprimeT <- - AprimeT
	for(int i = 0 ; i < PARAM_N * PARAM_D * PARAM_D * PARAM_K ; ++i)
		{
		AprimeT[i] = 2 * PARAM_Q - AprimeT[i];
		}
	
	
	// AprimeT <- AprimeT + G
	for(int i = 0 ; i < PARAM_D ; ++i)
		{
		poly_matrix AprimeT_iik = poly_matrix_element(AprimeT, PARAM_D * PARAM_K, i, i * PARAM_K);
		
		add_ring_gadget_vector(AprimeT_iik);
		}
	
	
	// A = (A_hat | -A'T) ( = (I | A_hat | -A'T) implicitly)
	for(int i = 0 ; i < PARAM_D ; ++i)
		{
		poly_matrix A_i0 = poly_matrix_element(A, PARAM_M - PARAM_D, i, 0);
		poly_matrix A_hat_i = poly_matrix_element(A_hat, PARAM_D, i, 0);
		
		memcpy(A_i0, A_hat_i, PARAM_D * PARAM_N * sizeof(scalar));
		}
	for(int i = 0 ; i < PARAM_D ; ++i)
		{
		poly_matrix A_id = poly_matrix_element(A, PARAM_M - PARAM_D, i, PARAM_D);
		poly_matrix AprimeT_i = poly_matrix_element(AprimeT, PARAM_D * PARAM_K, i, 0);
		
		memcpy(A_id, AprimeT_i, PARAM_N * PARAM_D * PARAM_K * sizeof(scalar));
		}
	
	// Reduce A's coefficients mod q
	freeze_poly(A, PARAM_N * PARAM_D * (PARAM_M - PARAM_D) - 1);
	
	
	end_timing = cpucycles_stop();
	timing_arith_KG += (end_timing - begin_timing);
	free(A_hat);
	free(AprimeT);
	}

/*
	Signs a message m of length m_len using the signing key (T, cplx_T, sch_comp) and the verification key A
*/
void Sign(poly_matrix nu, uint8_t *r, poly_matrix A, poly_matrix T, cplx_poly_matrix cplx_T, cplx_poly_matrix sch_comp, uint8_t *m, int m_len)
	{
	// Generate a random salt
	salt(r);
	
	// Construct a target u using the message and the salt
	scalar u_coeffs[PARAM_D * PARAM_N];
	poly_matrix u = u_coeffs;
	
	H(u, m, m_len, r);

	
	// Sample nu
	sample_pre(nu, A, T, cplx_T, sch_comp, u);
	}

/*
	Checks is the signature nu is valid for the message m, given the verification key A
*/
bool Verify(poly_matrix nu, uint8_t *r, poly_matrix A, uint8_t *m, int m_len)
	{
	// Construct a target u using the message and the salt
	scalar u_coeffs[PARAM_D * PARAM_N];
	poly_matrix u = u_coeffs;
	
	H(u, m, m_len, r);
	
	
	// Verify that A * nu = u mod q
	scalar prod_coeffs[PARAM_N * PARAM_D];
	poly_matrix prod = prod_coeffs;
	
	multiply_by_A(prod, A, (poly_matrix) nu);
	
	
	
	// Verify that nu has a small norm
	divide_by_2pow32(nu, PARAM_N * PARAM_M - 1);
	matrix_invntt(nu, PARAM_M, 1);
	freeze_poly(nu, PARAM_N * PARAM_M - 1);
	
	double_scalar norm_nu_squared = norm_squared((poly_matrix) nu, PARAM_M);
	double bound_squared = PARAM_T * PARAM_T * PARAM_ZETA * PARAM_ZETA * PARAM_N * PARAM_M;
	
	
	matrix_ntt(nu, PARAM_M, 1);
	
	if(!equals_poly(prod, u, PARAM_N * PARAM_D - 1) || (norm_nu_squared >= bound_squared))
		{
		return false;
		}
	
	return true;
	}
