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
	Generates the master public key (A,u) and the master secret key (T, cplx_T, sch_comp).
*/
void Setup(poly_matrix A, poly_matrix T, cplx_poly_matrix cplx_T, cplx_poly_matrix sch_comp, poly_matrix u)
	{

	// scalar A_hat_coeffs[PARAM_D * PARAM_D * PARAM_N], AprimeT_coeffs[PARAM_D * PARAM_D * PARAM_K * PARAM_N];

	scalar *A_hat_coeffs = malloc(PARAM_D * PARAM_D * PARAM_N * sizeof(scalar)), *AprimeT_coeffs = malloc(PARAM_D * PARAM_D * PARAM_K * PARAM_N * sizeof(scalar));
	poly_matrix A_hat = A_hat_coeffs, AprimeT = AprimeT_coeffs;
	
	// A_hat <- U(R_q^{d,d}) is considered to be in the CRT domain
	random_poly(A_hat, PARAM_N * PARAM_D * PARAM_D - 1);
		
	// T <- D_{R^{2d,dk},sigma}
	SampleR_matrix_centered((signed_poly_matrix) T, 2*PARAM_D, PARAM_D * PARAM_K, PARAM_SIGMA);

	// Compute the Schur complements (that are necessary to the Gaussian preimage sampling operation) + cplx_T
	construct_complex_private_key(cplx_T, sch_comp, T);
	
	// Add q to each component of T (so that T's coeffs are positive) and put it in the CRT domain
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


	random_poly(u, PARAM_N * PARAM_D - 1);

	
	free(A_hat);
	free(AprimeT);
	}






/*
	Generate the secret key of an identity id using the master secret key (T, cplx_T, sch_comp) and the master public key (A,u).
*/
void Extract(poly_matrix x, poly_matrix A, poly_matrix u, poly_matrix T, cplx_poly_matrix cplx_T, cplx_poly_matrix sch_comp, scalar *id)
	{
	// Compute id's inverse and put it in the CRT domain
	scalar id_inv_coeffs[PARAM_N];
	poly id_inv = id_inv_coeffs;
	
	invert_poly(id_inv, id, PARAM_N, 1);
	crt_representation(id_inv, LOG_R);
	
	// Use id to construct A_id
	construct_A_m(A, id);
	
	// Sample x
	// Use of sample_pre_target, which is sample_pre with a target u (needs id_inv as an argument)
	sample_pre_target(x, A, T, cplx_T, sch_comp, id_inv, u);
	
	// Deconstruct A_m
	deconstruct_A_m(A, id);
	}




/*
	Description de la fonction à compléter.
*/
void Encrypt(poly_matrix A, poly_matrix u, scalar *id, poly M, poly_matrix b, poly c)
{

	// Use id to construct A_id
	construct_A_m(A, id);


	scalar *s_coeffs = malloc(PARAM_D * PARAM_N * sizeof(scalar));
	poly_matrix s = s_coeffs;

	signed_scalar *e_0_coeffs = malloc((PARAM_M - PARAM_K) * PARAM_N * sizeof(signed_scalar));
	signed_scalar *e_1_coeffs = malloc(PARAM_K * PARAM_N * sizeof(signed_scalar));
	signed_scalar *e_prime_coeffs = malloc(PARAM_N * sizeof(signed_scalar));

	signed_poly_matrix e_0 = e_0_coeffs, e_1 = e_1_coeffs;
	signed_poly e_prime = e_prime_coeffs;

	// Sampling of s
	random_poly(s, PARAM_N * PARAM_D - 1);

	// Sampling of the errors

	// e_0 <- D_{R^{m-k,1},tau}
	SampleR_matrix_centered(e_0, PARAM_M - PARAM_K, 1, PARAM_TAU);

	// e_1 <- D_{R^{k,1},gamma}
	SampleR_matrix_centered(e_1, PARAM_K, 1, PARAM_GAMMA);

	// e_prime <- D_{R,tau}
	SampleR_centered(e_prime, PARAM_TAU);




	memset(b, 0, PARAM_M * PARAM_N * sizeof(scalar));

	signed_scalar *factor1_coeffs = malloc(PARAM_M * PARAM_N * sizeof(signed_scalar));
	signed_scalar *factor1_t_coeffs = malloc(PARAM_M * PARAM_N * sizeof(signed_scalar));
	signed_scalar *s_t_coeffs = malloc(PARAM_D * PARAM_N * sizeof(scalar));
	signed_scalar *e_0_t_coeffs = malloc((PARAM_M - PARAM_K) * PARAM_N * sizeof(signed_scalar));
	signed_scalar *e_1_t_coeffs = malloc(PARAM_K * PARAM_N * sizeof(signed_scalar));


	signed_poly_matrix factor1 = factor1_coeffs, factor1_t = factor1_t_coeffs, s_t = s_t_coeffs, e_0_t = e_0_t_coeffs, e_1_t = e_1_t_coeffs;

	transpose_signed_scalar_matrix2(s_t, s, PARAM_D, 1);
	transpose_signed_scalar_matrix(e_0_t, e_0, PARAM_M - PARAM_K, 1);
	transpose_signed_scalar_matrix(e_1_t, e_1, PARAM_K, 1);

	matrix_crt_representation((poly_matrix) s_t, 1, PARAM_D, LOG_R);
	mul_crt_poly_matrix((poly_matrix) factor1, (poly_matrix) s_t, A, 1, PARAM_D, PARAM_M, LOG_R);
	transpose_signed_scalar_matrix(factor1_t, factor1, 1, PARAM_M);


	scalar *factor2_coeffs = malloc(PARAM_M * PARAM_N * sizeof(scalar));
	scalar *factor2_t_coeffs = malloc(PARAM_M * PARAM_N * sizeof(scalar));	
	poly_matrix factor2 = factor2_coeffs, factor2_t = factor2_t_coeffs;

	for (int k = 0 ; k < PARAM_M - PARAM_K ; k++)
	{
		poly_matrix factor2_k = poly_matrix_element(factor2, PARAM_M, 0, k);
		poly_matrix e_0_t_k = (poly_matrix) poly_matrix_element(e_0_t, PARAM_M - PARAM_K, 0, k);

		factor2_k = e_0_t_k;
	}


	for (int k = PARAM_M - PARAM_K ; k < PARAM_M ; k++)
	{
		poly_matrix factor2_k = poly_matrix_element(factor2, PARAM_M, 0, k);
		poly_matrix e_1_t_k = (poly_matrix) poly_matrix_element(e_1_t, PARAM_K, 0, k);

		factor2_k = e_1_t_k;
	}

	transpose_signed_scalar_matrix((signed_scalar *) factor2_t, (signed_scalar *) factor2, 1, PARAM_M);

	add_to_poly_matrix(b, (poly_matrix) factor1_t, 1, PARAM_M);
	add_to_poly_matrix(b, factor2_t, 1, PARAM_M);


	scalar *prod_coeffs = malloc(PARAM_N * sizeof(scalar));
	poly_matrix prod = prod_coeffs;
	mul_crt_poly_matrix(prod, (poly_matrix) s_t, u, 1, PARAM_D, 1, LOG_R);

	c = poly_matrix_element(prod, 1,0,0);
	add_poly(c, c, (poly) e_prime, PARAM_N-1);

	for (int i=0 ; i <= PARAM_N ; i++)
		{
			M[i] = floor(PARAM_Q /2) * M[i];
		}

	add_poly(c, c, M, PARAM_N-1);

	free(s);
	free(e_0_t);
	free(e_1_t);
	free(e_prime);
	free(prod);


}



void Decrypt(poly_matrix x, poly_matrix b, poly c, poly M)
{
	scalar *res_coeffs = malloc(PARAM_N * sizeof(scalar));
	poly res = res_coeffs;

	scalar *b_t_coeffs = malloc(PARAM_N * PARAM_M * sizeof(scalar));
	poly_matrix b_t = b_t_coeffs;
	transpose_scalar_matrix(b_t, b, PARAM_M, 1);

	scalar *factor_coeffs = malloc(PARAM_N * sizeof(scalar));
	poly_matrix factor = factor_coeffs;

	matrix_crt_representation(b_t, 1, PARAM_M, LOG_R);

	for (int i = 0 ; i < PARAM_N * PARAM_M ; i++)
		{
		x[i] += PARAM_Q;
		}
	matrix_crt_representation(x, 1, PARAM_M, LOG_R);

	mul_crt_poly_matrix(factor, b_t, x, 1, PARAM_M, 1, LOG_R);

	matrix_coeffs_representation(factor, 1, 1, LOG_R);



	for (int i = 0 ; i < PARAM_N ; i++)
	{
		factor[i] = 2*PARAM_Q + (c[i] - factor[i]);
	}

	freeze_poly(factor, PARAM_N-1);

	res = factor;


	for (int i = 0 ; i < PARAM_N ; i++)
	{
		if (res[i] < floor(PARAM_Q/2)){
			M[i] = 1;
		}

		else {
			M[i] = 0;
		}
	}

}


/*bool Verify_target(poly_matrix nu, poly_matrix A, scalar *m, poly_matrix u)
	{
	// Verify that A * nu = u mod q
	scalar prod_coeffs[PARAM_N * PARAM_D];
	poly_matrix prod = prod_coeffs;
	
	construct_A_m(A, m);
	multiply_by_A(prod, A, (poly_matrix) nu);
	deconstruct_A_m(A, m);

	sub_poly(prod, prod, u, PARAM_N * PARAM_D - 1); // vérif. parameters
	
	
	if(!is_zero_poly(prod, PARAM_N * PARAM_D - 1))
		{
		
		printf("prod (CRT)\n");
		print_poly_matrix(prod, PARAM_D, 1);
		
		return false;
		}
	
	
	// Verify that nu has a small norm
	// matrix_coeffs_representation(nu, PARAM_M, 1, LOG_R);
	
	// double_scalar norm_nu_squared = norm_squared((poly_matrix) nu, PARAM_M);
	// double bound_squared = PARAM_T * PARAM_T * PARAM_ZETA * PARAM_ZETA * PARAM_N * PARAM_M;

	// matrix_crt_representation(nu, PARAM_M, 1, LOG_R);
	
	// if(norm_nu_squared >= bound_squared)
	// 	{
	// 	printf("norm^2 = %" PRIu64 " %c bound^2 = %f\n", norm_nu_squared, (norm_nu_squared < bound_squared)?'<':'>', bound_squared);
		
	// 	return false;
	// 	}
	
	return true;
	}*/